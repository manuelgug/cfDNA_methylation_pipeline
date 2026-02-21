/*
 * Module: FASTQC
 */
process FASTQC {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/qc/fastqc/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    def prefix  = meta.id
    def paired  = reads instanceof List ? "--paired" : ""
    def threads = task.cpus
    """
    fastqc \\
        --threads ${threads} \\
        --outdir . \\
        ${reads}
    """
}

/*
 * Module: TRIM_GALORE
 * Bisulfite-aware: retains methylation info while removing adapters.
 */
process TRIM_GALORE {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/trimmed/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fq.gz"),    emit: reads
    tuple val(meta), path("*trimming_report.txt"), emit: log

    script:
    def paired       = reads instanceof List
    def paired_flag  = paired ? "--paired" : ""
    def nextseq_flag = params.trim_nextseq > 0
                       ? "--nextseq ${params.trim_nextseq}" : ""
    """
    trim_galore \\
        --cores ${task.cpus} \\
        --quality 20 \\
        --phred33 \\
        --stringency 1 \\
        --length 20 \\
        --dont_gzip \\
        --gzip \\
        ${nextseq_flag} \\
        ${paired_flag} \\
        ${reads}
    """
}

/*
 * Module: BISMARK_ALIGN
 * Two-pass bisulfite alignment optimised for cfDNA (short fragments).
 */
process BISMARK_ALIGN {
    tag "${meta.id}"
    label 'process_long'

    publishDir "${params.outdir}/aligned/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    path  bismark_index

    output:
    tuple val(meta), path("*.bam"),                    emit: bam
    tuple val(meta), path("*_bismark_bt2_report.txt"), emit: report

    script:
    def paired      = reads instanceof List
    def reads_arg   = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : reads
    def mode_flag   = paired ? "--no_mixed --no_discordant" : ""
    def pbat_flag   = meta.library == 'WGBS' ? "--pbat" : ""
    """
    bismark \\
        --genome ${bismark_index} \\
        --bowtie2 \\
        --parallel ${Math.max(1, task.cpus.intdiv(4))} \\
        --score_min L,0,-0.6 \\
        --minins 0 \\
        --maxins 1000 \\
        --non_directional \\
        --nucleotide_coverage \\
        ${mode_flag} \\
        ${pbat_flag} \\
        ${reads_arg}
    """
}

/*
 * Module: PICARD_MARKDUPLICATES
 * For cfDNA: marks (does NOT remove) PCR duplicates — important for UMI-based
 * deduplication context and fragment-level analysis.
 */
process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/aligned/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.markdup.bam"),     emit: bam
    tuple val(meta), path("*.markdup.bam.bai"), emit: bai
    tuple val(meta), path("*.metrics.txt"),     emit: metrics

    script:
    def mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    """
    picard ${mem} MarkDuplicates \\
        I=${bam} \\
        O=${meta.id}.markdup.bam \\
        M=${meta.id}.markdup.metrics.txt \\
        VALIDATION_STRINGENCY=LENIENT \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=false \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500

    samtools index ${meta.id}.markdup.bam
    """
}

/*
 * Module: SAMTOOLS_FLAGSTAT
 */
process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/qc/flagstat/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat

    script:
    """
    samtools flagstat --threads ${task.cpus} ${bam} > ${meta.id}.flagstat
    """
}

/*
 * Module: BISMARK_METHYLATION_EXTRACTOR
 * Generates bismark coverage (CX_report) and bedGraph outputs.
 */
process BISMARK_METHYLATION_EXTRACTOR {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/methylation/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path  bismark_index

    output:
    tuple val(meta), path("*.bismark.cov.gz"),  emit: bismark_cov
    tuple val(meta), path("*.CX_report.txt.gz"), emit: cx_report
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt"),       emit: mbias

    script:
    def paired    = bam =~ /pe/ ? "--paired-end" : "--single-end"
    def cx_flag   = params.chg_context || params.chh_context ? "--CX_context" : ""
    """
    bismark_methylation_extractor \\
        --comprehensive \\
        --cytosine_report \\
        --genome_folder ${bismark_index} \\
        --parallel ${task.cpus} \\
        --buffer_size ${task.memory.toGiga() - 2}G \\
        --bedGraph \\
        --gzip \\
        ${cx_flag} \\
        ${paired} \\
        ${bam}
    """
}

/*
 * Module: COVERAGE_FILTER
 * Removes CpGs below minimum coverage threshold; outputs BED-style matrix.
 */
process COVERAGE_FILTER {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/methylation_filtered/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(cov_gz)
    val   min_cov

    output:
    tuple val(meta), path("*.filtered.cov.gz"), emit: filtered_cov

    script:
    """
    zcat ${cov_gz} \\
        | awk -v mc=${min_cov} '(\$5+\$6) >= mc' \\
        | gzip > ${meta.id}.filtered.cov.gz
    """
}

/*
 * Module: MULTIQC
 */
process MULTIQC {
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: params.publish_mode

    input:
    path qc_files

    output:
    path "multiqc_report.html",     emit: report
    path "multiqc_data/",           emit: data

    script:
    """
    multiqc . \\
        --force \\
        --title "Universal DX cfDNA Pipeline QC" \\
        --comment "Generated by UDX methylation pipeline v1.0" \\
        -o .
    """
}
