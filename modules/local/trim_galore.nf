process TRIM_GALORE {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/trimmed/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fq.gz"),      emit: reads
    tuple val(meta), path("*trimming_report.txt"),  emit: log

    script:
    def paired       = reads instanceof List && reads.size() > 1
    def paired_flag  = paired ? "--paired" : ""
    def nextseq_flag = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ""
    """
    trim_galore \\
        --cores ${task.cpus} \\
        --quality 20 \\
        --phred33 \\
        --length 20 \\
        --gzip \\
        ${nextseq_flag} \\
        ${paired_flag} \\
        ${reads}
    """
}
