process BISMARK_METHYLATION_EXTRACTOR {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/methylation/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path  bismark_index

    output:
    tuple val(meta), path("*.bismark.cov.gz"),       emit: bismark_cov
    tuple val(meta), path("*.CX_report.txt.gz"),     emit: cx_report
    tuple val(meta), path("*_splitting_report.txt"), emit: report
    tuple val(meta), path("*.M-bias.txt"),           emit: mbias

    script:
    def paired   = bam.name =~ /_pe\./ ? "--paired-end" : "--single-end"
    def cx_flag  = (params.chg_context || params.chh_context) ? "--CX_context" : ""
    def buf_gb   = task.memory ? task.memory.toGiga() - 2 : 6
    """
    bismark_methylation_extractor \\
        --comprehensive \\
        --cytosine_report \\
        --genome_folder ${bismark_index} \\
        --parallel ${task.cpus} \\
        --buffer_size ${buf_gb}G \\
        --bedGraph \\
        --gzip \\
        ${cx_flag} \\
        ${paired} \\
        ${bam}
    """
}
