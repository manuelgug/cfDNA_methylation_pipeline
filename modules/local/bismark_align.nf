process BISMARK_ALIGN {
    tag "${meta.id}"
    label 'process_long'

    publishDir "${params.outdir}/aligned/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    path  bismark_index

    output:
    tuple val(meta), path("*.bam"),                     emit: bam
    tuple val(meta), path("*_bismark_bt2_report.txt"),  emit: report

    script:
    def paired    = reads instanceof List && reads.size() > 1
    def reads_arg = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "${reads}"
    def mode_flag = paired ? "--no_mixed --no_discordant" : ""
    def par       = Math.max(1, task.cpus.intdiv(4))
    """
    bismark \\
        --genome ${bismark_index} \\
        --bowtie2 \\
        --parallel ${par} \\
        --non_directional \\
        --nucleotide_coverage \\
        ${mode_flag} \\
        ${reads_arg}
    """
}
