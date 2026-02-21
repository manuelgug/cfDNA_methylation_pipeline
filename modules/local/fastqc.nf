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
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}
    """
}
