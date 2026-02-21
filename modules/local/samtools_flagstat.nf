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
