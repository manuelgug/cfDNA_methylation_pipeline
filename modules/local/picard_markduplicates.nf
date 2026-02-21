process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/aligned/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.markdup.bam"),      emit: bam
    tuple val(meta), path("*.markdup.bam.bai"),  emit: bai
    tuple val(meta), path("*.metrics.txt"),      emit: metrics

    script:
    def mem_gb = task.memory ? task.memory.toGiga() - 1 : 4
    """
    picard -Xmx${mem_gb}g MarkDuplicates \\
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
