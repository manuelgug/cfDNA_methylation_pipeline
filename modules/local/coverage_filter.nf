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
        | awk -v mc=${min_cov} '(\$5 + \$6) >= mc' \\
        | gzip > ${meta.id}.filtered.cov.gz
    """
}
