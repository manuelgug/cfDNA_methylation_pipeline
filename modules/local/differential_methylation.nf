process DIFFERENTIAL_METHYLATION {
    tag "DMA_${method}"
    label 'process_high'

    publishDir "${params.outdir}/differential_methylation", mode: params.publish_mode

    input:
    path beta_matrix
    path metadata
    val  method
    val  delta_beta
    val  fdr_cutoff
    val  min_cpgs_dmr

    output:
    path "DMPs.csv.gz",           emit: dmp_table
    path "DMRs*.csv.gz",          emit: dmr_features
    path "dm_summary.json",       emit: summary_json
    path "results_table.csv.gz",  emit: results_table
    path "*.pdf",                 emit: plots

    script:
    """
    Rscript ${projectDir}/bin/differential_methylation.R \\
        --beta_matrix  ${beta_matrix} \\
        --metadata     ${metadata} \\
        --method       ${method} \\
        --delta_beta   ${delta_beta} \\
        --fdr          ${fdr_cutoff} \\
        --min_cpgs_dmr ${min_cpgs_dmr} \\
        --outdir       .

    cp DMPs.csv.gz results_table.csv.gz
    """
}
