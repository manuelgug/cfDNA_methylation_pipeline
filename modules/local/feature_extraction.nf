process FEATURE_EXTRACTION {
    tag "all_samples"
    label 'process_highmem'

    publishDir "${params.outdir}/features", mode: params.publish_mode

    input:
    path cov_files
    path samplesheet

    output:
    path "beta_matrix_variable.csv.gz",       emit: beta_matrix
    path "sample_metadata_qc.csv",            emit: sample_metadata
    path "sample_qc_metrics.csv",             emit: qc_summary
    path "feature_extraction_summary.json",   emit: summary

    script:
    """
    python3 ${projectDir}/bin/feature_extraction.py \\
        --cov_files   ${cov_files} \\
        --samplesheet ${samplesheet} \\
        --min_cov     ${params.min_coverage} \\
        --top_variable 100000 \\
        --output_dir  .
    """
}
