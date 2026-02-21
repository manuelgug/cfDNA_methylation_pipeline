process REPORT_GENERATE {
    tag "report"
    label 'process_low'

    publishDir "${params.outdir}/report", mode: params.publish_mode

    input:
    path multiqc_report
    path qc_summary
    path dm_results
    path metrics_json
    path roc_data

    output:
    path "udx_pipeline_report.html",  emit: html_report
    path "udx_pipeline_summary.json", emit: summary_json

    script:
    """
    python3 ${projectDir}/bin/generate_report.py \\
        --multiqc_report ${multiqc_report} \\
        --qc_summary     ${qc_summary} \\
        --dm_results     ${dm_results} \\
        --metrics_json   ${metrics_json} \\
        --roc_data       ${roc_data} \\
        --outdir         .
    """
}
