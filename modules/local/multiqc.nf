process MULTIQC {
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: params.publish_mode

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . \\
        --force \\
        --title "Universal DX cfDNA Pipeline QC" \\
        -o .
    """
}
