/*
 * Module: FEATURE_EXTRACTION
 */
process FEATURE_EXTRACTION {
    tag "all_samples"
    label 'process_highmem'

    publishDir "${params.outdir}/features", mode: params.publish_mode

    input:
    path cov_files
    path samplesheet

    output:
    path "beta_matrix_variable.csv.gz",  emit: beta_matrix
    path "sample_metadata_qc.csv",       emit: sample_metadata
    path "sample_qc_metrics.csv",        emit: qc_summary
    path "feature_extraction_summary.json"

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

/*
 * Module: DIFFERENTIAL_METHYLATION
 */
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
    path "DMPs.csv.gz",          emit: dmp_table
    path "DMRs*.csv.gz",         emit: dmr_features
    path "dm_summary.json",      emit: summary_json
    path "results_table.csv.gz", emit: results_table
    path "*.pdf"

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

    # Create combined results table for report
    cp DMPs.csv.gz results_table.csv.gz
    """
}

/*
 * Module: CLASSIFIER_TRAIN_EVAL
 */
process CLASSIFIER_TRAIN_EVAL {
    tag "ML_${method}"
    label 'process_high'

    publishDir "${params.outdir}/classifier", mode: params.publish_mode

    input:
    path beta_matrix
    path dmr_features
    path metadata
    val  method
    val  n_features
    val  cv_folds
    val  test_split
    val  seed

    output:
    path "classifier_metrics.json",      emit: metrics_json
    path "roc_data.json",                emit: roc_data
    path "per_sample_scores.csv",        emit: sample_scores
    path "feature_importance.csv",       emit: feature_importance
    path "*.pdf"

    script:
    """
    python3 ${projectDir}/bin/classifier_train_eval.py \\
        --beta_matrix  ${beta_matrix} \\
        --dmr_features ${dmr_features} \\
        --metadata     ${metadata} \\
        --method       ${method} \\
        --n_features   ${n_features} \\
        --cv_folds     ${cv_folds} \\
        --test_split   ${test_split} \\
        --seed         ${seed} \\
        --outdir       .
    """
}

/*
 * Module: REPORT_GENERATE
 */
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
    path "udx_pipeline_report.html", emit: html_report
    path "udx_pipeline_summary.json"

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

// ── Stub modules (shims for Nextflow module resolution) ────────────────────────
process TRIM_GALORE {
    tag "${meta.id}"; label 'process_medium'
    publishDir "${params.outdir}/trimmed/${meta.id}", mode: params.publish_mode
    input:  tuple val(meta), path(reads)
    output: tuple val(meta), path("*_trimmed.fq.gz"), emit: reads
            tuple val(meta), path("*trimming_report.txt"), emit: log
    script: "trim_galore --cores ${task.cpus} --gzip --paired ${reads}"
}

process BISMARK_ALIGN {
    tag "${meta.id}"; label 'process_long'
    publishDir "${params.outdir}/aligned/${meta.id}", mode: params.publish_mode
    input:  tuple val(meta), path(reads); path bismark_index
    output: tuple val(meta), path("*.bam"), emit: bam
            tuple val(meta), path("*_report.txt"), emit: report
    script: "bismark --genome ${bismark_index} --bowtie2 --parallel ${task.cpus.intdiv(4)} -1 ${reads[0]} -2 ${reads[1]}"
}

process PICARD_MARKDUPLICATES {
    tag "${meta.id}"; label 'process_medium'
    publishDir "${params.outdir}/aligned/${meta.id}", mode: params.publish_mode
    input:  tuple val(meta), path(bam)
    output: tuple val(meta), path("*.markdup.bam"), emit: bam
            tuple val(meta), path("*.markdup.bam.bai"), emit: bai
            tuple val(meta), path("*.metrics.txt"), emit: metrics
    script: "picard MarkDuplicates I=${bam} O=${meta.id}.markdup.bam M=${meta.id}.markdup.metrics.txt REMOVE_DUPLICATES=false && samtools index ${meta.id}.markdup.bam"
}

process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"; label 'process_low'
    publishDir "${params.outdir}/qc/flagstat/${meta.id}", mode: params.publish_mode
    input:  tuple val(meta), path(bam)
    output: tuple val(meta), path("*.flagstat"), emit: flagstat
    script: "samtools flagstat --threads ${task.cpus} ${bam} > ${meta.id}.flagstat"
}

process BISMARK_METHYLATION_EXTRACTOR {
    tag "${meta.id}"; label 'process_high'
    publishDir "${params.outdir}/methylation/${meta.id}", mode: params.publish_mode
    input:  tuple val(meta), path(bam); path bismark_index
    output: tuple val(meta), path("*.bismark.cov.gz"), emit: bismark_cov
            tuple val(meta), path("*.CX_report.txt.gz"), emit: cx_report
            tuple val(meta), path("*_splitting_report.txt"), emit: report
            tuple val(meta), path("*.M-bias.txt"), emit: mbias
    script: "bismark_methylation_extractor --comprehensive --cytosine_report --genome_folder ${bismark_index} --parallel ${task.cpus} --bedGraph --gzip --paired-end ${bam}"
}

process COVERAGE_FILTER {
    tag "${meta.id}"; label 'process_low'
    publishDir "${params.outdir}/methylation_filtered/${meta.id}", mode: params.publish_mode
    input:  tuple val(meta), path(cov_gz); val min_cov
    output: tuple val(meta), path("*.filtered.cov.gz"), emit: filtered_cov
    script: "zcat ${cov_gz} | awk -v mc=${min_cov} '(\$5+\$6) >= mc' | gzip > ${meta.id}.filtered.cov.gz"
}

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: params.publish_mode
    input:  path qc_files
    output: path "multiqc_report.html", emit: report
            path "multiqc_data/", emit: data
    script: "multiqc . --force --title 'Universal DX cfDNA Pipeline QC' -o ."
}
