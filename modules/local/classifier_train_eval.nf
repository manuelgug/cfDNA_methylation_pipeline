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
    path "classifier_metrics.json",  emit: metrics_json
    path "roc_data.json",            emit: roc_data
    path "per_sample_scores.csv",    emit: sample_scores
    path "feature_importance.csv",   emit: feature_importance
    path "*.pdf",                    emit: plots

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
