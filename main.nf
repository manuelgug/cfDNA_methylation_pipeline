#!/usr/bin/env nextflow

/*
================================================================================
  Universal DX  |  cfDNA Methylation & NGS Analytical Pipeline  v1.0.0
================================================================================
*/

nextflow.enable.dsl = 2

// ── Parameter defaults (must use params.key = val syntax in DSL2 scripts) ────
params.input         = null
params.genome        = null
params.bismark_index = null

params.min_coverage  = 10
params.min_mapq      = 20
params.trim_nextseq  = 20
params.cpg_context   = true
params.chg_context   = false
params.chh_context   = false

params.dm_method     = "methylkit"
params.delta_beta    = 0.20
params.fdr_cutoff    = 0.05
params.min_cpgs_dmr  = 3

params.ml_method     = "rf"
params.n_features    = 500
params.cv_folds      = 5
params.test_split    = 0.20
params.random_seed   = 42

params.outdir        = "results"
params.publish_mode  = "copy"

params.max_memory    = "64.GB"
params.max_cpus      = 16
params.max_time      = "72.h"

// ── Include modules ───────────────────────────────────────────────────────────
include { VALIDATE_SAMPLESHEET          } from './modules/local/validate_samplesheet.nf'
include { FASTQC                        } from './modules/local/fastqc.nf'
include { TRIM_GALORE                   } from './modules/local/trim_galore.nf'
include { FASTQC as FASTQC_TRIM         } from './modules/local/fastqc.nf'
include { BISMARK_ALIGN                 } from './modules/local/bismark_align.nf'
include { PICARD_MARKDUPLICATES         } from './modules/local/picard_markduplicates.nf'
include { SAMTOOLS_FLAGSTAT             } from './modules/local/samtools_flagstat.nf'
include { BISMARK_METHYLATION_EXTRACTOR } from './modules/local/bismark_methyl_extractor.nf'
include { COVERAGE_FILTER               } from './modules/local/coverage_filter.nf'
include { MULTIQC                       } from './modules/local/multiqc.nf'
include { FEATURE_EXTRACTION            } from './modules/local/feature_extraction.nf'
include { DIFFERENTIAL_METHYLATION      } from './modules/local/differential_methylation.nf'
include { CLASSIFIER_TRAIN_EVAL         } from './modules/local/classifier_train_eval.nf'
include { REPORT_GENERATE               } from './modules/local/report_generate.nf'

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║     Universal DX  |  cfDNA Methylation Pipeline v1.0        ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  Input samplesheet : ${params.input}
    ║  Bismark index     : ${params.bismark_index}
    ║  Output dir        : ${params.outdir}
    ║  DM method         : ${params.dm_method}
    ║  ML method         : ${params.ml_method}
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()

    // 1. Validate samplesheet
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
    VALIDATE_SAMPLESHEET(ch_samplesheet)

    ch_reads = VALIDATE_SAMPLESHEET.out.validated_csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id          : row.sample_id,
                patient     : row.patient_id,
                condition   : row.condition,
                sample_type : row.sample_type,
                library     : row.library_type
            ]
            def reads = (row.read2 && row.read2.trim())
                ? [ file(row.read1), file(row.read2) ]
                : [ file(row.read1) ]
            [ meta, reads ]
        }

    // 2. Raw QC
    FASTQC(ch_reads)

    // 3. Adapter trimming
    TRIM_GALORE(ch_reads)
    FASTQC_TRIM(TRIM_GALORE.out.reads)

    // 4. Bisulfite alignment
    ch_index = Channel.value( file(params.bismark_index, checkIfExists: true) )
    BISMARK_ALIGN(TRIM_GALORE.out.reads, ch_index)

    // 5. Mark duplicates (flag, do NOT remove – cfDNA best practice)
    PICARD_MARKDUPLICATES(BISMARK_ALIGN.out.bam)

    // 6. Alignment QC
    SAMTOOLS_FLAGSTAT(PICARD_MARKDUPLICATES.out.bam)

    // 7. Methylation extraction
    BISMARK_METHYLATION_EXTRACTOR(
        PICARD_MARKDUPLICATES.out.bam,
        ch_index
    )

    // 8. Coverage filter
    COVERAGE_FILTER(
        BISMARK_METHYLATION_EXTRACTOR.out.bismark_cov,
        params.min_coverage
    )

    // 9. MultiQC
    ch_multiqc_inputs = FASTQC.out.zip
        .mix( FASTQC_TRIM.out.zip )
        .mix( TRIM_GALORE.out.log )
        .mix( BISMARK_ALIGN.out.report )
        .mix( PICARD_MARKDUPLICATES.out.metrics )
        .mix( SAMTOOLS_FLAGSTAT.out.flagstat )
        .mix( BISMARK_METHYLATION_EXTRACTOR.out.report )
        .collect()

    MULTIQC(ch_multiqc_inputs)

    // 10. Feature extraction → β-value matrix
    ch_coverage_files = COVERAGE_FILTER.out.filtered_cov.collect()
    FEATURE_EXTRACTION(
        ch_coverage_files,
        VALIDATE_SAMPLESHEET.out.validated_csv
    )

    // 11. Differential methylation
    DIFFERENTIAL_METHYLATION(
        FEATURE_EXTRACTION.out.beta_matrix,
        FEATURE_EXTRACTION.out.sample_metadata,
        params.dm_method,
        params.delta_beta,
        params.fdr_cutoff,
        params.min_cpgs_dmr
    )

    // 12. ML classifier
    CLASSIFIER_TRAIN_EVAL(
        FEATURE_EXTRACTION.out.beta_matrix,
        DIFFERENTIAL_METHYLATION.out.dmr_features,
        FEATURE_EXTRACTION.out.sample_metadata,
        params.ml_method,
        params.n_features,
        params.cv_folds,
        params.test_split,
        params.random_seed
    )

    // 13. Report
    REPORT_GENERATE(
        MULTIQC.out.report,
        FEATURE_EXTRACTION.out.qc_summary,
        DIFFERENTIAL_METHYLATION.out.results_table,
        CLASSIFIER_TRAIN_EVAL.out.metrics_json,
        CLASSIFIER_TRAIN_EVAL.out.roc_data
    )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n✅  Pipeline complete → ${params.outdir}"
        : "\n❌  Pipeline failed. See .nextflow.log" )
}
