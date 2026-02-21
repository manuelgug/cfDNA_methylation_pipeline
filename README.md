# Universal DX | cfDNA Methylation Pipeline

> **End-to-end Nextflow DSL2 pipeline for liquid-biopsy cancer early detection**  
> QC → Alignment → Methylation Extraction → Feature Engineering → Differential Methylation → ML Classifier → Report

---

## Overview

This pipeline implements the complete analytical stack for cfDNA methylation-based
colorectal cancer (CRC) detection from NGS data:

```
Raw FASTQ reads
    │
    ▼
[FastQC] Raw QC
    │
    ▼
[Trim Galore] Adapter trimming (bisulfite-aware)
    │
    ▼
[FastQC] Post-trim QC
    │
    ▼
[Bismark] Bisulfite alignment → hg38
    │
    ▼
[Picard] Mark duplicates (flagged, not removed — cfDNA best practice)
    │
    ▼
[samtools flagstat] Alignment QC
    │
    ▼
[Bismark extractor] CpG methylation extraction (CX report)
    │
    ▼
[Coverage filter] Min 10× CpG sites
    │
    ▼
[MultiQC] Aggregate QC report
    │
    ▼
[Feature extraction] β-matrix + variance filter + QC metrics
    │
    ▼
[Differential methylation] DMPs + DMRs (methylKit/DSS + dmrseq)
    │
    ├──▶ Volcano plots, PCA, heatmaps
    │
    ▼
[ML Classifier] Ensemble (RF + LR + XGB)
    │
    ├──▶ ROC/PRC curves, calibration, SHAP importance
    │
    ▼
[Report] Self-contained HTML + JSON summary
```

---

## Quick Start

### 1. Requirements

- Nextflow ≥ 23.04
- Docker or Singularity
- Java 11+

### 2. Prepare samplesheet

Create a CSV with the following columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `patient_id` | Patient / donor identifier |
| `condition` | `CRC` \| `control` \| `adenoma` \| `polyp` |
| `sample_type` | `cfDNA` \| `tissue` |
| `library_type` | `WGBS` \| `targeted` \| `hybrid_capture` \| `RRBS` |
| `read1` | Absolute path to R1 FASTQ (.gz) |
| `read2` | Absolute path to R2 FASTQ (.gz) — omit for SE |

See `tests/samplesheet_test.csv` for an example.

### 3. Prepare Bismark genome index

```bash
bismark_genome_preparation --bowtie2 /path/to/hg38/
```

### 4. Run

```bash
# Local with Docker
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv \
    --bismark_index /path/to/bismark_index/ \
    --genome /path/to/hg38.fa \
    --outdir results/

# AWS Batch
nextflow run main.nf \
    -profile aws \
    --input s3://bucket/samplesheet.csv \
    --bismark_index s3://bucket/bismark_index/ \
    --outdir s3://bucket/results/

# Test run
nextflow run main.nf -profile docker,test
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_coverage` | 10 | Minimum CpG read depth |
| `--dm_method` | `methylkit` | DM method: `methylkit` or `DSS` |
| `--delta_beta` | 0.20 | \|Δβ\| threshold for DMP/DMR calling |
| `--fdr_cutoff` | 0.05 | FDR threshold |
| `--ml_method` | `ensemble` | Classifier: `rf`, `xgb`, `logistic`, `ensemble` |
| `--n_features` | 500 | CpG features for classifier |
| `--cv_folds` | 5 | Cross-validation folds |
| `--test_split` | 0.20 | Hold-out fraction |

---

## Output Structure

```
results/
├── pipeline_info/          # Nextflow execution timeline, trace, DAG
├── qc/
│   ├── fastqc/             # Per-sample FastQC (raw + trimmed)
│   └── flagstat/           # Alignment statistics
├── trimmed/                # Trimmed FASTQ logs
├── aligned/                # BAM files + duplication metrics
├── methylation/            # Bismark coverage + CX reports
├── methylation_filtered/   # Coverage-filtered CpG tables
├── multiqc/                # Aggregated QC report
├── features/
│   ├── beta_matrix_full.csv.gz       # Full β-value matrix
│   ├── beta_matrix_variable.csv.gz   # Top-variance CpGs
│   ├── sample_qc_metrics.csv         # Per-sample QC stats
│   └── feature_extraction_summary.json
├── differential_methylation/
│   ├── DMPs.csv.gz          # Differential methylation positions
│   ├── DMRs.csv.gz          # Differential methylation regions
│   ├── volcano_DMPs.pdf     # Volcano plot
│   ├── pca_m_values.pdf     # PCA of M-values
│   ├── heatmap_top1000cpgs.pdf
│   └── dm_summary.json
├── classifier/
│   ├── classifier_metrics.json  # Full performance metrics + CIs
│   ├── roc_data.json
│   ├── per_sample_scores.csv
│   ├── feature_importance.csv
│   ├── roc_curve_holdout.pdf
│   ├── prc_curve_holdout.pdf
│   ├── calibration_holdout.pdf
│   └── shap_summary.pdf
└── report/
    ├── udx_pipeline_report.html  # ← Main deliverable
    └── udx_pipeline_summary.json
```

---

## Scientific Rationale

### Why M-values for DM testing?
β-values are bounded [0,1] and heteroscedastic — inappropriate for linear models.
M-values (logit transform) are approximately normally distributed and homoscedastic,
making them better suited for parametric DM tests.

### Why not remove duplicates in cfDNA?
cfDNA fragments are short (~160 bp) and often share true start/end positions due to
biology (apoptotic cleavage at nucleosome linkers), not just PCR artefacts. Removing
all "duplicates" discards genuine biological signal. Best practice is to flag and
account for duplication in downstream analysis.

### Classifier ensemble rationale
- **Random Forest**: robust to feature correlation, provides permutation importance
- **Logistic Regression (L1/L2)**: interpretable, good calibration, handles high-dimensional sparse data
- **XGBoost**: captures non-linear interactions, handles class imbalance well
- **Soft voting ensemble**: reduces variance while preserving calibration

### Performance metrics for clinical context
AUROC alone is insufficient for screening tests with low prevalence (~3.4% for CRC).
The pipeline reports:
- Sensitivity at fixed 95% and 90% specificity (FDA/CLSI convention)
- AUPRC (better for imbalanced classes)
- PPV/NPV adjusted to screening population prevalence
- Bootstrap 95% confidence intervals (n=1000)

---

## Regulatory Considerations

This pipeline is designed with IVD/LDT development in mind:
- All parameters versioned in `nextflow.config`
- Execution trace with task hash, CPU/memory, exit codes
- Containerized execution for reproducibility
- JSON metric exports for downstream validation documentation
- Samplesheet validation step catches input errors before compute starts

---

## Dependencies

### Bioinformatics tools
| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Raw / post-trim QC |
| Trim Galore | 0.6.10 | Adapter trimming |
| Bismark | 0.24.0 | Bisulfite alignment |
| Picard | 2.27.5 | Duplicate marking |
| samtools | 1.18 | BAM handling |
| MultiQC | 1.19 | QC aggregation |

### R packages
`methylKit`, `DSS`, `dmrseq`, `genomation`, `GenomicRanges`, `pheatmap`, `ggplot2`

### Python packages
`numpy`, `pandas`, `scipy`, `scikit-learn`, `xgboost`, `shap`, `matplotlib`

---

## Citation

If using this pipeline in published work, please cite the underlying tools (Bismark, methylKit, dmrseq, etc.) and reference this pipeline repository.

---
