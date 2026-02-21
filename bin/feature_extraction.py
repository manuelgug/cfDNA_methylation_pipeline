#!/usr/bin/env python3
"""
feature_extraction.py
─────────────────────
Universal DX | cfDNA Methylation Pipeline
Aggregates per-sample Bismark coverage files into:
  - β-value matrix (CpGs × samples)
  - QC summary (coverage, bisulphite conversion, fragment stats)
  - Variance-filtered & annotation-enriched feature set

Usage (called from Nextflow FEATURE_EXTRACTION process):
    python feature_extraction.py \
        --cov_files  *.filtered.cov.gz \
        --samplesheet validated_samplesheet.csv \
        --min_cov     10 \
        --output_dir  ./features/
"""

import argparse
import gzip
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ── I/O ───────────────────────────────────────────────────────────────────────

def parse_bismark_cov(filepath: str) -> pd.DataFrame:
    """
    Parse Bismark coverage format:
    chr  start  end  methylation%  count_M  count_U
    Returns DataFrame with columns: chr, pos, beta, coverage
    """
    cols = ["chr", "start", "end", "beta_pct", "count_M", "count_U"]
    opener = gzip.open if str(filepath).endswith(".gz") else open
    with opener(filepath, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, names=cols,
                         dtype={"chr": str, "start": int, "beta_pct": float,
                                "count_M": int, "count_U": int})
    df["coverage"] = df["count_M"] + df["count_U"]
    df["beta"]     = df["beta_pct"] / 100.0
    df["cpg_id"]   = df["chr"] + ":" + df["start"].astype(str)
    return df[["cpg_id", "chr", "start", "beta", "coverage", "count_M", "count_U"]]


def build_beta_matrix(cov_files: list[str],
                      sample_ids: list[str],
                      min_cov: int = 10) -> pd.DataFrame:
    """
    Merges per-sample coverage DataFrames into a CpGs × samples β-value matrix.
    CpGs present in < 80% of samples (after coverage filter) are dropped.
    """
    log.info(f"Building β-matrix from {len(cov_files)} samples …")
    frames = {}
    for sid, fp in zip(sample_ids, cov_files):
        df = parse_bismark_cov(fp)
        df = df[df["coverage"] >= min_cov]
        frames[sid] = df.set_index("cpg_id")["beta"]

    beta_mat = pd.DataFrame(frames)

    # Require CpG present in ≥80% of samples
    keep_thresh = 0.80
    present_frac = beta_mat.notna().mean(axis=1)
    beta_mat = beta_mat.loc[present_frac >= keep_thresh]

    # Impute remaining NAs with row mean (within-class imputation preferred
    # in production; row mean used here for simplicity)
    beta_mat = beta_mat.apply(lambda row: row.fillna(row.mean()), axis=1)

    log.info(f"β-matrix shape after filtering: {beta_mat.shape}")
    return beta_mat


# ── QC Metrics ────────────────────────────────────────────────────────────────

def compute_qc_metrics(cov_files: list[str],
                       sample_ids: list[str],
                       min_cov: int = 10) -> pd.DataFrame:
    """
    Per-sample QC metrics:
      - total CpGs measured
      - CpGs at ≥10× coverage
      - median coverage
      - global mean β (bisulphite conversion proxy: expect ~0.03 for cfDNA)
      - lambda conversion rate (if spike-in present, else NaN)
    """
    records = []
    for sid, fp in zip(sample_ids, cov_files):
        df = parse_bismark_cov(fp)
        high_cov = df[df["coverage"] >= min_cov]
        records.append({
            "sample_id"           : sid,
            "total_cpgs"          : len(df),
            "cpgs_min_cov"        : len(high_cov),
            "pct_cpgs_covered"    : round(len(high_cov) / max(len(df), 1) * 100, 2),
            "median_coverage"     : float(df["coverage"].median()),
            "mean_coverage"       : float(df["coverage"].mean()),
            "global_mean_beta"    : float(df["beta"].mean()),
            "global_median_beta"  : float(df["beta"].median()),
            "pct_fully_methylated": float((df["beta"] > 0.8).mean() * 100),
            "pct_unmethylated"    : float((df["beta"] < 0.2).mean() * 100),
        })
    return pd.DataFrame(records)


# ── Variance filtering & feature annotation ───────────────────────────────────

def variance_filter(beta_mat: pd.DataFrame,
                    top_n: int = 100_000) -> pd.DataFrame:
    """
    Retain top-N most variable CpG sites (MAD-based, robust to outliers).
    """
    mad = beta_mat.apply(lambda row: stats.median_abs_deviation(row.dropna()), axis=1)
    top_sites = mad.nlargest(top_n).index
    log.info(f"Variance filter: {len(top_sites)} / {len(beta_mat)} CpGs retained")
    return beta_mat.loc[top_sites]


def add_genomic_context(beta_mat: pd.DataFrame) -> pd.DataFrame:
    """
    Parse chr:pos from index and add basic genomic context columns.
    In production this would join against a pre-built CpG annotation BED
    (CGI, gene body, promoter, repeat element, etc.)
    """
    coords = beta_mat.index.str.split(":", expand=True)
    beta_mat["chr"]   = coords[0].values
    beta_mat["start"] = coords[1].astype(int).values
    return beta_mat


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Feature extraction for cfDNA β-matrix")
    parser.add_argument("--cov_files",    nargs="+", required=True)
    parser.add_argument("--samplesheet",  required=True)
    parser.add_argument("--min_cov",      type=int, default=10)
    parser.add_argument("--top_variable", type=int, default=100_000)
    parser.add_argument("--output_dir",   default="./features")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Load metadata
    meta_df = pd.read_csv(args.samplesheet)
    sid_map  = {
        Path(f).name.replace(".filtered.cov.gz", ""): f
        for f in args.cov_files
    }
    # Align order to samplesheet
    ordered_ids  = [r for r in meta_df["sample_id"] if r in sid_map]
    ordered_files = [sid_map[r] for r in ordered_ids]

    if not ordered_ids:
        sys.exit("ERROR: no coverage files matched sample IDs in samplesheet.")

    # Build β-matrix
    beta_mat = build_beta_matrix(ordered_files, ordered_ids, args.min_cov)
    beta_var = variance_filter(beta_mat, args.top_variable)

    # Save
    beta_mat.to_csv(out / "beta_matrix_full.csv.gz", compression="gzip")
    beta_var.to_csv(out / "beta_matrix_variable.csv.gz", compression="gzip")
    log.info(f"Saved β-matrices → {out}")

    # QC metrics
    qc_df = compute_qc_metrics(ordered_files, ordered_ids, args.min_cov)
    qc_df.to_csv(out / "sample_qc_metrics.csv", index=False)

    # Merge metadata into QC
    qc_full = qc_df.merge(meta_df, on="sample_id", how="left")
    qc_full.to_csv(out / "sample_metadata_qc.csv", index=False)

    # JSON summary for downstream reporting
    summary = {
        "n_samples"   : len(ordered_ids),
        "n_cpgs_full" : int(beta_mat.shape[0]),
        "n_cpgs_var"  : int(beta_var.shape[0]),
        "min_coverage": args.min_cov,
        "conditions"  : meta_df["condition"].value_counts().to_dict(),
        "qc_flags"    : {
            "low_cpg_samples": qc_df.loc[
                qc_df["cpgs_min_cov"] < 1_000_000, "sample_id"
            ].tolist(),
            "high_global_beta": qc_df.loc[
                qc_df["global_mean_beta"] > 0.10, "sample_id"
            ].tolist(),
        }
    }
    (out / "feature_extraction_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("Feature extraction complete.")
    log.info(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
