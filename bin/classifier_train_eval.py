#!/usr/bin/env python3
"""
classifier_train_eval.py
─────────────────────────
Universal DX | cfDNA Methylation Pipeline

Trains and evaluates ML classifiers for CRC early detection from cfDNA
methylation features. Implements:

  - Feature selection (DMR-guided + variance filter)
  - Cross-validated training (RF, XGBoost, Logistic Regression, Ensemble)
  - Performance evaluation: AUROC, AUPRC, sensitivity at fixed specificity
  - Bootstrap confidence intervals
  - Calibration curves
  - SHAP-based feature importance
  - Regulatory-grade output (per-sample scores, performance metrics JSON)

Usage:
    python classifier_train_eval.py \
        --beta_matrix   beta_matrix_variable.csv.gz \
        --dmr_features  DMRs.csv.gz \
        --metadata      sample_metadata_qc.csv \
        --method        ensemble \
        --n_features    500 \
        --cv_folds      5 \
        --test_split    0.20 \
        --seed          42 \
        --outdir        ./classifier_results/
"""

from __future__ import annotations

import argparse
import json
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import Pipeline
from sklearn.model_selection import (StratifiedKFold, StratifiedShuffleSplit,
                                     cross_val_predict)
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import (roc_auc_score, average_precision_score,
                              roc_curve, precision_recall_curve,
                              confusion_matrix, classification_report)
from sklearn.inspection import permutation_importance

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    warnings.warn("xgboost not installed; falling back to GradientBoosting.")

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RNG = np.random.default_rng(42)


# ── Helpers ───────────────────────────────────────────────────────────────────

def bootstrap_ci(y_true, y_score, metric_fn, n_boot=1000, ci=0.95):
    """Return (point_est, lower, upper) with bootstrap CI."""
    scores = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = RNG.integers(0, n, size=n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        scores.append(metric_fn(y_true[idx], y_score[idx]))
    scores = np.array(scores)
    alpha = (1 - ci) / 2
    return (metric_fn(y_true, y_score),
            float(np.quantile(scores, alpha)),
            float(np.quantile(scores, 1 - alpha)))


def sensitivity_at_specificity(y_true, y_score, target_spec=0.95):
    """Clinical metric: sensitivity at 95% specificity."""
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    specificities = 1 - fpr
    idx = np.searchsorted(specificities[::-1], target_spec)
    idx = len(specificities) - idx - 1
    idx = max(0, min(idx, len(tpr) - 1))
    return float(tpr[idx]), float(thresholds[idx])


def ppv_npv(y_true, y_pred, prevalence=None):
    """Positive / Negative predictive value, optionally prevalence-adjusted."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    ppv = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    npv = tn / (tn + fn) if (tn + fn) > 0 else float("nan")
    if prevalence is not None:
        # Bayes adjustment for population prevalence
        ppv = (ppv * prevalence) / (ppv * prevalence + (1 - ppv) * (1 - prevalence))
        npv = (npv * (1 - prevalence)) / (npv * (1 - prevalence) + (1 - npv) * prevalence)
    return ppv, npv


# ── Feature selection ─────────────────────────────────────────────────────────

def select_features(beta_mat: pd.DataFrame,
                    dmr_df: pd.DataFrame | None,
                    n_features: int,
                    labels: pd.Series) -> list[str]:
    """
    Tiered feature selection:
      1. DMR-anchored CpGs (if available) — biologically motivated
      2. Supplemented with top-variance sites to reach n_features
    """
    selected: set[str] = set()

    if dmr_df is not None and not dmr_df.empty:
        # Extract CpG IDs overlapping significant DMRs
        dmr_cpgs = set()
        for _, row in dmr_df.iterrows():
            try:
                chrom  = str(row.get("chr", row.get("seqnames", "")))
                start  = int(row.get("start", 0))
                end    = int(row.get("end",   start))
                matching = [c for c in beta_mat.index
                            if c.startswith(chrom + ":") and
                               start <= int(c.split(":")[1]) <= end]
                dmr_cpgs.update(matching)
            except Exception:
                continue
        selected.update(list(dmr_cpgs)[:n_features])
        log.info(f"DMR-anchored CpGs: {len(selected)}")

    # Top-variance supplement
    remaining = n_features - len(selected)
    if remaining > 0:
        X = beta_mat.loc[~beta_mat.index.isin(selected)]
        mad = X.apply(lambda r: float(np.median(np.abs(r - np.median(r)))), axis=1)
        top_var = mad.nlargest(remaining).index.tolist()
        selected.update(top_var)

    return [f for f in selected if f in beta_mat.index][:n_features]


# ── Model factory ─────────────────────────────────────────────────────────────

def build_model(method: str, seed: int) -> Pipeline:
    if method == "rf":
        clf = RandomForestClassifier(n_estimators=500, max_depth=None,
                                     min_samples_leaf=3,
                                     class_weight="balanced",
                                     n_jobs=-1, random_state=seed)
    elif method == "xgb":
        if HAS_XGB:
            clf = xgb.XGBClassifier(n_estimators=300, max_depth=4,
                                     learning_rate=0.05, subsample=0.8,
                                     colsample_bytree=0.8,
                                     use_label_encoder=False,
                                     eval_metric="logloss",
                                     random_state=seed, n_jobs=-1)
        else:
            clf = GradientBoostingClassifier(n_estimators=300, max_depth=4,
                                              learning_rate=0.05, subsample=0.8,
                                              random_state=seed)
    elif method == "logistic":
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=0.1,
                                  class_weight="balanced", max_iter=1000,
                                  random_state=seed)
    elif method == "ensemble":
        rf_clf  = RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                          n_jobs=-1, random_state=seed)
        lr_clf  = LogisticRegression(penalty="l2", C=0.1, class_weight="balanced",
                                      solver="lbfgs", max_iter=1000, random_state=seed)
        gb_base = (xgb.XGBClassifier(n_estimators=200, random_state=seed, n_jobs=-1)
                   if HAS_XGB else
                   GradientBoostingClassifier(n_estimators=200, random_state=seed))
        clf = VotingClassifier(
            estimators=[("rf", rf_clf), ("lr", lr_clf), ("gb", gb_base)],
            voting="soft",
            n_jobs=-1
        )
    else:
        raise ValueError(f"Unknown method: {method}")

    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", CalibratedClassifierCV(clf, cv=3, method="isotonic"))
    ])


# ── Evaluation ────────────────────────────────────────────────────────────────

def evaluate(y_true, y_score, y_pred, label=""):
    auc_roc, auc_lo, auc_hi = bootstrap_ci(y_true, y_score, roc_auc_score)
    auc_pr                  = average_precision_score(y_true, y_score)
    sens_95, thresh_95      = sensitivity_at_specificity(y_true, y_score, 0.95)
    sens_90, _              = sensitivity_at_specificity(y_true, y_score, 0.90)
    ppv, npv                = ppv_npv(y_true, y_pred, prevalence=0.034)  # CRC ~3.4% in screening pop.

    metrics = {
        "label"         : label,
        "n_samples"     : int(len(y_true)),
        "n_pos"         : int(y_true.sum()),
        "auroc"         : round(auc_roc, 4),
        "auroc_ci_low"  : round(auc_lo,  4),
        "auroc_ci_high" : round(auc_hi,  4),
        "auprc"         : round(auc_pr,  4),
        "sensitivity_at_95spec" : round(sens_95, 4),
        "sensitivity_at_90spec" : round(sens_90, 4),
        "ppv_prevalence_adjusted": round(ppv, 4) if not np.isnan(ppv) else None,
        "npv_prevalence_adjusted": round(npv, 4) if not np.isnan(npv) else None,
    }
    log.info(f"[{label}] AUROC={auc_roc:.3f} [{auc_lo:.3f}–{auc_hi:.3f}]  "
             f"Sens@95spec={sens_95:.3f}  AUPRC={auc_pr:.3f}")
    return metrics


# ── Plots ─────────────────────────────────────────────────────────────────────

def plot_roc(y_true, y_score, label, outdir, tag=""):
    fpr, tpr, _ = roc_curve(y_true, y_score)
    auc_val     = roc_auc_score(y_true, y_score)
    fig, ax     = plt.subplots(figsize=(5, 5))
    ax.plot(fpr, tpr, lw=2, label=f"AUROC = {auc_val:.3f}")
    ax.plot([0,1],[0,1], "--", color="grey", lw=1)
    ax.axvline(0.05, color="red", linestyle=":", alpha=0.7, label="95% specificity")
    ax.set_xlabel("False Positive Rate (1 – Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title(f"ROC Curve  |  {label}")
    ax.legend(loc="lower right")
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1.01])
    plt.tight_layout()
    plt.savefig(outdir / f"roc_curve{tag}.pdf", dpi=150)
    plt.close()
    return {"fpr": fpr.tolist(), "tpr": tpr.tolist(), "auroc": float(auc_val)}


def plot_prc(y_true, y_score, label, outdir, tag=""):
    prec, rec, _ = precision_recall_curve(y_true, y_score)
    auprc        = average_precision_score(y_true, y_score)
    prevalence   = y_true.mean()
    fig, ax      = plt.subplots(figsize=(5, 5))
    ax.plot(rec, prec, lw=2, label=f"AUPRC = {auprc:.3f}")
    ax.axhline(prevalence, color="grey", linestyle="--", label=f"Prevalence = {prevalence:.3f}")
    ax.set_xlabel("Recall (Sensitivity)")
    ax.set_ylabel("Precision (PPV)")
    ax.set_title(f"Precision–Recall Curve  |  {label}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"prc_curve{tag}.pdf", dpi=150)
    plt.close()


def plot_calibration(y_true, y_score, label, outdir, tag=""):
    frac_pos, mean_pred = calibration_curve(y_true, y_score, n_bins=10)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(mean_pred, frac_pos, "s-", label="Model")
    ax.plot([0,1],[0,1],"--", color="grey", label="Perfectly calibrated")
    ax.set_xlabel("Mean Predicted Probability"); ax.set_ylabel("Fraction Positives")
    ax.set_title(f"Calibration Curve  |  {label}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"calibration{tag}.pdf", dpi=150)
    plt.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta_matrix",  required=True)
    parser.add_argument("--dmr_features", default=None)
    parser.add_argument("--metadata",     required=True)
    parser.add_argument("--method",       default="ensemble")
    parser.add_argument("--n_features",   type=int, default=500)
    parser.add_argument("--cv_folds",     type=int, default=5)
    parser.add_argument("--test_split",   type=float, default=0.20)
    parser.add_argument("--seed",         type=int, default=42)
    parser.add_argument("--outdir",       default="./classifier_results")
    args = parser.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load
    log.info("Loading data …")
    beta_mat = pd.read_csv(args.beta_matrix, index_col=0)
    meta     = pd.read_csv(args.metadata)
    meta     = meta[meta["sample_id"].isin(beta_mat.columns)]
    beta_mat = beta_mat[meta["sample_id"].tolist()]

    dmr_df = None
    if args.dmr_features and Path(args.dmr_features).exists():
        dmr_df = pd.read_csv(args.dmr_features)

    # Binary label: CRC=1, control=0; exclude non-binary classes from main model
    label_map = {"CRC": 1, "control": 0}
    meta_bin  = meta[meta["condition"].isin(label_map)].copy()
    meta_bin["label"] = meta_bin["condition"].map(label_map)
    beta_bin  = beta_mat[meta_bin["sample_id"].tolist()]

    log.info(f"Samples: {meta_bin['label'].value_counts().to_dict()}")

    # ── Feature selection
    feat_ids = select_features(beta_bin, dmr_df, args.n_features, meta_bin["label"])
    X        = beta_bin.loc[feat_ids].T.values.astype(np.float32)
    y        = meta_bin["label"].values
    log.info(f"Feature matrix: {X.shape}")

    # ── Hold-out split (stratified)
    sss = StratifiedShuffleSplit(n_splits=1, test_size=args.test_split,
                                  random_state=args.seed)
    train_idx, test_idx = next(sss.split(X, y))
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    # ── Cross-validated OOF predictions (training set)
    log.info(f"CV ({args.cv_folds}-fold) training with method='{args.method}' …")
    model_cv   = build_model(args.method, args.seed)
    cv         = StratifiedKFold(n_splits=args.cv_folds, shuffle=True,
                                  random_state=args.seed)
    oof_scores = cross_val_predict(model_cv, X_train, y_train,
                                   cv=cv, method="predict_proba")[:, 1]
    oof_preds  = (oof_scores >= 0.5).astype(int)
    cv_metrics = evaluate(y_train, oof_scores, oof_preds, label="CV OOF")

    # ── Final model on full training set
    log.info("Training final model on full training set …")
    model_final = build_model(args.method, args.seed)
    model_final.fit(X_train, y_train)

    # ── Hold-out evaluation
    test_scores = model_final.predict_proba(X_test)[:, 1]
    test_preds  = (test_scores >= 0.5).astype(int)
    ho_metrics  = evaluate(y_test, test_scores, test_preds, label="Hold-out")

    # ── Per-sample scores table
    all_scores = model_final.predict_proba(X)[:, 1]
    score_df   = pd.DataFrame({
        "sample_id"       : meta_bin["sample_id"].values,
        "condition"       : meta_bin["condition"].values,
        "label"           : y,
        "pred_score"      : all_scores,
        "pred_label"      : (all_scores >= 0.5).astype(int),
        "split"           : ["train" if i in train_idx else "test"
                             for i in range(len(y))]
    })
    score_df.to_csv(out / "per_sample_scores.csv", index=False)

    # ── Plots
    roc_data = plot_roc(y_test, test_scores, "Hold-out set", out, "_holdout")
    plot_prc(y_test, test_scores, "Hold-out set", out, "_holdout")
    plot_calibration(y_test, test_scores, "Hold-out set", out, "_holdout")
    plot_roc(y_train, oof_scores, "Cross-validation OOF", out, "_cv")

    # ── Feature importance (permutation, then SHAP if available)
    perm_imp = permutation_importance(
        model_final, X_test, y_test, n_repeats=20,
        random_state=args.seed, scoring="roc_auc"
    )
    feat_imp_df = pd.DataFrame({
        "feature"   : feat_ids,
        "importance": perm_imp.importances_mean,
        "std"       : perm_imp.importances_std
    }).sort_values("importance", ascending=False)
    feat_imp_df.to_csv(out / "feature_importance.csv", index=False)

    # Top-30 importance barplot
    top30 = feat_imp_df.head(30)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top30["feature"][::-1], top30["importance"][::-1],
            xerr=top30["std"][::-1], capsize=3, color="#1f77b4", alpha=0.8)
    ax.set_xlabel("Permutation Importance (ΔAUROC)")
    ax.set_title("Top-30 CpG Features by Permutation Importance")
    plt.tight_layout()
    plt.savefig(out / "feature_importance_top30.pdf", dpi=150)
    plt.close()

    # SHAP (if available and model is RF)
    if HAS_SHAP and args.method in ("rf", "xgb"):
        try:
            log.info("Computing SHAP values …")
            explainer  = shap.TreeExplainer(
                model_final.named_steps["clf"].calibrated_classifiers_[0].estimator
            )
            shap_vals  = explainer.shap_values(
                model_final.named_steps["scaler"].transform(X_test)
            )
            sv = shap_vals[1] if isinstance(shap_vals, list) else shap_vals
            shap.summary_plot(sv, X_test, feature_names=feat_ids,
                               show=False, max_display=30)
            plt.savefig(out / "shap_summary.pdf", bbox_inches="tight", dpi=150)
            plt.close()
        except Exception as e:
            log.warning(f"SHAP failed: {e}")

    # ── Final JSON output
    metrics_out = {
        "pipeline_version"      : "1.0.0",
        "method"                : args.method,
        "n_features_selected"   : len(feat_ids),
        "n_train"               : int(len(train_idx)),
        "n_test"                : int(len(test_idx)),
        "cross_validation"      : cv_metrics,
        "holdout"               : ho_metrics,
        "classification_report" : classification_report(y_test, test_preds,
                                    target_names=["control","CRC"],
                                    output_dict=True),
    }
    (out / "classifier_metrics.json").write_text(
        json.dumps(metrics_out, indent=2, default=str)
    )
    (out / "roc_data.json").write_text(json.dumps(roc_data, indent=2))

    log.info(f"✅  Classifier pipeline complete. Results → {out}")


if __name__ == "__main__":
    main()
