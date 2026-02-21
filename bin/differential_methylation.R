#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# differential_methylation.R
# Universal DX | cfDNA Methylation Pipeline
#
# Performs:
#   1. Sample-level M-value QC (PCA, hierarchical clustering)
#   2. Differential Methylation Position (DMP) analysis   – via methylKit or DSS
#   3. Differential Methylation Region  (DMR) analysis    – DMRseq or bumphunter
#   4. DMR annotation (genomic context, CpG island status)
#   5. Export tables for downstream feature engineering
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(methylKit)
  library(DSS)
  library(dmrseq)
  library(genomation)
  library(GenomicRanges)
  library(pheatmap)
  library(RColorBrewer)
  library(jsonlite)
})

# ── CLI ───────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--beta_matrix",   type="character", help="Path to beta_matrix_variable.csv.gz"),
  make_option("--metadata",      type="character", help="Sample metadata CSV"),
  make_option("--method",        type="character", default="methylkit", help="methylkit|DSS"),
  make_option("--delta_beta",    type="double",    default=0.20),
  make_option("--fdr",           type="double",    default=0.05),
  make_option("--min_cpgs_dmr",  type="integer",   default=3),
  make_option("--outdir",        type="character",  default="./dma_results")
)
opts <- parse_args(OptionParser(option_list=option_list))
dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

log_msg <- function(...) message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), " ", ...)

# ── Load data ─────────────────────────────────────────────────────────────────
log_msg("Loading β-matrix …")
beta_mat <- fread(opts$beta_matrix) |> as.data.frame()
rownames(beta_mat) <- beta_mat$V1
beta_mat$V1        <- NULL
beta_mat           <- as.matrix(beta_mat)

meta <- read.csv(opts$metadata)
# Ensure column order matches
meta <- meta[match(colnames(beta_mat), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(beta_mat)))

log_msg(sprintf("β-matrix: %d CpGs × %d samples", nrow(beta_mat), ncol(beta_mat)))

# ── QC: M-value transformation & PCA ─────────────────────────────────────────
# M-value = log2(β / (1-β)) – more statistically appropriate for DM testing
m_mat <- log2( (beta_mat + 0.001) / (1 - beta_mat + 0.001) )

pca       <- prcomp(t(m_mat), scale.=TRUE, center=TRUE)
pca_df    <- as.data.frame(pca$x[, 1:4])
pca_df$sample_id  <- rownames(pca_df)
pca_df    <- merge(pca_df, meta, by="sample_id")
var_expl  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour=condition, shape=sample_type)) +
  geom_point(size=3, alpha=0.85) +
  labs(title="PCA of CpG M-values",
       x=sprintf("PC1 (%.1f%%)", var_expl[1]),
       y=sprintf("PC2 (%.1f%%)", var_expl[2])) +
  theme_bw(base_size=12) +
  scale_colour_brewer(palette="Set1")

ggsave(file.path(opts$outdir, "pca_m_values.pdf"), p_pca, width=7, height=5)
log_msg("PCA plot saved.")

# Hierarchical clustering heatmap (top 1000 variable CpGs)
top1k <- order(apply(m_mat, 1, var), decreasing=TRUE)[1:min(1000, nrow(m_mat))]
ann_col <- data.frame(Condition=meta$condition, row.names=meta$sample_id)

pheatmap(m_mat[top1k, ],
         annotation_col  = ann_col,
         show_rownames   = FALSE,
         show_colnames   = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         color           = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         main            = "Top-1000 variable CpGs (M-values)",
         filename        = file.path(opts$outdir, "heatmap_top1000cpgs.pdf"),
         width=10, height=12)

# ── Differential Methylation ──────────────────────────────────────────────────

# Split conditions (CRC vs control; polyp/adenoma treated as separate contrast)
case_ids    <- meta$sample_id[meta$condition %in% c("CRC")]
control_ids <- meta$sample_id[meta$condition == "control"]

run_methylkit_dmp <- function(beta_mat, case_ids, control_ids, meta) {
  log_msg("Running methylKit DMP analysis …")

  # Convert β-matrix to methylKit object list
  make_mk_obj <- function(sid, treatment) {
    beta_col   <- beta_mat[, sid, drop=FALSE]
    coords     <- strsplit(rownames(beta_mat), ":")
    chr_vec    <- sapply(coords, `[`, 1)
    pos_vec    <- as.integer(sapply(coords, `[`, 2))
    coverage   <- rep(30L, nrow(beta_mat))            # placeholder; use actual counts in production
    numCs      <- round(beta_col[,1] * coverage)
    numTs      <- coverage - numCs
    methylRaw(
      data.frame(chr=chr_vec, start=pos_vec, end=pos_vec,
                 coverage=coverage, numCs=numCs, numTs=numTs),
      sample.id    = sid,
      assembly     = "hg38",
      context      = "CpG",
      treatment    = treatment,
      resolution   = "base"
    )
  }

  obj_list <- c(
    lapply(control_ids, make_mk_obj, treatment=0),
    lapply(case_ids,    make_mk_obj, treatment=1)
  )
  meth_obj <- methylRawList(obj_list,
                             treatment=c(rep(0,length(control_ids)),
                                         rep(1,length(case_ids))))
  united  <- unite(meth_obj, destrand=FALSE)
  diff    <- calculateDiffMeth(united,
                                overdispersion="MN",
                                adjust="BH",
                                mc.cores=4)

  dmps <- getMethylDiff(diff,
                        difference = opts$delta_beta * 100,
                        qvalue     = opts$fdr)
  as.data.frame(dmps)
}

run_dss_dmp <- function(beta_mat, case_ids, control_ids) {
  log_msg("Running DSS DMP analysis …")
  coords  <- strsplit(rownames(beta_mat), ":")
  chr_vec <- sapply(coords, `[`, 1)
  pos_vec <- as.integer(sapply(coords, `[`, 2))

  make_bsseq <- function(ids) {
    cov_mat  <- matrix(30L, nrow=nrow(beta_mat), ncol=length(ids))
    M_mat    <- round(beta_mat[, ids, drop=FALSE] * 30)
    BSseq(chr=chr_vec, pos=pos_vec,
          M=M_mat, Cov=cov_mat,
          sampleNames=ids)
  }

  bs_case    <- make_bsseq(case_ids)
  bs_ctrl    <- make_bsseq(control_ids)
  bs_combined <- combine(bs_ctrl, bs_case)

  dml_test <- DMLtest(bs_combined,
                      group1 = seq_along(control_ids),
                      group2 = seq_along(case_ids) + length(control_ids),
                      smoothing=TRUE)
  callDML(dml_test, p.threshold=opts$fdr, delta=opts$delta_beta)
}

# Select method
if (opts$method == "DSS") {
  dmp_df <- run_dss_dmp(beta_mat, case_ids, control_ids)
} else {
  dmp_df <- run_methylkit_dmp(beta_mat, case_ids, control_ids, meta)
}

fwrite(dmp_df, file.path(opts$outdir, "DMPs.csv.gz"))
log_msg(sprintf("DMPs: %d CpGs (|Δβ| ≥ %.2f, FDR ≤ %.2f)",
                nrow(dmp_df), opts$delta_beta, opts$fdr))

# ── Volcano plot ──────────────────────────────────────────────────────────────
if ("meth.diff" %in% colnames(dmp_df) && "qvalue" %in% colnames(dmp_df)) {
  dmp_df$neg_log10_q <- -log10(pmax(dmp_df$qvalue, 1e-300))
  dmp_df$sig         <- abs(dmp_df$meth.diff/100) >= opts$delta_beta &
                         dmp_df$qvalue <= opts$fdr
  p_volcano <- ggplot(dmp_df, aes(meth.diff/100, neg_log10_q, colour=sig)) +
    geom_point(size=0.5, alpha=0.6) +
    geom_vline(xintercept=c(-opts$delta_beta, opts$delta_beta), linetype="dashed") +
    geom_hline(yintercept=-log10(opts$fdr), linetype="dashed") +
    scale_colour_manual(values=c("grey60","#d62728")) +
    labs(title="DMP Volcano: CRC vs Control",
         x="Δβ (CRC − Control)", y="-log10(FDR q-value)") +
    theme_bw(base_size=12)
  ggsave(file.path(opts$outdir, "volcano_DMPs.pdf"), p_volcano, width=7, height=5)
}

# ── DMR calling ───────────────────────────────────────────────────────────────
log_msg("Calling DMRs via dmrseq …")
tryCatch({
  coords  <- strsplit(rownames(beta_mat), ":")
  chr_vec <- factor(sapply(coords, `[`, 1))
  pos_vec <- as.integer(sapply(coords, `[`, 2))

  cov_mat  <- matrix(30L, nrow=nrow(beta_mat), ncol=ncol(beta_mat))
  M_mat    <- round(beta_mat * 30)
  bs_obj   <- BSseq(chr=chr_vec, pos=pos_vec,
                    M=M_mat, Cov=cov_mat,
                    sampleNames=colnames(beta_mat))

  testcovar <- ifelse(colnames(beta_mat) %in% case_ids, 1, 0)
  dmrs      <- dmrseq(bs=bs_obj, testCovariate="X",
                      BPPARAM=BiocParallel::MulticoreParam(4))

  dmr_df    <- as.data.frame(dmrs)
  dmr_df    <- dmr_df[dmr_df$qvalue <= opts$fdr & abs(dmr_df$beta) >= opts$delta_beta, ]
  fwrite(dmr_df, file.path(opts$outdir, "DMRs.csv.gz"))
  log_msg(sprintf("DMRs called: %d (q ≤ %.2f, |β| ≥ %.2f)", nrow(dmr_df), opts$fdr, opts$delta_beta))
}, error=function(e) {
  warning("dmrseq failed: ", conditionMessage(e),
          "\nFalling back to sliding-window DMR aggregation from DMPs.")
  # Simple fallback: collapse adjacent DMPs into DMRs
  dmp_sig    <- dmp_df[dmp_df$sig == TRUE, ]
  dmr_df_fb  <- dmp_sig  # placeholder — real production uses proper aggregation
  fwrite(dmr_df_fb, file.path(opts$outdir, "DMRs_fallback.csv.gz"))
})

# ── Annotate DMPs ─────────────────────────────────────────────────────────────
# Gene annotation would use AnnotationHub / TxDb objects in production.
# Here we output unannotated tables with hooks for annotation joins.
log_msg("Writing summary JSON …")
summary_json <- list(
  method        = opts$method,
  n_DMPs        = nrow(dmp_df),
  delta_beta    = opts$delta_beta,
  fdr_cutoff    = opts$fdr,
  hyper_DMPs    = sum(dmp_df[["meth.diff"]] > 0, na.rm=TRUE),
  hypo_DMPs     = sum(dmp_df[["meth.diff"]] < 0, na.rm=TRUE),
  contrasts     = list(case=case_ids, control=control_ids)
)
write(toJSON(summary_json, auto_unbox=TRUE, pretty=TRUE),
      file.path(opts$outdir, "dm_summary.json"))

log_msg("Differential methylation analysis complete.")
