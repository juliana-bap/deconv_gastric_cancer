#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 04_qc_filtering_mad.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Filter low-quality cells using fixed cutoffs (unified)
#              + MAD-based outlier removal as intra-sample safety net.
#              Fixed cutoffs are applied FIRST (same for all datasets),
#              then MAD removes extreme outliers within each sample.
# Input: Seurat objects with QC metrics (from script 03)
# Output: Filtered Seurat objects + filtering summary CSV
# Usage: Rscript 04_qc_filtering_mad.R <config_path>
# Author: Juliana Pinto
# Date: 2026
#
# Config parameters used:
#   qc_min_features, qc_max_features, qc_min_counts, qc_max_percent_mt
#   remove_doublets
#   mad_threshold (default: 5 — conservative, only catches extreme outliers)
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 04_qc_filtering_mad.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- MAD threshold (default 5 if not in config) ----
if (!exists("mad_threshold")) mad_threshold <- 5

# ---- Helper: flag outliers by MAD ----
is_outlier_mad <- function(x, nmads = mad_threshold, type = "both") {
  med <- median(x)
  mad_val <- mad(x, constant = 1.4826)

  if (mad_val == 0) return(rep(FALSE, length(x)))

  if (type == "both") {
    return(x < (med - nmads * mad_val) | x > (med + nmads * mad_val))
  } else if (type == "upper") {
    return(x > (med + nmads * mad_val))
  } else if (type == "lower") {
    return(x < (med - nmads * mad_val))
  }
}

# ---- Create output directory ----
dir.create(seurat_qc_filtered_dir, recursive = TRUE, showWarnings = FALSE)

# ---- List QC Seurat objects ----
seurat_files <- list.files(
  seurat_qc_metrics_dir,
  pattern = "_seurat_qc_metrics.rds$",
  full.names = TRUE
)

cat("\n=============================\n")
cat("Script: 04_qc_filtering_mad.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Samples found:", length(seurat_files), "\n")
cat("Filtering strategy: fixed cutoffs + MAD safety net\n")
cat("-----------------------------\n")
cat("  Fixed cutoffs:\n")
cat("    nFeature_RNA: >", qc_min_features, "and <", qc_max_features, "\n")
cat("    nCount_RNA:   >", qc_min_counts, "\n")
cat("    percent.mt:   <", qc_max_percent_mt, "\n")
cat("    Remove doublets:", remove_doublets, "\n")
cat("  MAD threshold:", mad_threshold, "(intra-sample outlier removal)\n")

# ---- Storage for summary ----
filter_summary <- list()

# ---- Loop over samples ----
total_samples <- length(seurat_files)

for (i in seq_along(seurat_files)) {

  file <- seurat_files[i]
  samp <- sub("_seurat_qc_metrics.rds", "", basename(file))

  cat("\n=============================\n")
  cat(sprintf("Sample [%d/%d]: %s\n", i, total_samples, samp))

  seu <- readRDS(file)
  n_before <- ncol(seu)
  cat("Cells before filtering:", n_before, "\n")

  # ---- Step 1: Apply fixed cutoffs ----
  seu <- subset(
    seu,
    subset =
      nFeature_RNA > qc_min_features &
      nFeature_RNA < qc_max_features &
      nCount_RNA > qc_min_counts &
      percent.mt < qc_max_percent_mt
  )

  n_after_fixed <- ncol(seu)
  cat("Cells after fixed cutoffs:", n_after_fixed, "\n")

  # ---- Step 2: MAD-based outlier removal (intra-sample) ----
  # Applied on the surviving cells after fixed cutoffs
  # Only flags extreme outliers (default 5 MADs)
  df <- seu@meta.data

  outlier_features <- is_outlier_mad(df$nFeature_RNA, type = "both")
  outlier_counts   <- is_outlier_mad(df$nCount_RNA, type = "both")
  outlier_mt       <- is_outlier_mad(df$percent.mt, type = "upper")

  mad_outliers <- outlier_features | outlier_counts | outlier_mt
  n_mad_removed <- sum(mad_outliers)

  if (n_mad_removed > 0) {
    cat("MAD outliers removed:", n_mad_removed, "\n")
    cat("  - nFeature_RNA:", sum(outlier_features), "\n")
    cat("  - nCount_RNA:", sum(outlier_counts), "\n")
    cat("  - percent.mt:", sum(outlier_mt), "\n")

    cells_keep <- rownames(df)[!mad_outliers]
    seu <- subset(seu, cells = cells_keep)
  } else {
    cat("MAD outliers removed: 0\n")
  }

  n_after_mad <- ncol(seu)

  # ---- Step 3: Remove doublets ----
  n_doublets_removed <- 0

  if (remove_doublets && "scDblFinder.class" %in% colnames(seu@meta.data)) {
    n_before_doublet <- ncol(seu)
    seu <- subset(seu, subset = scDblFinder.class == "singlet")
    n_doublets_removed <- n_before_doublet - ncol(seu)
    cat("Doublets removed:", n_doublets_removed, "\n")
  }

  n_after <- ncol(seu)
  n_removed <- n_before - n_after
  pct_removed <- round(100 * n_removed / n_before, 1)

  cat("Cells after filtering:", n_after, "\n")
  cat("Total removed:", n_removed, "(", pct_removed, "%)\n")

  # ---- Store summary ----
  filter_summary[[samp]] <- data.frame(
    sample = samp,
    cells_before = n_before,
    cells_after_fixed = n_after_fixed,
    mad_outliers = n_mad_removed,
    cells_after_mad = n_after_mad,
    doublets_removed = n_doublets_removed,
    cells_after = n_after,
    cells_removed = n_removed,
    pct_removed = pct_removed
  )

  # ---- Save ----
  outfile <- file.path(seurat_qc_filtered_dir, paste0(samp, "_seurat_filtered.rds"))
  saveRDS(seu, file = outfile)
  cat("Saved:", outfile, "\n")

  rm(seu, df)
  gc()
}

# ---- Save filtering summary ----
summary_df <- do.call(rbind, filter_summary)

totals <- data.frame(
  sample = "TOTAL",
  cells_before = sum(summary_df$cells_before),
  cells_after_fixed = sum(summary_df$cells_after_fixed),
  mad_outliers = sum(summary_df$mad_outliers),
  cells_after_mad = sum(summary_df$cells_after_mad),
  doublets_removed = sum(summary_df$doublets_removed),
  cells_after = sum(summary_df$cells_after),
  cells_removed = sum(summary_df$cells_removed),
  pct_removed = round(100 * sum(summary_df$cells_removed) / sum(summary_df$cells_before), 1)
)

summary_with_total <- rbind(summary_df, totals)

summary_path <- file.path(seurat_qc_filtered_dir, "QC_filtering_summary.csv")
write.csv(summary_with_total, summary_path, row.names = FALSE)

cat("\n=============================\n")
cat("Filtering summary:\n")
print(summary_with_total)
cat("\nSummary saved:", summary_path, "\n")
cat("Total cells:", totals$cells_before, "->", totals$cells_after,
    "(", totals$pct_removed, "% removed)\n")
cat("\nDone! Filtered objects saved to:", seurat_qc_filtered_dir, "\n")
