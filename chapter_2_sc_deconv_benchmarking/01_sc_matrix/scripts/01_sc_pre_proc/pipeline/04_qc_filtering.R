#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 04_qc_filtering.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Filter low-quality cells using fixed cutoffs (unified across datasets)
# Input: Seurat objects with QC metrics (from script 03)
# Output: Filtered Seurat objects + filtering summary CSV
# Usage: Rscript 04_qc_filtering.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 04_qc_filtering.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Create output directory ----
dir.create(seurat_qc_filtered_dir, recursive = TRUE, showWarnings = FALSE)

# ---- List QC Seurat objects ----
seurat_files <- list.files(
  seurat_qc_metrics_dir,
  pattern = "_seurat_qc_metrics.rds$",
  full.names = TRUE
)

# ---- Exclude samples (from config) ----
sample_names_in_files <- sub("_seurat_qc_metrics.rds", "", basename(seurat_files))

# Biological exclusion (e.g., metastases, non-primary tumor)
if (exists("samples_to_exclude_bio") && length(samples_to_exclude_bio) > 0) {
  matched_bio <- samples_to_exclude_bio[samples_to_exclude_bio %in% sample_names_in_files]
  if (length(matched_bio) > 0) {
    cat("Excluding samples (biological):", paste(matched_bio, collapse = ", "), "\n")
    seurat_files <- seurat_files[!sample_names_in_files %in% matched_bio]
    sample_names_in_files <- sub("_seurat_qc_metrics.rds", "", basename(seurat_files))
  }
}

# Quality exclusion (e.g., failed QC, low cell count)
if (exists("samples_to_exclude_qc") && length(samples_to_exclude_qc) > 0) {
  matched_qc <- samples_to_exclude_qc[samples_to_exclude_qc %in% sample_names_in_files]
  if (length(matched_qc) > 0) {
    cat("Excluding samples (quality):", paste(matched_qc, collapse = ", "), "\n")
    seurat_files <- seurat_files[!sample_names_in_files %in% matched_qc]
  }
}

cat("\n=============================\n")
cat("Script: 04_qc_filtering.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Samples to process:", length(seurat_files), "\n")
cat("Filtering strategy: fixed cutoffs\n")
cat("-----------------------------\n")
cat("  nFeature_RNA: >", qc_min_features, "and <", qc_max_features, "\n")
cat("  nCount_RNA:   >", qc_min_counts, "\n")
cat("  percent.mt:   <", qc_max_percent_mt, "\n")
cat("  Remove doublets:", remove_doublets, "\n")

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

  # ---- Apply fixed cutoffs ----
  seu <- subset(
    seu,
    subset =
      nFeature_RNA > qc_min_features &
      nFeature_RNA < qc_max_features &
      nCount_RNA > qc_min_counts &
      percent.mt < qc_max_percent_mt
  )

  n_after_qc <- ncol(seu)
  cat("Cells after QC cutoffs:", n_after_qc, "\n")

  # ---- Remove doublets ----
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
    cells_after_qc = n_after_qc,
    doublets_removed = n_doublets_removed,
    cells_after = n_after,
    cells_removed = n_removed,
    pct_removed = pct_removed
  )

  # ---- Save ----
  outfile <- file.path(seurat_qc_filtered_dir, paste0(samp, "_seurat_filtered.rds"))
  saveRDS(seu, file = outfile)
  cat("Saved:", outfile, "\n")

  rm(seu)
  gc()
}

# ---- Save filtering summary ----
summary_df <- do.call(rbind, filter_summary)

# Add totals row
totals <- data.frame(
  sample = "TOTAL",
  cells_before = sum(summary_df$cells_before),
  cells_after_qc = sum(summary_df$cells_after_qc),
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
