#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01b_patch_mt_GSE275648.R
# Project: Gastric Cancer Deconvolution
# Dataset: GSE275648
# Description: Patches percent.mt in post-QC Seurat objects using pre-computed
#              values saved during import (01_import_GSE275648.R).
#
#              WHY THIS EXISTS:
#              The authors removed MT genes from the count matrix before GEO
#              deposit, so script 03 (PercentageFeatureSet) returns 0 for all
#              cells. The real percent.mt values were pre-computed by the authors
#              and are stored in _orig_meta.rds (metadata_dir). This script
#              overwrites the zeros with those real values so the information
#              propagates through all downstream steps (merge, integration, etc.)
#
#              RUN ORDER: after 03_qc_metrics.R, before 04_filter.R
#              Modifies objects in seurat_qc_metrics_dir in place.
#
# Usage: Rscript 01b_patch_mt_GSE275648.R <config_path>
# Author: Juliana Pinto
# ==============================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01b_patch_mt_GSE275648.R <config_path>")
}

config_path <- normalizePath(args[1])
source(config_path)

library(Seurat)

cat("\n=============================\n")
cat("Script: 01b_patch_mt_GSE275648.R\n")
cat("Dataset:", dataset_id, "\n")
cat("QC metrics dir:", seurat_qc_metrics_dir, "\n")
cat("Metadata dir:", metadata_dir, "\n")

qc_files <- list.files(seurat_qc_metrics_dir,
                       pattern = "_seurat_qc_metrics\\.rds$",
                       full.names = TRUE)

if (length(qc_files) == 0) {
  stop("No _seurat_qc_metrics.rds files found in: ", seurat_qc_metrics_dir,
       "\nRun 03_qc_metrics.R first.")
}

cat("Samples found:", length(qc_files), "\n")

patched <- 0
skipped <- 0

for (f in qc_files) {

  samp <- sub("_seurat_qc_metrics\\.rds$", "", basename(f))

  cat("\n-----------------------------\n")
  cat("Sample:", samp, "\n")

  orig_meta_file <- file.path(metadata_dir, paste0(samp, "_orig_meta.rds"))

  if (!file.exists(orig_meta_file)) {
    cat("  _orig_meta.rds not found — skipping\n")
    skipped <- skipped + 1
    next
  }

  seu <- readRDS(f)
  orig_meta <- readRDS(orig_meta_file)

  # Match on cell barcodes (both use original GEO barcodes at this stage)
  common_cells <- intersect(colnames(seu), rownames(orig_meta))
  n_matched <- length(common_cells)
  n_total   <- ncol(seu)

  if (n_matched == 0) {
    cat("  WARNING: no barcode overlap — skipping (check if barcodes changed)\n")
    skipped <- skipped + 1
    rm(seu, orig_meta)
    gc()
    next
  }

  if (n_matched < n_total) {
    cat("  WARNING:", n_total - n_matched, "cells have no match in orig_meta",
        "(expected if CreateSeuratObject min.cells/min.features filtered some)\n")
  }

  # Overwrite percent.mt with real values; unmatched cells keep 0
  mt_patch <- rep(0, n_total)
  names(mt_patch) <- colnames(seu)
  mt_patch[common_cells] <- orig_meta[common_cells, "percent.mt"]
  seu$percent.mt <- mt_patch

  cat("  percent.mt patched:", n_matched, "/", n_total, "cells\n")
  cat("  Range:", round(min(seu$percent.mt), 2), "–",
      round(max(seu$percent.mt), 2), "%\n")
  cat("  Mean:", round(mean(seu$percent.mt), 2), "%\n")

  saveRDS(seu, f)
  cat("  Saved:", f, "\n")

  patched <- patched + 1
  rm(seu, orig_meta)
  gc()
}

cat("\n=============================\n")
cat("Done.\n")
cat("Patched:", patched, "samples\n")
cat("Skipped:", skipped, "samples\n")
