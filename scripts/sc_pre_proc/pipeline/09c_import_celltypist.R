#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 09c_import_celltypist.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Import CellTypist annotations (CSV from 09b) into the
#              Seurat object. Exploration and comparison with SingleR
#              is done in the annotation notebook.
# Input: Annotated Seurat object (from 09a) + CellTypist CSV (from 09b)
# Output: Seurat object with all annotations (.rds)
# Usage: Rscript 09c_import_celltypist.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 09c_import_celltypist.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load annotated Seurat object ----
cat("\n=============================\n")
cat("Script: 09c_import_celltypist.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading:", basename(annotated_seurat_path), "\n")

seu <- readRDS(annotated_seurat_path)

cat("Dimensions:", nrow(seu), "genes x", ncol(seu), "cells\n")
cat("Existing annotations:", paste(grep("singler|celltypist", colnames(seu@meta.data), value = TRUE), collapse = ", "), "\n")

# ---- Load CellTypist CSV ----
celltypist_csv <- file.path(annotation_output_dir, paste0(dataset_id, "_celltypist_annotations.csv"))

if (!file.exists(celltypist_csv)) {
  stop("CellTypist CSV not found: ", celltypist_csv, "\n  Run 09b first.")
}

cat("\n=============================\n")
cat("Loading CellTypist annotations:", basename(celltypist_csv), "\n")

ct_df <- read.csv(celltypist_csv, row.names = "cell_id", stringsAsFactors = FALSE)

# Match cells
common_cells <- intersect(colnames(seu), rownames(ct_df))
cat("Cells matched:", length(common_cells), "/", ncol(seu), "\n")

if (length(common_cells) == 0) {
  stop("No matching cell IDs. Check barcode formats between Seurat and h5ad.")
}

# Add to Seurat
seu$celltypist_label <- ct_df[colnames(seu), "celltypist_label"]
seu$celltypist_conf_score <- as.numeric(ct_df[colnames(seu), "celltypist_conf_score"])

cat("CellTypist labels:", length(unique(na.omit(seu$celltypist_label))), "cell types\n")
cat("Mean confidence:", round(mean(seu$celltypist_conf_score, na.rm = TRUE), 3), "\n")

# ---- Print summary ----
cat("\n=============================\n")
cat("CellTypist label counts:\n")
print(sort(table(seu$celltypist_label), decreasing = TRUE))

# ---- Save ----
cat("\n=============================\n")
cat("Saving object with all annotations...\n")

saveRDS(seu, annotated_seurat_path)
cat("Saved:", annotated_seurat_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(annotated_seurat_path) / 1e6))

cat("\n=============================\n")
cat("Done! Metadata columns now include:\n")
cat("  singler_hpca\n")
cat("  singler_blueprint\n")
cat("  celltypist_label\n")
cat("  celltypist_conf_score\n")
cat("Next: explore and compare annotations in the notebook\n")

rm(seu)
gc()
