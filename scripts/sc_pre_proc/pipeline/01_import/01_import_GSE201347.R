#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE201347.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE201347
# Description: Import pre-built Seurat object (~18GB), extract count matrices,
#              convert gene symbols to Ensembl IDs, split by sample,
#              save individual sparse RDS files
# Input: data/sc_reference/raw/GSE201347/GSE201347.rds
# Output: data/sc_reference/raw/GSE201347/GSE201347_FIXED/<sample>.rds
# Note: GEO metadata is downloaded separately (01_download_metadata_GSE201347.R)
# Usage: Rscript 01_import_GSE201347.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# Ensembl version: 109 (GRCh38)
# Note: Runs on server due to memory requirements
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE201347.R <config_path>\n  No config file provided.")
}
if (!file.exists(args[1])) {
  stop("Config file not found: ", args[1])
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load functions ----
source(utils_path)

# ---- Libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ---- Load Seurat object ----
cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading Seurat object:", basename(seurat_rds_path), "\n")
cat("File size:", round(file.size(seurat_rds_path) / 1e9, 1), "GB\n")

seu <- readRDS(seurat_rds_path)

cat("Object dimensions:", dim(seu), "\n")

# ---- Identify samples ----
if (!sample_column %in% colnames(seu@meta.data)) {
  cat("Available metadata columns:", paste(colnames(seu@meta.data), collapse = ", "), "\n")
  stop("Sample column '", sample_column, "' not found in metadata. Update config.")
}

samples <- unique(seu@meta.data[[sample_column]])
cat("Number of samples:", length(samples), "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")

# ---- Extract count matrix ----
counts <- GetAssayData(seu, slot = "counts")
cat("Count matrix dimensions:", dim(counts), "\n")

# ---- Check gene ID type and convert if needed ----
if (is_ensembl(rownames(counts))) {
  cat("-> Ensembl IDs detected\n")
  counts <- clean_ensembl_ids(counts)
} else {
  cat("-> Gene symbols detected, converting to Ensembl\n")
  mapping_table <- readRDS(mapping_table_path)
  counts <- convert_to_ensembl(counts, mapping_table)
  rm(mapping_table)
}

cat("Count matrix after conversion:", dim(counts), "\n")

# ---- Free Seurat object, keep only metadata ----
cell_metadata <- seu@meta.data
rm(seu)
gc()

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Split by sample and save ----
for (s in samples) {

  cat("\n=============================\n")
  cat("Sample:", s, "\n")

  cells <- rownames(cell_metadata)[cell_metadata[[sample_column]] == s]
  # Keep only cells present in the count matrix
  cells <- intersect(cells, colnames(counts))
  mat <- counts[, cells]

  cat("Dimension:", dim(mat), "\n")
  cat("Final dimension (Ensembl):", dim(mat), "\n")

  outfile <- file.path(input_dir, paste0(s, ".rds"))
  saveRDS(mat, outfile)
  cat("Saved:", outfile, "\n")

  rm(mat)
  gc()
}

rm(counts, cell_metadata)
gc()

cat("\n=============================\n")
cat("Done! All samples saved to:", input_dir, "\n")
