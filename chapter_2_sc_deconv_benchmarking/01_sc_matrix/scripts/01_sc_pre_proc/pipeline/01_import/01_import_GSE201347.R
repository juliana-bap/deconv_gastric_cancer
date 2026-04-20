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
# NOTE: Do NOT load Seurat before readRDS.
# The RDS may have been saved with a different Seurat version,
# causing S4 class validation errors on load.
# Loading only Matrix allows R to read the object as a generic S4,
# and we extract counts/metadata directly via @ slots.
suppressPackageStartupMessages({
  library(Matrix)
})

# ---- Load Seurat object (without Seurat library) ----
cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading Seurat object:", basename(seurat_rds_path), "\n")
cat("File size:", round(file.size(seurat_rds_path) / 1e9, 1), "GB\n")

seu <- readRDS(seurat_rds_path)

cat("Object loaded successfully\n")

# ---- Extract count matrix (handle v3/v4 and v5 structures) ----
counts <- tryCatch(
  seu@assays$RNA@counts,
  error = function(e1) {
    tryCatch(
      seu@assays$RNA@layers$counts,
      error = function(e2) {
        stop("Could not extract count matrix. Tried @counts (v3/v4) and @layers$counts (v5).")
      }
    )
  }
)

cat("Count matrix dimensions:", nrow(counts), "x", ncol(counts), "\n")

# ---- Extract cell metadata ----
cell_metadata <- seu@meta.data

cat("Available metadata columns:", paste(colnames(cell_metadata), collapse = ", "), "\n")

# ---- Identify samples ----
if (!sample_column %in% colnames(cell_metadata)) {
  stop("Sample column '", sample_column, "' not found in metadata. Update config.")
}

samples <- unique(cell_metadata[[sample_column]])
cat("Number of samples:", length(samples), "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")

# ---- Free Seurat object early (keep only counts + metadata) ----
rm(seu)
gc()

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
