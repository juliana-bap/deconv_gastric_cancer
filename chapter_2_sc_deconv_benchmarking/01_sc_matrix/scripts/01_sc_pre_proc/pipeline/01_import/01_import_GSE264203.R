#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE264203.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE264203
# Description: Import HDF5 10x file, convert gene symbols to Ensembl IDs,
#              save as sparse RDS, and download GEO metadata
# Input: data/sc_reference/raw/GSE264203/GSE264203.h5
# Output: data/sc_reference/raw/GSE264203/GSE264203_FIXED/<sample>.rds (one per barcode suffix)
#         data/sc_reference/metadata/metadata_GSE264203.rds
# Usage: Rscript 01_import_GSE264203.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# Note: Single H5 file (not split by sample)
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE264203.R <config_path>\n  No config file provided.")
}
if (!file.exists(args[1])) {
  stop("Config file not found: ", args[1])
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load functions ----
source(utils_path)

# ---- Libraries ----
library(Seurat)
library(Matrix)
library(GEOquery)

# ---- Mapping table ----
mapping_table <- readRDS(mapping_table_path)

# ---- Input ----
h5_path <- file.path(raw_base_dir, paste0(dataset_id, ".h5"))

cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading H5 file:", basename(h5_path), "\n")
cat("File size:", round(file.size(h5_path) / 1e6, 1), "MB\n")

# ---- Read H5 ----
cat("\n=============================\n")
cat("Sample:", basename(h5_path), "\n")

counts <- Read10X_h5(filename = h5_path)

# Handle multiple modalities (keep Gene Expression only)
if (is.list(counts)) {
  cat("-> Multiple modalities detected:", paste(names(counts), collapse = ", "), "\n")
  cat("-> Using 'Gene Expression' matrix\n")
  counts <- counts[["Gene Expression"]]
}

cat("Dimension after import:", dim(counts), "\n")

# ---- Check gene ID type and convert if needed ----
if (is_ensembl(rownames(counts))) {
  cat("-> Ensembl IDs detected\n")
  counts <- clean_ensembl_ids(counts)
} else {
  cat("-> Gene symbols detected, converting to Ensembl\n")
  counts <- convert_to_ensembl(counts, mapping_table)
}

cat("Total dimension (Ensembl):", dim(counts), "\n")

# ---- Split by barcode suffix ----
suffixes <- sub(".*-", "", colnames(counts))
unique_suffixes <- sort(unique(suffixes))
cat("Barcode suffixes found:", length(unique_suffixes), "\n")

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Map suffixes to sample names (from config, if available) ----
# suffix_to_sample allows renaming: c("1"="GC1", "2"="GC2", ...)
# If NULL, uses generic names: GSE264203_sample_1, GSE264203_sample_2, ...

for (s in unique_suffixes) {

  cat("\n=============================\n")

  if (!is.null(suffix_to_sample) && s %in% names(suffix_to_sample)) {
    sample_name <- suffix_to_sample[[s]]
  } else {
    sample_name <- paste0(dataset_id, "_sample_", s)
  }

  cat("Suffix -", s, "->", sample_name, "\n")

  cells <- colnames(counts)[suffixes == s]
  mat <- counts[, cells]

  cat("Dimension:", dim(mat), "\n")
  cat("Final dimension (Ensembl):", dim(mat), "\n")

  outfile <- file.path(input_dir, paste0(sample_name, ".rds"))
  saveRDS(mat, outfile)
  cat("Saved:", outfile, "\n")

  rm(mat)
  gc()
}

rm(counts)
gc()

# ---- Download GEO metadata ----
cat("\n=============================\n")
cat("Downloading GEO metadata for", dataset_id, "...\n")

dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse[[1]])

saveRDS(metadata, metadata_path)
cat("Metadata saved:", metadata_path, "\n")
