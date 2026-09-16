#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE308231.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE308231
# Description: Read per-sample CellRanger filtered 10x files (prefixed triplets
#              in a single RAW directory), keep only gastric carcinoma in situ
#              samples (Ca_1–3), save as sparse RDS with Ensembl IDs as rownames.
# Input:  data/sc_reference/raw/GSE308231/GSE308231_RAW/GSM*_Ca_*_matrix.mtx.gz
#         (+ matching barcodes and features files)
# Output: data/sc_reference/raw/GSE308231/GSE308231_FIXED/<SAMPLE>.rds
#         data/sc_reference/metadata/metadata_GSE308231.rds
# Usage:  Rscript 01_import_GSE308231.R <config_path>
# Notes:  Gene IDs are already Ensembl (col1 of features file) — no conversion.
#         Peritoneal metastasis samples (F_14, F_15, F_16) are excluded.
# Author: Juliana Pinto
# Date: 2026-09-16
# Ensembl version: 109 (GRCh38)
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE308231.R <config_path>\n  No config file provided.")
}
if (!file.exists(args[1])) {
  stop("Config file not found: ", args[1])
}
source(normalizePath(args[1]))

# ---- Libraries ----
library(Matrix)
library(GEOquery)

cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Input dir:", raw_data_dir, "\n")

# ---- Find Ca_ sample matrix files (exclude peritoneal metastases F_) ----
all_matrix_files <- list.files(raw_data_dir,
                               pattern = "matrix\\.mtx\\.gz$",
                               full.names = TRUE)
matrix_files <- all_matrix_files[grepl("_Ca_[0-9]+_matrix\\.mtx\\.gz$",
                                       all_matrix_files)]

if (length(matrix_files) == 0) {
  stop("No Ca_ matrix files found in: ", raw_data_dir)
}

cat("Samples to import (Ca_ only):", length(matrix_files), "\n")
cat(paste("-", basename(matrix_files)), sep = "\n")
cat("\n")

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Import each sample ----
for (mat_f in matrix_files) {

  # Extract sample name from filename: GSMid_Ca_N_matrix.mtx.gz -> Ca_N
  prefix  <- sub("_matrix\\.mtx\\.gz$", "", basename(mat_f))
  samp_code <- sub("^[^_]+_", "", prefix)        # strip leading GSMid_
  samp_id   <- paste0(dataset_id, "_", samp_code)

  dir_f  <- dirname(mat_f)
  bc_f   <- file.path(dir_f, paste0(prefix, "_barcodes.tsv.gz"))
  feat_f <- file.path(dir_f, paste0(prefix, "_features.tsv.gz"))

  cat("=============================\n")
  cat("Sample:", samp_id, "\n")

  # ---- Read barcodes ----
  barcodes <- read.table(gzfile(bc_f), header = FALSE,
                         stringsAsFactors = FALSE)[, 1]
  cat("  Barcodes:", length(barcodes), "\n")

  # ---- Read features ----
  features <- read.table(gzfile(feat_f), header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
  if (ncol(features) >= 3) {
    colnames(features)[1:3] <- c("ensembl_id", "gene_symbol", "feature_type")
    gene_rows <- which(features$feature_type == "Gene Expression")
  } else {
    colnames(features)[1:2] <- c("ensembl_id", "gene_symbol")
    gene_rows <- seq_len(nrow(features))
  }
  cat("  Features total:", nrow(features),
      "| Gene Expression:", length(gene_rows), "\n")

  # ---- Read matrix ----
  mat <- readMM(gzfile(mat_f))   # rows = features, cols = barcodes
  mat <- mat[gene_rows, , drop = FALSE]
  rownames(mat) <- features$ensembl_id[gene_rows]
  colnames(mat) <- barcodes
  cat("  Matrix dimensions:", nrow(mat), "genes x", ncol(mat), "cells\n")

  # ---- Save ----
  out_file <- file.path(input_dir, paste0(samp_id, ".rds"))
  saveRDS(mat, out_file)
  cat("  Saved:", out_file, "\n")

  rm(mat, barcodes, features)
  gc()
}

# ---- Download GEO metadata ----
cat("\n=============================\n")
cat("Downloading GEO metadata...\n")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
gse      <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse[[1]])
saveRDS(metadata, metadata_path)
cat("Metadata saved:", metadata_path, "\n")

cat("\n=============================\n")
cat("Import complete.\n")
cat("Output directory:", input_dir, "\n")
cat("Files saved:", length(matrix_files), "samples\n")
