#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE249874.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE249874
# Description: Read CellRanger aggr raw (unfiltered) MTX, remove empty droplets
#              by min UMI pre-filter, split by barcode suffix into per-sample
#              sparse matrices, save as RDS with Ensembl IDs as rownames.
# Input:  data/sc_reference/raw/GSE249874/GSE249874_raw_feature_*.tsv.gz
#         data/sc_reference/raw/GSE249874/GSE249874_raw_feature_matrix.mtx.gz
# Output: data/sc_reference/raw/GSE249874/GSE249874_FIXED/<SAMPLE>.rds
#         data/sc_reference/metadata/metadata_GSE249874.rds
# Usage:  Rscript 01_import_GSE249874.R <config_path>
# Note:   Reads ~120M barcodes into memory. Requires 8–16 GB RAM.
#         Recommend running on HPC for comfortable execution.
# Author: Juliana Pinto
# Date: 2026-09-14
# Ensembl version: 109 (GRCh38)
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE249874.R <config_path>\n  No config file provided.")
}
if (!file.exists(args[1])) {
  stop("Config file not found: ", args[1])
}
source(normalizePath(args[1]))

# ---- Libraries ----
library(Matrix)
library(GEOquery)

# ---- Sample name mapping (barcode suffix → clean sample ID) ----
# Barcode suffix -N corresponds to sampleN in the CellRanger aggr output.
# Sample groups:
#   GC-HP-N (1-3): H. pylori-negative gastric adenocarcinoma
#   GC-HP-P (4-6): H. pylori-positive gastric adenocarcinoma
#   GS-HP-N (7-9): H. pylori-negative non-atrophic gastritis
#   GS-HP-P (10-12): H. pylori-positive non-atrophic gastritis
#   IM-HP-N (13-15): H. pylori-negative intestinal metaplasia
#   IM-HP-P (16-18): H. pylori-positive intestinal metaplasia
sample_names <- c(
  "1"  = "GC-HP-N_1",
  "2"  = "GC-HP-N_2",
  "3"  = "GC-HP-N_3",
  "4"  = "GC-HP-P_1",
  "5"  = "GC-HP-P_2",
  "6"  = "GC-HP-P_3",
  "7"  = "GS-HP-N_1",
  "8"  = "GS-HP-N_2",
  "9"  = "GS-HP-N_3",
  "10" = "GS-HP-P_1",
  "11" = "GS-HP-P_2",
  "12" = "GS-HP-P_3",
  "13" = "IM-HP-N_1",
  "14" = "IM-HP-N_2",
  "15" = "IM-HP-N_3",
  "16" = "IM-HP-P_1",
  "17" = "IM-HP-P_2",
  "18" = "IM-HP-P_3"
)

cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Format:  CellRanger aggr (raw/unfiltered), 18 samples\n")
cat("Pre-filter min counts:", prefilter_min_counts, "\n")

# ---- Read features ----
cat("\nReading features...\n")
features <- read.table(
  gzfile(feature_file), header = FALSE, sep = "\t",
  stringsAsFactors = FALSE,
  col.names = c("ensembl_id", "gene_symbol", "feature_type")
)
gene_features <- features[features$feature_type == "Gene Expression", ]
cat("  Features total:", nrow(features), "| Gene Expression:", nrow(gene_features), "\n")

# ---- Read barcodes ----
cat("Reading barcodes (360 MB compressed — may take a few minutes)...\n")
barcodes <- read.table(gzfile(barcode_file), header = FALSE,
                       stringsAsFactors = FALSE)[, 1]
cat("  Total barcodes:", length(barcodes), "\n")

# Extract sample suffix (-1 to -18 → "1" to "18")
barcode_suffix <- sub(".*-([0-9]+)$", "\\1", barcodes)

# ---- Read MTX matrix ----
cat("Reading MTX matrix (1.3 GB compressed — may take several minutes)...\n")
mat <- readMM(gzfile(matrix_file))   # rows = features, cols = barcodes
rownames(mat) <- features$ensembl_id
colnames(mat) <- barcodes
cat("  Raw matrix dimensions:", nrow(mat), "features x", ncol(mat), "barcodes\n")

# Subset to Gene Expression features only
mat <- mat[gene_features$ensembl_id, , drop = FALSE]
cat("  After Gene Expression filter:", nrow(mat), "genes x", ncol(mat), "barcodes\n")

# ---- Pre-filter: remove empty droplets ----
cat("\nPre-filtering barcodes (min counts =", prefilter_min_counts, ")...\n")
cell_counts    <- Matrix::colSums(mat)
keep           <- cell_counts >= prefilter_min_counts
mat            <- mat[, keep, drop = FALSE]
barcode_suffix <- barcode_suffix[keep]
cat("  Retained:", sum(keep), "/", length(keep), "barcodes\n")

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Split by sample and save ----
cat("\nSplitting by sample and saving...\n")

for (sfx in names(sample_names)) {
  samp_code <- sample_names[[sfx]]
  samp_id   <- paste0(dataset_id, "_", samp_code)

  cells_i <- which(barcode_suffix == sfx)

  if (length(cells_i) == 0) {
    cat("  WARNING: no cells for suffix -", sfx, "(", samp_id, ")\n")
    next
  }

  mat_i <- mat[, cells_i, drop = FALSE]

  # Normalise barcode suffix to -1 (convention for single-sample objects)
  colnames(mat_i) <- sub("-[0-9]+$", "-1", colnames(mat_i))

  out_file <- file.path(input_dir, paste0(samp_id, ".rds"))
  saveRDS(mat_i, out_file)
  cat("  Saved:", samp_id, "->", ncol(mat_i), "cells x", nrow(mat_i), "genes\n")

  rm(mat_i)
  gc()
}

# ---- Download GEO metadata ----
cat("\nDownloading GEO metadata...\n")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
gse      <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse[[1]])
saveRDS(metadata, metadata_path)
cat("Metadata saved:", metadata_path, "\n")

cat("\n=============================\n")
cat("Import complete.\n")
cat("Output directory:", input_dir, "\n")
cat("Files saved:", length(sample_names), "samples\n")
