#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE275648.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE275648
# Description: Import per-sample Seurat objects (.rds.gz) from GEO RAW archive,
#              extract count matrices, convert gene symbols → Ensembl IDs using
#              the pipeline mapping table, and save as sparse RDS per sample.
# Input:  data/sc_reference/raw/GSE275648/GSE275648_RAW/<GSM>_<sample>.rds.gz
#           Each file is a gzip-compressed Seurat object (gene symbols as rownames).
# Output: data/sc_reference/raw/GSE275648/GSE275648_FIXED/<GSM>_<sample>.rds
#         data/sc_reference/metadata/metadata_GSE275648.rds
# Usage:  Rscript 01_import_GSE275648.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# Notes:  - Format differs from pilot (GSE163558): GEO RAW contains Seurat
#           objects with gene SYMBOLS as rownames, not 10x MTX with Ensembl IDs.
#         - Gene symbols are converted to Ensembl IDs via convert_to_ensembl()
#           (ensembl109 mapping table, same as the rest of the pipeline).
#         - Genes with no Ensembl match (~hgnc_symbol=="") are dropped — these
#           are typically pseudogenes or non-standard symbols.
#         - readRDS(gzfile(...)) fails for these files; decompression is done
#           via system("gunzip -c") to a tempfile before readRDS.
#         - All 11 samples are included (7 GC + 4 adjacent normal).
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE275648.R <config_path>\n  No config file provided.")
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

# ---- Input: list .rds.gz files ----
gz_files <- list.files(raw_base_dir, pattern = "\\.rds\\.gz$", full.names = TRUE)

cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Number of samples found:", length(gz_files), "\n")
cat("Samples:\n")
for (f in gz_files) cat(" -", basename(f), "\n")

if (length(gz_files) == 0) {
  stop("No .rds.gz files found in: ", raw_base_dir)
}

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load mapping table (Ensembl 109) ----
mapping_table <- readRDS(mapping_table_path)
mapping_table$ensembl_gene_id <- sub("\\..*", "", mapping_table$ensembl_gene_id)
cat("\nMapping table loaded:", nrow(mapping_table), "entries\n")

# ---- Import per sample ----
for (gz_path in gz_files) {

  sample_name <- sub("\\.rds\\.gz$", "", basename(gz_path))

  cat("\n=============================\n")
  cat("Sample:", sample_name, "\n")

  # readRDS(gzfile(...)) and pipe() both fail for these files.
  # Use system2() to decompress to a tempfile — this correctly redirects
  # gunzip stdout to a file (unlike system() with shell redirection).
  tmp <- tempfile(fileext = ".rds")
  system2("gunzip", args = c("-c", gz_path), stdout = tmp)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Decompression failed for: ", basename(gz_path))
  }
  obj <- readRDS(tmp)
  file.remove(tmp)

  cat("Object class:", class(obj), "\n")
  cat("Dimensions (genes x cells):", nrow(obj), "x", ncol(obj), "\n")

  # Extract count matrix (Seurat v5 API; falls back to v4 slot= if needed)
  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )

  # Ensure sparse (dgCMatrix)
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  }

  # ---- Filter cells by percent.mt (from original Seurat meta.data) ----
  # MT genes were removed from the count matrix by the authors before GEO
  # submission, so percent.mt cannot be recomputed downstream (script 03
  # gives 0 for all cells). Filtering is applied here at import using the
  # pre-computed values stored in the original Seurat meta.data.
  # Uses qc_max_percent_mt from the config (same variable as script 04),
  # so the threshold is controlled in one place.
  if ("percent.mt" %in% colnames(obj@meta.data)) {
    mt_vals <- obj@meta.data[colnames(counts), "percent.mt"]
    keep_cells <- !is.na(mt_vals) & mt_vals <= qc_max_percent_mt
    n_removed <- sum(!keep_cells)
    cat("MT filter (percent.mt <=", qc_max_percent_mt, "%):",
        n_removed, "cells removed,", sum(keep_cells), "kept\n")
    counts <- counts[, keep_cells, drop = FALSE]
  } else {
    cat("WARNING: percent.mt not in original meta.data — MT filter skipped\n")
  }

  cat("Dimensions after MT filter:", nrow(counts), "genes x", ncol(counts), "cells\n")
  cat("Gene IDs (first 5):", paste(head(rownames(counts), 5), collapse = ", "), "\n")

  # Detect if already Ensembl (precaution)
  if (is_ensembl(rownames(counts))) {
    cat("-> Ensembl IDs already present; running clean_ensembl_ids()\n")
    counts <- clean_ensembl_ids(counts)
  } else {
    cat("-> Gene symbols detected; converting to Ensembl via convert_to_ensembl()\n")
    counts <- convert_to_ensembl(counts, mapping_table)
  }

  cat("Final dimension:", nrow(counts), "genes x", ncol(counts), "cells\n")

  # Save count matrix
  outfile <- file.path(input_dir, paste0(sample_name, ".rds"))
  saveRDS(counts, outfile)
  cat("Saved:", outfile, "\n")

  # Save original cell metadata for documentation (filtered cells only).
  # percent.mt was pre-computed by authors before they removed MT genes from
  # the count matrix. Saved to metadata_dir (not input_dir) to avoid being
  # picked up by rds_files in the config — not consumed by downstream scripts.
  orig_meta_cols <- intersect("percent.mt", colnames(obj@meta.data))
  if (length(orig_meta_cols) > 0) {
    orig_meta <- obj@meta.data[colnames(counts), orig_meta_cols, drop = FALSE]
    meta_outfile <- file.path(metadata_dir, paste0(sample_name, "_orig_meta.rds"))
    saveRDS(orig_meta, meta_outfile)
    cat("percent.mt range (post-filter):", round(min(orig_meta$percent.mt), 2),
        "–", round(max(orig_meta$percent.mt), 2), "\n")
  }

  rm(obj, counts)
  gc()
}

# ---- Download GEO metadata ----
cat("\n=============================\n")
cat("Downloading GEO metadata for", dataset_id, "...\n")

dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse[[1]])

saveRDS(metadata, metadata_path)
cat("Metadata saved:", metadata_path, "\n")
cat("Metadata columns available:", paste(colnames(metadata), collapse = ", "), "\n")
