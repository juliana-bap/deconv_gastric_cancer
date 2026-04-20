#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE291080.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE291080
# Description: Import 10x Genomics data, use Ensembl IDs from features.tsv,
#              save as sparse RDS, and download GEO metadata.
#              Excludes samples listed in config (data integrity issues).
# Input: data/sc_reference/raw/GSE291080/<sample_dirs>/
# Output: data/sc_reference/raw/GSE291080/GSE291080_FIXED/*.rds
#         data/sc_reference/metadata/metadata_GSE291080.rds
# Usage: Rscript 01_import_GSE291080.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE291080.R <config_path>\n  No config file provided.")
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

# ---- Input ----
sample_dirs <- list.dirs(raw_base_dir, recursive = FALSE, full.names = TRUE)

# Exclude the output directory (FIXED) from sample list
sample_dirs <- sample_dirs[basename(sample_dirs) != basename(input_dir)]

# Exclude samples with known data integrity issues (from config)
sample_dirs <- sample_dirs[!basename(sample_dirs) %in% samples_to_exclude_import]

cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Number of samples found:", length(sample_dirs), "\n")
cat("Excluded samples:", paste(samples_to_exclude_import, collapse = ", "), "\n")

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Import count matrices ----
failed_samples <- character()

for (d in sample_dirs) {

  sample_name <- basename(d)
  cat("\n=============================\n")
  cat("Sample:", sample_name, "\n")

  # Read with error handling (some samples may have issues)
  counts <- tryCatch(
    Read10X(data.dir = d, gene.column = 1),
    error = function(e) {
      cat("FAILED -", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(counts)) {
    failed_samples <- c(failed_samples, sample_name)
    next
  }

  cat("Dimension after import:", dim(counts), "\n")

  # Verify Ensembl IDs
  if (is_ensembl(rownames(counts))) {
    cat("-> Ensembl IDs detected\n")
  } else {
    cat("WARNING: Gene IDs do not appear to be Ensembl format\n")
  }

  # Clean Ensembl IDs (strip version, deduplicate)
  counts <- clean_ensembl_ids(counts)

  cat("Final dimension (Ensembl):", dim(counts), "\n")

  # Save as RDS
  outfile <- file.path(input_dir, paste0(sample_name, ".rds"))
  saveRDS(counts, outfile)
  cat("Saved:", outfile, "\n")

  rm(counts)
  gc()
}

# ---- Report ----
if (length(failed_samples) > 0) {
  cat("\n=============================\n")
  cat("WARNING: Failed samples during import:\n")
  cat(paste("-", failed_samples), sep = "\n")
  cat("\n")
}

# ---- Download GEO metadata ----
cat("\n=============================\n")
cat("Downloading GEO metadata for", dataset_id, "...\n")

dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse[[1]])

saveRDS(metadata, metadata_path)
cat("Metadata saved:", metadata_path, "\n")
