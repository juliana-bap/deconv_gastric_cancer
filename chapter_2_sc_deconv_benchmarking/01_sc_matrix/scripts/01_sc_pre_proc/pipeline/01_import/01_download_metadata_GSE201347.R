#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_download_metadata_GSE201347.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE201347
# Description: Download GEO metadata (runs locally, requires internet)
# Output: data/sc_reference/metadata/metadata_GSE201347.rds
# Usage: Rscript 01_download_metadata_GSE201347.R
# Author: Juliana Pinto
# Date: 2026
# Note: Run locally, then transfer metadata_GSE201347.rds to the server
# ==============================

library(GEOquery)

# ---- Paths (local) ----
# here::here() detects the repo root automatically — no need to edit
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' required. Run scripts/setup/install_r_deps.R first.")
}
project_root <- here::here()
dataset_id <- "GSE201347"
metadata_dir <- file.path(project_root, "data", "sc_reference", "metadata")
metadata_path <- file.path(metadata_dir, paste0("metadata_", dataset_id, ".rds"))

# ---- Download ----
cat("\n=============================\n")
cat("Downloading GEO metadata for", dataset_id, "...\n")

dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse[[1]])

cat("Samples in metadata:", nrow(metadata), "\n")
cat("Columns:", paste(colnames(metadata), collapse = ", "), "\n")

saveRDS(metadata, metadata_path)
cat("Metadata saved:", metadata_path, "\n")
cat("Transfer to server: data/sc_reference/metadata/metadata_GSE201347.rds\n")
