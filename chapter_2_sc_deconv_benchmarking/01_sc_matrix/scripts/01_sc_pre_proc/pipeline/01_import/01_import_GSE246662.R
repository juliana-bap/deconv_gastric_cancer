#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_and_fix_GSE246662.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE246662
# Description: Import CSV matrices, fix orientation, convert gene symbols
#              to Ensembl IDs, save as sparse RDS, and download GEO metadata
# Input: data/sc_reference/raw/GSE246662/*.csv.gz
# Output: data/sc_reference/raw/GSE246662/GSE246662_FIXED/*.rds
#         data/sc_reference/metadata/metadata_GSE246662.rds
# Usage: Rscript 01_import_and_fix_GSE246662.R <config_path>
# Author: Juliana Pinto
# Date: 23-02-2026
# Ensembl version: 109 (GRCh38)
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_and_fix_GSE246662.R <config_path>\n  No config file provided.")
}
if (!file.exists(args[1])) {
  stop("Config file not found: ", args[1])
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load functions ----
source(utils_path)

# ---- Libraries ----
library(Matrix)
library(GEOquery)

# ---- Mapping table ----
mapping_table <- readRDS(mapping_table_path)

# ---- Input files ----
csvs <- list.files(raw_base_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Number of files found:", length(csvs), "\n")

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Import and fix count matrices ----
for (path in csvs) {

  cat("\n=============================\n")
  cat("Sample:", basename(path), "\n")

  df <- read.csv(path, row.names = 1, check.names = FALSE)

  # Fix orientation (genes in rows, cells in columns)
  df <- fix_orientation(df)
  cat("Dimension after orientation:", dim(df), "\n")

  # Convert gene symbols to Ensembl IDs
  df <- convert_to_ensembl(df, mapping_table)
  cat("Final dimension (Ensembl):", dim(df), "\n")

  # Save as sparse RDS
  outfile <- file.path(input_dir, sub("\\.csv\\.gz$", ".rds", basename(path)))
  saveRDS(df, outfile)
  cat("Saved:", outfile, "\n")

  rm(df)
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
