#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 02_create_seurat_object.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Create Seurat objects from sparse RDS matrices (output of script 01)
# Input: config file
# Output: Seurat .rds objects
# Author: Juliana Pinto
# Date: 23-02-2026
# ==============================

library(Seurat)

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Provide config file path. Example:
       Rscript 02_create_seurat_object.R config/config_GSE246662.R")
}

config_path <- normalizePath(args[1])
source(config_path)

cat("\n=============================\n")
cat("Script: 02_create_seurat_object.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Input dir:", input_dir, "\n")
cat("Output dir:", seurat_raw_dir, "\n")
cat("Number of RDS files found:", length(rds_files), "\n")

# ---- Create output directory ----
dir.create(seurat_raw_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read mapping table ----
mapping_table <- readRDS(mapping_table_path)

# Remove version if exists
mapping_table$ensembl_gene_id <- sub("\\..*", "", mapping_table$ensembl_gene_id)

# Create named vector: Ensembl -> Symbol
ens_to_symbol <- setNames(
  mapping_table$hgnc_symbol,
  mapping_table$ensembl_gene_id
)

# ---- Read GEO metadata ----
if (file.exists(metadata_path)) {
  geo_metadata <- readRDS(metadata_path)
  cat("GEO metadata loaded:", nrow(geo_metadata), "samples,", ncol(geo_metadata), "columns\n")
} else {
  geo_metadata <- NULL
  cat("WARNING: GEO metadata not found at:", metadata_path, "\n")
  cat("Continuing without sample-level metadata.\n")
}

total_samples <- length(rds_files)

for (i in seq_along(rds_files)) {

  path <- rds_files[i]
  samp <- sub("\\.rds$", "", basename(path))

  cat("\n=============================\n")
  cat(sprintf("Sample [%d/%d]: %s\n", i, total_samples, samp))

  # Read sparse matrix from RDS
  counts_mat <- readRDS(path)
  cat("Input matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

  seu <- CreateSeuratObject(
    counts = counts_mat,
    project = samp,
    min.cells = min_cells,
    min.features = min_features
  )

  cat("After filtering (min.cells=", min_cells, ", min.features=", min_features, "):",
      nrow(seu), "genes x", ncol(seu), "cells\n")

  # ---- Add gene symbol as feature metadata ----
  gene_symbols <- ens_to_symbol[rownames(seu)]
  n_mapped <- sum(!is.na(gene_symbols))
  gene_symbols[is.na(gene_symbols)] <- ""
  seu[["RNA"]]@meta.data$gene_symbol <- gene_symbols

  cat("Gene symbols mapped:", n_mapped, "/", nrow(seu), "\n")

  # ---- Add dataset and sample identifiers ----
  seu$dataset_id <- dataset_id

  # Extract GSM ID if filename starts with GSM, otherwise use full sample name
  if (grepl("^GSM", samp)) {
    gsm_id <- sub("_.*", "", samp)
  } else {
    gsm_id <- samp
  }
  seu$sample_id <- gsm_id

  cat("dataset_id:", dataset_id, "\n")
  cat("sample_id:", gsm_id, "\n")

  # ---- Add GEO metadata (sample-level -> all cells) ----
  if (!is.null(geo_metadata) && gsm_id %in% rownames(geo_metadata)) {

    sample_meta <- geo_metadata[gsm_id, ]

    # Select biologically relevant columns (pattern: ":ch1" suffix)
    bio_cols <- grep(":ch1$", colnames(sample_meta), value = TRUE)

    # Also include title and source_name
    extra_cols <- intersect(c("title", "source_name_ch1"), colnames(sample_meta))
    bio_cols <- unique(c(extra_cols, bio_cols))

    for (col in bio_cols) {
      # Clean column name: "tissue:ch1" -> "geo_tissue"
      clean_name <- sub(":ch1$", "", col)
      clean_name <- sub("_ch1$", "", clean_name)
      clean_name <- paste0("geo_", clean_name)

      seu[[clean_name]] <- as.character(sample_meta[[col]])
    }

    cat("GEO metadata added:", length(bio_cols), "columns\n")

  } else if (!is.null(geo_metadata)) {
    cat("WARNING: GSM", gsm_id, "not found in GEO metadata\n")
  }

  # ---- Save ----
  outfile <- file.path(seurat_raw_dir, paste0(samp, "_seurat_raw.rds"))
  saveRDS(seu, file = outfile)
  cat("Saved:", outfile, "\n")

  rm(seu, counts_mat)
  gc()
}

cat("\n=============================\n")
cat("Done! All", total_samples, "Seurat objects saved to:", seurat_raw_dir, "\n")

rm(list = ls())
gc()
