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

message("Working directory: ", getwd())
message("Input dir: ", input_dir)
message("Exists? ", dir.exists(input_dir))

# ---- Create output directory ----
dir.create(seurat_raw_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Loop over samples ----
message("Number of RDS files found: ", length(rds_files))
print(rds_files)

# ---- Read mapping ----
mapping_table <- readRDS(mapping_table_path)

# remove version if exists
mapping_table$ensembl_gene_id <- sub("\\..*", "", mapping_table$ensembl_gene_id)

# Create named vector: Ensembl -> Symbol
ens_to_symbol <- setNames(
  mapping_table$hgnc_symbol,
  mapping_table$ensembl_gene_id
)

for (path in rds_files) {

  samp <- sub("_.*", "", basename(path))
  message("Processing: ", samp)

  # Read sparse matrix from RDS
  counts_mat <- readRDS(path)

  seu <- CreateSeuratObject(
    counts = counts_mat,
    project = samp,
    min.cells = min_cells,
    min.features = min_features
  )

  # ---- Add gene symbol as feature metadata ----
  gene_symbols <- ens_to_symbol[rownames(seu)]
  gene_symbols[is.na(gene_symbols)] <- ""
  seu[["RNA"]]@meta.data$gene_symbol <- gene_symbols

  message("Gene symbols added as feature metadata.")

  saveRDS(
    seu,
    file = file.path(seurat_raw_dir, paste0(samp, "_seurat_raw.rds"))
  )

  rm(seu, counts_mat)
  gc()
}

rm(list = ls())
gc()
