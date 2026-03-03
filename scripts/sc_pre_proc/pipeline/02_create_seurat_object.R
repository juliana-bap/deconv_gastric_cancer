#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 02_create_seurat_object.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Create Seurat objects from fixed CSV matrices
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
message("Number of CSV files found: ", length(csvs))
print(csvs)

# ---- Read mapping ----
mapping_table <- readRDS(mapping_table_path)

# remove version if exists
mapping_table$ensembl_gene_id <- sub("\\..*", "", mapping_table$ensembl_gene_id)

# Create named vector: Ensembl → Symbol
ens_to_symbol <- setNames(
  mapping_table$hgnc_symbol,
  mapping_table$ensembl_gene_id
)

for (path in csvs) {
  
  samp <- sub("_.*", "", basename(path))
  message("Processing: ", samp)
  
  df <- read.csv(path, check.names = FALSE)
  
  rownames(df) <- make.unique(df[[1]])
  df[[1]] <- NULL
  
  counts_mat <- as.matrix(df)
  colnames(counts_mat) <- make.unique(colnames(counts_mat))
  
  seu <- CreateSeuratObject(
    counts = counts_mat,
    project = samp,
    min.cells = min_cells,
    min.features = min_features
  )
  
  # ---- Add gene symbol as feature metadata ----
  
  # Match in correct order
  gene_symbols <- ens_to_symbol[rownames(seu)]
  
  # Replace NA with empty string
  gene_symbols[is.na(gene_symbols)] <- ""
  
  # Add as feature metadata
  seu[["RNA"]]@meta.data$gene_symbol <- gene_symbols
  
  message("Gene symbols added as feature metadata.")

message("Gene symbols added as feature metadata.")
  
  saveRDS(
    seu,
    file = file.path(seurat_raw_dir, paste0(samp, "_seurat_raw.rds"))
  )
  
  rm(seu, df, counts_mat)
  gc()
}

rm(list = ls())
gc()