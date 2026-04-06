#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 05_merge.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Merge filtered Seurat objects into a single object
#              and join layers for downstream analysis
# Input: Filtered Seurat objects (from script 04)
# Output: Merged Seurat object (.rds)
# Usage: Rscript 05_merge.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 05_merge.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- List filtered Seurat objects ----
seurat_files <- list.files(
  seurat_qc_filtered_dir,
  pattern = "_seurat_filtered.rds$",
  full.names = TRUE
)

cat("\n=============================\n")
cat("Script: 05_merge.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Filtered objects found:", length(seurat_files), "\n")

if (length(seurat_files) == 0) {
  stop("No filtered Seurat objects found in: ", seurat_qc_filtered_dir)
}

# ---- Exclude samples (from config) ----
sample_names_in_files <- sub("_seurat_filtered.rds", "", basename(seurat_files))

# Biological exclusion
if (exists("samples_to_exclude_bio") && length(samples_to_exclude_bio) > 0) {
  matched_bio <- samples_to_exclude_bio[samples_to_exclude_bio %in% sample_names_in_files]
  if (length(matched_bio) > 0) {
    cat("Excluding samples (biological):", paste(matched_bio, collapse = ", "), "\n")
    seurat_files <- seurat_files[!sample_names_in_files %in% matched_bio]
    sample_names_in_files <- sub("_seurat_filtered.rds", "", basename(seurat_files))
  }
}

# Quality exclusion
if (exists("samples_to_exclude_qc") && length(samples_to_exclude_qc) > 0) {
  matched_qc <- samples_to_exclude_qc[samples_to_exclude_qc %in% sample_names_in_files]
  if (length(matched_qc) > 0) {
    cat("Excluding samples (quality):", paste(matched_qc, collapse = ", "), "\n")
    seurat_files <- seurat_files[!sample_names_in_files %in% matched_qc]
  }
}

cat("Samples to merge:", length(seurat_files), "\n")

# ---- Load all objects ----
cat("\nLoading objects...\n")

seu_list <- list()

for (file in seurat_files) {

  samp <- sub("_seurat_filtered.rds", "", basename(file))
  cat("  Loading:", samp, "\n")

  seu_list[[samp]] <- readRDS(file)
}

cat("All objects loaded:", length(seu_list), "samples\n")

# ---- Print cell counts per sample ----
cat("\nCells per sample:\n")
for (s in names(seu_list)) {
  cat(sprintf("  %-30s %d cells\n", s, ncol(seu_list[[s]])))
}

total_cells <- sum(sapply(seu_list, ncol))
cat("  Total:", total_cells, "cells\n")

# ---- Merge ----
cat("\nMerging objects...\n")

if (length(seu_list) == 1) {
  merged_seurat <- seu_list[[1]]
} else {
  merged_seurat <- Reduce(function(x, y) merge(x, y), seu_list)
}

cat("Merged dimensions:", nrow(merged_seurat), "genes x", ncol(merged_seurat), "cells\n")

# ---- Free individual objects ----
rm(seu_list)
gc()

# ---- Join layers ----
cat("Joining layers...\n")
merged_seurat <- JoinLayers(merged_seurat)
cat("Layers joined successfully\n")

# ---- Rebuild gene symbol mapping (lost during merge) ----
cat("\nRebuilding gene symbol mapping...\n")
mapping_table <- readRDS(mapping_table_path)
mapping_table$ensembl_gene_id <- sub("[.].*", "", mapping_table$ensembl_gene_id)

ens_to_symbol <- setNames(
  mapping_table$hgnc_symbol,
  mapping_table$ensembl_gene_id
)

# Supplement with known immunoglobulin genes missing from Ensembl mapping
ig_supplement <- c(
  "ENSG00000211890" = "IGHG3",   "ENSG00000211897" = "IGHG1",
  "ENSG00000211892" = "IGHG4",   "ENSG00000211896" = "IGHG2",
  "ENSG00000211899" = "IGHM",    "ENSG00000211893" = "IGHA2",
  "ENSG00000211895" = "IGHA1",   "ENSG00000211891" = "IGHD",
  "ENSG00000211898" = "IGHE",    "ENSG00000211651" = "IGKV1-5",
  "ENSG00000239951" = "JCHAIN",  "ENSG00000275385" = "IGHGP",
  "ENSG00000211677" = "IGLC2",   "ENSG00000211679" = "IGLC3",
  "ENSG00000211662" = "IGLV3-21","ENSG00000211598" = "IGKV4-1",
  "ENSG00000254709" = "IGLL5"
)

# Only add if not already in mapping
new_ig <- ig_supplement[!names(ig_supplement) %in% names(ens_to_symbol) |
                         ens_to_symbol[names(ig_supplement)] == ""]
ens_to_symbol[names(new_ig)] <- new_ig

gene_symbols <- ens_to_symbol[rownames(merged_seurat)]
gene_symbols[is.na(gene_symbols)] <- ""
merged_seurat[["RNA"]]@meta.data$gene_symbol <- gene_symbols

n_mapped <- sum(gene_symbols != "")
n_ig <- sum(rownames(merged_seurat) %in% names(ig_supplement))
cat("Gene symbols mapped:", n_mapped, "/", nrow(merged_seurat), "\n")
cat("  (including", n_ig, "supplemented IG genes)\n")
cat("Unmapped (Ensembl only):", nrow(merged_seurat) - n_mapped, "\n")

rm(mapping_table, ens_to_symbol, ig_supplement, new_ig)

# ---- Verify metadata ----
cat("\nMetadata columns:", paste(colnames(merged_seurat@meta.data), collapse = ", "), "\n")
cat("Samples in merged object:", length(unique(merged_seurat$orig.ident)), "\n")

if ("dataset_id" %in% colnames(merged_seurat@meta.data)) {
  cat("Dataset IDs:", paste(unique(merged_seurat$dataset_id), collapse = ", "), "\n")
}

# ---- Save ----
dir.create(dirname(merged_seurat_path), recursive = TRUE, showWarnings = FALSE)

cat("\nSaving merged object...\n")
saveRDS(merged_seurat, merged_seurat_path)
cat("Saved:", merged_seurat_path, "\n")

cat(sprintf("File size: %.1f MB\n", file.size(merged_seurat_path) / 1e6))

cat("\n=============================\n")
cat("Done! Merged object:", ncol(merged_seurat), "cells,", nrow(merged_seurat), "genes\n")

rm(merged_seurat)
gc()
