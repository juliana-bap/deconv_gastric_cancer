#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 06_normalize_hvg.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Normalize merged Seurat object and identify
#              Highly Variable Genes (HVGs). Saves top HVG plot.
# Input: Merged Seurat object (from script 05)
# Output: Normalized Seurat object with HVGs (.rds) + HVG plot (PDF)
# Usage: Rscript 06_normalize_hvg.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 06_normalize_hvg.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load merged object ----
cat("\n=============================\n")
cat("Script: 06_normalize_hvg.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading merged object:", basename(merged_seurat_path), "\n")

seu <- readRDS(merged_seurat_path)

cat("Dimensions:", nrow(seu), "genes x", ncol(seu), "cells\n")
cat("Samples:", length(unique(seu$orig.ident)), "\n")

# ---- Normalization ----
cat("\n=============================\n")
cat("Normalizing...\n")
cat("  Method:", norm_method, "\n")
cat("  Scale factor:", scale_factor, "\n")

seu <- NormalizeData(
  seu,
  normalization.method = norm_method,
  scale.factor = scale_factor
)

cat("Normalization complete\n")

# ---- Highly Variable Genes ----
cat("\n=============================\n")
cat("Finding variable features...\n")
cat("  Method:", hvg_method, "\n")
cat("  Number of features:", n_hvg, "\n")

seu <- FindVariableFeatures(
  seu,
  selection.method = hvg_method,
  nfeatures = n_hvg
)

cat("Variable features identified:", length(VariableFeatures(seu)), "\n")

# ---- Top HVGs ----
top10 <- head(VariableFeatures(seu), 10)
top20 <- head(VariableFeatures(seu), 20)

# Map Ensembl IDs to gene symbols for display
# Note: feature metadata rownames are numeric, so use named vector lookup
if ("gene_symbol" %in% colnames(seu[["RNA"]]@meta.data)) {
  gene_meta <- seu[["RNA"]]@meta.data
  symbol_lookup <- setNames(gene_meta$gene_symbol, rownames(seu))

  top10_symbols <- symbol_lookup[top10]
  top10_symbols[is.na(top10_symbols) | top10_symbols == ""] <- top10[is.na(top10_symbols) | top10_symbols == ""]

  top20_symbols <- symbol_lookup[top20]
  top20_symbols[is.na(top20_symbols) | top20_symbols == ""] <- top20[is.na(top20_symbols) | top20_symbols == ""]

  cat("\nTop 10 HVGs (symbol - Ensembl):\n")
  for (i in seq_along(top10)) {
    cat(sprintf("  %2d. %-15s %s\n", i, top10_symbols[i], top10[i]))
  }
} else {
  top10_symbols <- top10
  top20_symbols <- top20

  cat("\nTop 10 HVGs:\n")
  for (i in seq_along(top10)) {
    cat(sprintf("  %2d. %s\n", i, top10[i]))
  }
}

# ---- Save HVG plot ----
cat("\n=============================\n")
cat("Generating HVG plot...\n")

dir.create(dirname(norm_seurat_path), recursive = TRUE, showWarnings = FALSE)

plot_path <- file.path(
  dirname(norm_seurat_path),
  paste0(dataset_id, "_HVG_plot.pdf")
)

plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(plot_path, width = 12, height = 6)
print(plot1 + plot2 + plot_annotation(title = paste("HVGs -", dataset_id)))
dev.off()

cat("HVG plot saved:", plot_path, "\n")

# ---- Save normalized object ----
cat("\n=============================\n")
cat("Saving normalized object...\n")

saveRDS(seu, norm_seurat_path)
cat("Saved:", norm_seurat_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(norm_seurat_path) / 1e6))

cat("\n=============================\n")
cat("Done!", ncol(seu), "cells,", nrow(seu), "genes,",
    length(VariableFeatures(seu)), "HVGs\n")

rm(seu)
gc()
