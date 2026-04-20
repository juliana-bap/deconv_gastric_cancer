#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 07a_dimred.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: ScaleData + PCA + ElbowPlot.
#              Run once per dataset. Use ElbowPlot to choose dims_to_use
#              before running 07b_clustering.R.
# Input: Normalized Seurat object with HVGs (from script 06)
# Output: Seurat object with PCA (.rds) + ElbowPlot (PDF)
# Usage: Rscript 07a_dimred.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 07a_dimred.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load normalized object ----
cat("\n=============================\n")
cat("Script: 07a_dimred.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading:", basename(norm_seurat_path), "\n")

seu <- readRDS(norm_seurat_path)

cat("Dimensions:", nrow(seu), "genes x", ncol(seu), "cells\n")
cat("HVGs:", length(VariableFeatures(seu)), "\n")

# ---- Scale data ----
# ScaleData is required for PCA (centers and scales genes).
# This is standard within-dataset preprocessing, NOT the same as
# scaling before integration (which Luecken et al. 2021 cautions about).
# The integration step will decide on scaling independently.
cat("\n=============================\n")
cat("Scaling data (required for PCA)...\n")

seu <- ScaleData(seu, verbose = FALSE)

cat("Scaling complete\n")

# ---- PCA ----
cat("\n=============================\n")
cat("Running PCA...\n")
cat("  n_pcs:", n_pcs, "\n")

seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)

cat("PCA complete\n")

# Print variance explained by top PCs
stdev <- seu[["pca"]]@stdev
var_explained <- (stdev^2 / sum(stdev^2)) * 100
cat("\nVariance explained:\n")
cat(sprintf("  PC1:       %.1f%%\n", var_explained[1]))
cat(sprintf("  PC1-PC%d: %.1f%%\n", max(dims_to_use), sum(var_explained[dims_to_use])))
cat(sprintf("  PC1-PC%d: %.1f%%\n", n_pcs, sum(var_explained)))

# ---- Save ElbowPlot ----
dir.create(dirname(pca_seurat_path), recursive = TRUE, showWarnings = FALSE)

elbow_path <- file.path(dirname(pca_seurat_path), paste0(dataset_id, "_elbow_plot.pdf"))

pdf(elbow_path, width = 8, height = 5)
print(
  ElbowPlot(seu, ndims = n_pcs, reduction = "pca") +
    geom_vline(xintercept = max(dims_to_use), linetype = "dashed", color = "red") +
    ggtitle(paste("ElbowPlot -", dataset_id, "- current cutoff at PC", max(dims_to_use)))
)
dev.off()

cat("ElbowPlot saved:", elbow_path, "\n")

# ---- Save object ----
cat("\n=============================\n")
cat("Saving PCA object...\n")

saveRDS(seu, pca_seurat_path)
cat("Saved:", pca_seurat_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(pca_seurat_path) / 1e6))

cat("\n=============================\n")
cat("Done! PCA computed with", n_pcs, "PCs\n")
cat("Next: check ElbowPlot, adjust dims_to_use in config, then run 07b_clustering.R\n")

rm(seu)
gc()
