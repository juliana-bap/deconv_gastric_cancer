#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 07b_clustering.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: FindNeighbors, FindClusters (Leiden), UMAP.
#              Can be re-run with different dims_to_use or resolutions
#              without repeating PCA (script 07a).
# Input: Seurat object with PCA (from script 07a)
# Output: Clustered Seurat object (.rds) + UMAP plots (PDF)
# Usage: Rscript 07b_clustering.R <config_path>
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
  stop("Usage: Rscript 07b_clustering.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load PCA object ----
cat("\n=============================\n")
cat("Script: 07b_clustering.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading:", basename(pca_seurat_path), "\n")

seu <- readRDS(pca_seurat_path)

cat("Dimensions:", nrow(seu), "genes x", ncol(seu), "cells\n")

# ---- Find Neighbors ----
cat("\n=============================\n")
cat("Finding neighbors...\n")
cat("  dims:", min(dims_to_use), "to", max(dims_to_use), "\n")

seu <- FindNeighbors(seu, dims = dims_to_use)

cat("Neighbors computed\n")

# ---- Clustering ----
cat("\n=============================\n")
cat("Clustering...\n")
cat("  Algorithm:", cluster_algorithm, "(", ifelse(cluster_algorithm == 4, "Leiden", "Louvain"), ")\n")
cat("  Resolutions:", paste(cluster_resolutions, collapse = ", "), "\n")

seu <- FindClusters(
  seu,
  algorithm = cluster_algorithm,
  resolution = cluster_resolutions
)

# Print cluster counts per resolution
cat("\nClusters per resolution:\n")
for (res in cluster_resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  if (col_name %in% colnames(seu@meta.data)) {
    n_clusters <- length(unique(seu@meta.data[[col_name]]))
    cat(sprintf("  res %.2f: %d clusters\n", res, n_clusters))
  }
}

# ---- UMAP ----
cat("\n=============================\n")
cat("Running UMAP...\n")
cat("  dims:", min(dims_to_use), "to", max(dims_to_use), "\n")

seu <- RunUMAP(
  seu,
  dims = dims_to_use,
  reduction = "pca",
  verbose = FALSE
)

cat("UMAP complete\n")

# ---- Save UMAP plots ----
dir.create(dirname(clustered_seurat_path), recursive = TRUE, showWarnings = FALSE)

umap_path <- file.path(dirname(clustered_seurat_path), paste0(dataset_id, "_umap_clusters.pdf"))

plots <- list()
for (res in cluster_resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  if (col_name %in% colnames(seu@meta.data)) {
    plots[[as.character(res)]] <- DimPlot(seu, group.by = col_name, label = TRUE) +
      ggtitle(paste0("res = ", res)) +
      NoLegend()
  }
}

# UMAP colored by sample
plots[["sample"]] <- DimPlot(seu, group.by = "orig.ident") +
  ggtitle("By sample")

pdf(umap_path, width = 14, height = 10)
print(wrap_plots(plots, ncol = 2))
dev.off()

cat("UMAP plots saved:", umap_path, "\n")

# ---- Save object ----
cat("\n=============================\n")
cat("Saving clustered object...\n")

saveRDS(seu, clustered_seurat_path)
cat("Saved:", clustered_seurat_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(clustered_seurat_path) / 1e6))

cat("\n=============================\n")
cat("Done!", ncol(seu), "cells,", nrow(seu), "genes\n")
cat("Clusters:", paste(sapply(cluster_resolutions, function(r) {
  col <- paste0("RNA_snn_res.", r)
  paste0("res", r, "=", length(unique(seu@meta.data[[col]])))
}), collapse = ", "), "\n")

rm(seu)
gc()
