#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 08_rogue.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Compute ROGUE scores per cluster x sample.
#              ROGUE measures cluster purity — high scores (close to 1)
#              indicate homogeneous clusters, low scores suggest mixed
#              cell types or poor clustering.
# Input: Clustered Seurat object (from script 07b)
# Output: ROGUE scores table (CSV) + boxplot + heatmap (PDF)
# Usage: Rscript 08_rogue.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ROGUE)
  library(tibble)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 08_rogue.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load clustered object ----
cat("\n=============================\n")
cat("Script: 08_rogue.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading:", basename(clustered_seurat_path), "\n")

seu <- readRDS(clustered_seurat_path)

cat("Dimensions:", nrow(seu), "genes x", ncol(seu), "cells\n")
cat("Samples:", length(unique(seu$orig.ident)), "\n")

# ---- Select resolution ----
rogue_cluster_col <- paste0("RNA_snn_res.", rogue_resolution)

if (!rogue_cluster_col %in% colnames(seu@meta.data)) {
  stop("Column '", rogue_cluster_col, "' not found. Available: ",
       paste(grep("RNA_snn_res", colnames(seu@meta.data), value = TRUE), collapse = ", "))
}

n_clusters <- length(unique(seu@meta.data[[rogue_cluster_col]]))
cat("Resolution:", rogue_resolution, "(", n_clusters, "clusters )\n")

# ---- Downsampling (optional) ----
Idents(seu) <- rogue_cluster_col

if (exists("rogue_downsample") && !is.null(rogue_downsample)) {
  cat("\n=============================\n")
  cat("Downsampling to", rogue_downsample, "cells per cluster...\n")
  seu_ds <- subset(x = seu, downsample = rogue_downsample)
  cat("Cells after downsampling:", ncol(seu_ds), "(original:", ncol(seu), ")\n")
} else {
  cat("\n=============================\n")
  cat("No downsampling - using all", ncol(seu), "cells\n")
  seu_ds <- seu
}

# ---- Extract count matrix ----
cat("\n=============================\n")
cat("Extracting count matrix...\n")

counts_mat <- LayerData(seu_ds, assay = "RNA", layer = "counts")

# Remove genes with zero expression (reduces matrix size before dense conversion)
gene_sums <- Matrix::rowSums(counts_mat)
counts_mat <- counts_mat[gene_sums > 0, ]
cat("Matrix after removing zero-expression genes:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

# Convert to dense matrix (ROGUE requires it)
cat("Converting to dense matrix...\n")
counts_mat <- as.matrix(counts_mat)
cat(sprintf("Memory: %.1f MB\n", object.size(counts_mat) / 1e6))

# ---- Filter genes with SE ----
cat("\n=============================\n")
cat("Filtering genes with ROGUE SE method...\n")

expr_entropy <- SE_fun(counts_mat)

cat("Genes after SE filtering:", nrow(expr_entropy), "\n")

# ---- Global ROGUE score ----
cat("\n=============================\n")
cat("Computing global ROGUE score...\n")

global_rogue <- CalculateRogue(expr_entropy, platform = "UMI")
cat("Global ROGUE:", round(global_rogue, 4), "\n")

# ---- Compute ROGUE per cluster x sample ----
cat("\n=============================\n")
cat("Computing ROGUE scores per cluster x sample...\n")
cat("  Cluster column:", rogue_cluster_col, "\n")
cat("  Sample column: orig.ident\n")

rogue_scores <- rogue(
  counts_mat,
  labels = seu_ds@meta.data[[rogue_cluster_col]],
  samples = seu_ds$orig.ident,
  platform = "UMI",
  span = 0.6
)

cat("ROGUE computation complete\n")

# ---- Output directory ----
dir.create(rogue_output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Save ROGUE scores ----
rogue_csv_path <- file.path(rogue_output_dir, paste0(dataset_id, "_rogue_scores_res", rogue_resolution, ".csv"))
write.csv(rogue_scores, rogue_csv_path, row.names = TRUE)
cat("ROGUE scores saved:", rogue_csv_path, "\n")

# ---- Print summary per cluster ----
cat("\n=============================\n")
cat("ROGUE score summary per cluster:\n")

rogue_mat <- as.matrix(rogue_scores)
for (cl in colnames(rogue_mat)) {
  vals <- rogue_mat[, cl]
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    cat(sprintf("  Cluster %s: mean = %.3f (range: %.3f - %.3f)\n",
                cl, mean(vals), min(vals), max(vals)))
  }
}

overall_mean <- mean(rogue_mat, na.rm = TRUE)
cat(sprintf("\nOverall mean ROGUE: %.3f\n", overall_mean))

# ---- Boxplot: ROGUE built-in ----
box_path <- file.path(rogue_output_dir, paste0(dataset_id, "_rogue_boxplot_res", rogue_resolution, ".pdf"))
pdf(box_path, width = 10, height = 6)
print(rogue.boxplot(rogue_scores))
dev.off()
cat("Boxplot saved:", box_path, "\n")

# ---- Heatmap: ROGUE cluster x sample ----
rogue_long <- as.data.frame(as.table(rogue_mat))
colnames(rogue_long) <- c("sample", "cluster", "rogue")
rogue_long$rogue <- as.numeric(rogue_long$rogue)
rogue_long <- rogue_long[!is.na(rogue_long$rogue), ]

p_heat <- ggplot(rogue_long, aes(x = cluster, y = sample, fill = rogue)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rogue)), size = 3) +
  scale_fill_gradient2(low = "red", mid = "white", high = "steelblue",
                       midpoint = 0.5, limits = c(0, 1)) +
  theme_minimal() +
  ggtitle(paste("ROGUE heatmap -", dataset_id, "(res", rogue_resolution, ")")) +
  xlab("Cluster") + ylab("Sample")

heat_path <- file.path(rogue_output_dir, paste0(dataset_id, "_rogue_heatmap_res", rogue_resolution, ".pdf"))
pdf(heat_path, width = 12, height = 5)
print(p_heat)
dev.off()
cat("Heatmap saved:", heat_path, "\n")

# ---- Done ----
cat("\n=============================\n")
cat("Done! ROGUE analysis complete\n")
cat("  Resolution:", rogue_resolution, "\n")
cat("  Clusters:", n_clusters, "\n")
cat("  Global ROGUE:", round(global_rogue, 4), "\n")
cat("  Overall mean ROGUE:", round(overall_mean, 3), "\n")
cat("  Interpretation:\n")
cat("    > 0.8 = pure, homogeneous cluster\n")
cat("    0.5-0.8 = moderate purity\n")
cat("    < 0.5 = potentially mixed cluster\n")

rm(seu, seu_ds, counts_mat, expr_entropy)
gc()
