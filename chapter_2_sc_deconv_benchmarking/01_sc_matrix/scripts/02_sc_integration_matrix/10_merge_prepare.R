#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 10_merge_prepare.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_integration_matrix
# Description: Load 5 annotated Seurat objects, filter intestinal metaplasia
#              samples, join sample_sheet metadata, merge all datasets,
#              select batch-aware HVGs, and run PCA as the unintegrated
#              baseline for the integration benchmark.
# Input:  5 annotated Seurat objects (output of script 09c per dataset)
#         sample_sheet.csv (metadata with patient_id, disease_group, etc.)
# Output: Merged Seurat object with PCA (.rds) + summary CSV
# Usage:  Rscript 10_merge_prepare.R <config_integration_path>
# Author: Juliana Pinto
# Date:   2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 10_merge_prepare.R <config_integration_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

cat("\n=============================\n")
cat("Script: 10_merge_prepare.R\n")
cat("Config:", config_path, "\n")
cat("Datasets:", paste(datasets, collapse = ", "), "\n")

# ---- Load sample sheet ----
cat("\n=============================\n")
cat("Loading sample sheet:", sample_sheet_path, "\n")

sample_sheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)

cat("Sample sheet rows:", nrow(sample_sheet), "\n")
cat("Datasets in sheet:", paste(unique(sample_sheet$dataset_id), collapse = ", "), "\n")
cat("Disease groups:", paste(unique(sample_sheet$disease_group), collapse = ", "), "\n")

# ---- Determine samples to exclude ----
samples_to_exclude <- sample_sheet$sample_id[
  sample_sheet$disease_group %in% disease_groups_to_exclude
]

cat("\nSamples excluded (disease_groups_to_exclude =",
    paste(disease_groups_to_exclude, collapse = ", "), "):\n")
for (s in samples_to_exclude) {
  dg <- sample_sheet$disease_group[sample_sheet$sample_id == s]
  cat(sprintf("  %s  [%s]\n", s, dg))
}
cat("Total excluded:", length(samples_to_exclude), "samples\n")

# ---- Load annotated objects ----
cat("\n=============================\n")
cat("Loading annotated Seurat objects...\n")

seu_list <- list()

for (ds in datasets) {

  obj_path <- annotated_object_paths[[ds]]

  if (!file.exists(obj_path)) {
    stop(sprintf("Annotated object not found for %s:\n  %s", ds, obj_path))
  }

  cat(sprintf("  Loading %s... ", ds))
  seu <- readRDS(obj_path)
  cat(sprintf("%d cells, %d genes\n", ncol(seu), nrow(seu)))

  # Verify dataset_id is set in metadata
  if (!"dataset_id" %in% colnames(seu@meta.data)) {
    cat(sprintf("    Adding dataset_id = '%s' to metadata\n", ds))
    seu$dataset_id <- ds
  }

  seu_list[[ds]] <- seu
  rm(seu)
  gc()
}

# ---- Join sample_sheet metadata ----
cat("\n=============================\n")
cat("Joining sample_sheet metadata to each object...\n")

# Columns to add from sample_sheet (beyond what's already in Seurat metadata)
sheet_cols <- c("patient_id", "condition", "disease_group",
                "hp_status", "treatment_status", "sex", "age",
                "histology")

for (ds in datasets) {

  seu <- seu_list[[ds]]
  cell_barcodes <- colnames(seu)

  sheet_subset <- sample_sheet[sample_sheet$dataset_id == ds, ]
  row_idx <- match(trimws(as.character(seu$orig.ident)), sheet_subset$sample_id)
  for (col in sheet_cols) {
    seu[[col]] <- sheet_subset[[col]][row_idx]
  }

  n_na <- sum(is.na(seu$disease_group))
  if (n_na > 0) {
    warning(sprintf("%s: %d cells have NA disease_group after join — check orig.ident vs sample_id", ds, n_na))
  }

  cat(sprintf("  %s: metadata joined (%d cells, disease_groups: %s)\n",
              ds,
              ncol(seu),
              paste(unique(seu$disease_group), collapse = "/")))

  seu_list[[ds]] <- seu
  rm(seu)
}

# ---- Filter intestinal metaplasia ----
cat("\n=============================\n")
cat("Filtering disease_group =", paste(disease_groups_to_exclude, collapse = "/"), "\n")

cell_counts_before  <- sapply(seu_list, ncol)
missing_samples_log <- setNames(vector("list", length(datasets)), datasets)

for (ds in datasets) {
  seu <- seu_list[[ds]]

  # Warn about samples absent from the sample_sheet (likely a naming error)
  missing_samples <- setdiff(
    unique(trimws(as.character(seu$orig.ident))),
    sample_sheet$sample_id[sample_sheet$dataset_id == ds]
  )
  missing_samples_log[[ds]] <- missing_samples
  if (length(missing_samples) > 0) {
    warning(sprintf(
      "%s: %d sample(s) not found in sample_sheet — will be excluded: %s\n  Check for typos or missing entries.",
      ds, length(missing_samples), paste(missing_samples, collapse = ", ")
    ))
  }

  # Remove IM samples AND cells from samples absent from the sample_sheet (NA disease_group)
  keep <- !is.na(seu$disease_group) & !seu$disease_group %in% disease_groups_to_exclude
  if (!all(keep)) {
    seu_list[[ds]] <- subset(seu, cells = colnames(seu)[keep])
  }
  rm(seu)
}

cell_counts_after <- sapply(seu_list, ncol)

cat("\nCell counts (before / after IM filter):\n")
total_before <- 0L
total_after  <- 0L

for (ds in datasets) {
  before <- cell_counts_before[[ds]]
  after  <- cell_counts_after[[ds]]
  removed <- before - after
  total_before <- total_before + before
  total_after  <- total_after  + after
  cat(sprintf("  %-12s  %6d  →  %6d  (-%d cells)\n", ds, before, after, removed))
}
cat(sprintf("  %-12s  %6d  →  %6d  (-%d cells)\n",
            "TOTAL", total_before, total_after, total_before - total_after))

# ---- Batch-aware HVG selection ----
# FindVariableFeatures is run per dataset to avoid batch-specific genes
# dominating the HVG list. SelectIntegrationFeatures then picks features
# that are variable in multiple batches.
cat("\n=============================\n")
cat("Batch-aware HVG selection (SelectIntegrationFeatures)...\n")
cat("  nfeatures:", n_hvg_integration, "\n")

for (ds in datasets) {
  cat(sprintf("  FindVariableFeatures: %s... ", ds))
  seu_list[[ds]] <- FindVariableFeatures(
    seu_list[[ds]],
    selection.method = "vst",
    nfeatures        = n_hvg_integration,
    verbose          = FALSE
  )
  cat(sprintf("%d HVGs\n", length(VariableFeatures(seu_list[[ds]]))))
}

hvg_features <- SelectIntegrationFeatures(
  object.list = seu_list,
  nfeatures   = n_hvg_integration,
  verbose     = FALSE
)

cat(sprintf("Batch-aware HVGs selected: %d\n", length(hvg_features)))

# ---- Merge all datasets ----
cat("\n=============================\n")
cat("Merging", length(seu_list), "objects...\n")

merged <- Reduce(function(x, y) merge(x, y), seu_list)

cat("Merged dimensions:", nrow(merged), "genes x", ncol(merged), "cells\n")
cat("Datasets present:", paste(unique(merged$dataset_id), collapse = ", "), "\n")
cat("Disease groups:", paste(unique(merged$disease_group), collapse = ", "), "\n")

rm(seu_list)
gc()

# ---- Join layers ----
cat("Joining layers...\n")
merged <- JoinLayers(merged)

# ---- Apply batch-aware HVGs to merged object ----
VariableFeatures(merged) <- hvg_features
cat("Variable features set on merged object:", length(VariableFeatures(merged)), "\n")

# ---- Normalize ----
# RC normalization is consistent with the per-dataset pipeline (deconvolution reference)
cat("\n=============================\n")
cat("Normalizing (RC, scale factor 10000)...\n")
merged <- NormalizeData(
  merged,
  normalization.method = "RC",
  scale.factor         = 10000,
  verbose              = FALSE
)

# ---- Scale data ----
# Scale only on HVGs for PCA (ScaleData on all genes would be memory-prohibitive)
cat("Scaling data (HVGs only)...\n")
merged <- ScaleData(merged, features = hvg_features, verbose = FALSE)

# ---- PCA (unintegrated baseline) ----
cat("\n=============================\n")
cat("Running PCA (unintegrated baseline)...\n")
cat("  npcs:", n_pcs_integration, "\n")

merged <- RunPCA(
  merged,
  features = hvg_features,
  npcs     = n_pcs_integration,
  verbose  = FALSE
)

stdev       <- merged[["pca"]]@stdev
var_exp     <- (stdev^2 / sum(stdev^2)) * 100
dims_to_use <- dims_integration

cat(sprintf("  Variance explained by PC1:        %.1f%%\n", var_exp[1]))
cat(sprintf("  Variance explained by PC1-PC%d:  %.1f%%\n",
            max(dims_to_use), sum(var_exp[dims_to_use])))

# ---- Save elbow plot ----
dir.create(dirname(merged_integration_path), recursive = TRUE, showWarnings = FALSE)

elbow_path <- file.path(
  dirname(merged_integration_path),
  "unintegrated_elbow_plot.pdf"
)

pdf(elbow_path, width = 8, height = 5)
print(
  ElbowPlot(merged, ndims = n_pcs_integration, reduction = "pca") +
    geom_vline(xintercept = max(dims_to_use), linetype = "dashed", color = "red") +
    ggtitle(paste("ElbowPlot — unintegrated baseline — cutoff PC", max(dims_to_use)))
)
dev.off()
cat("ElbowPlot saved:", elbow_path, "\n")

# ---- UMAP (unintegrated — visual baseline for benchmark) ----
cat("\n=============================\n")
cat("Running UMAP (unintegrated)...\n")

merged <- RunUMAP(
  merged,
  dims      = dims_to_use,
  reduction = "pca",
  verbose   = FALSE
)

umap_path <- file.path(dirname(merged_integration_path), "umap_unintegrated.pdf")

dataset_colors <- c(
  GSE163558 = "#E41A1C",
  GSE275648 = "#377EB8",
  GSE308231 = "#4DAF4A",
  GSE201347 = "#FF7F00",
  GSE249874 = "#984EA3"
)

use_raster <- requireNamespace("ggrastr", quietly = TRUE)
cat(sprintf("Rasterized UMAP: %s\n", if (use_raster) "yes (ggrastr)" else "no"))

pdf(umap_path, width = 16, height = 7)
p1 <- DimPlot(merged, group.by = "dataset_id",    reduction = "umap",
              pt.size = 0.5, raster = use_raster,
              raster.dpi = if (use_raster) c(1024, 1024) else NULL,
              cols = dataset_colors) +
        ggtitle("Unintegrated — by dataset") +
        guides(colour = guide_legend(override.aes = list(size = 3)))
p2 <- DimPlot(merged, group.by = "disease_group", reduction = "umap",
              pt.size = 0.5, raster = use_raster,
              raster.dpi = if (use_raster) c(1024, 1024) else NULL) +
        ggtitle("Unintegrated — by disease_group") +
        guides(colour = guide_legend(override.aes = list(size = 3)))
print(p1 | p2)
dev.off()
cat("UMAP plot saved:", umap_path, "\n")

# ---- Summary CSV ----
summary_df <- data.frame(
  dataset_id              = datasets,
  cells_before            = as.integer(cell_counts_before),
  cells_after             = as.integer(cell_counts_after),
  cells_removed           = as.integer(cell_counts_before - cell_counts_after),
  samples_not_in_sheet    = sapply(missing_samples_log, function(x) {
                              if (length(x) == 0) "" else paste(x, collapse = "; ")
                            }),
  stringsAsFactors = FALSE
)

summary_path <- file.path(dirname(merged_integration_path), "merge_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
cat("Summary CSV saved:", summary_path, "\n")

# ---- Save merged object ----
cat("\n=============================\n")
cat("Saving merged object...\n")

saveRDS(merged, merged_integration_path)
cat("Saved:", merged_integration_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(merged_integration_path) / 1e6))

cat("\n=============================\n")
cat("Done!\n")
cat("  Total cells:", ncol(merged), "\n")
cat("  Genes:", nrow(merged), "\n")
cat("  HVGs:", length(VariableFeatures(merged)), "\n")
cat("  Datasets:", length(unique(merged$dataset_id)), "\n")
cat("Next: run 11_integration_benchmark.py to compare integration methods\n")

rm(merged)
gc()
