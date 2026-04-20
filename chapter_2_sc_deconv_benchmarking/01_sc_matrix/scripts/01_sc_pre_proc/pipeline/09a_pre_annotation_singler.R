#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 09a_pre_annotation_singler.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Automatic pre-annotation with SingleR using two references.
#              This is a preliminary annotation to understand cell composition
#              per dataset BEFORE integration. Final annotation will be done
#              on the integrated object.
# Input: Clustered Seurat object (from script 07b)
# Output: Seurat object with SingleR columns + plots (PDF)
# Usage: Rscript 09a_pre_annotation_singler.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(reticulate)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 09a_pre_annotation_singler.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load clustered object ----
cat("\n=============================\n")
cat("Script: 09a_pre_annotation_singler.R\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading:", basename(clustered_seurat_path), "\n")

seu <- readRDS(clustered_seurat_path)

cat("Dimensions:", nrow(seu), "genes x", ncol(seu), "cells\n")
cat("Samples:", length(unique(seu$orig.ident)), "\n")

# ---- Set cluster identity ----
annotation_cluster_col <- paste0("RNA_snn_res.", annotation_resolution)

if (!annotation_cluster_col %in% colnames(seu@meta.data)) {
  stop("Column '", annotation_cluster_col, "' not found.")
}

Idents(seu) <- annotation_cluster_col
n_clusters <- length(unique(Idents(seu)))
cat("Resolution:", annotation_resolution, "(", n_clusters, "clusters )\n")

# ---- Swap rownames to gene symbols (SingleR references use symbols, not Ensembl) ----
cat("\n=============================\n")
cat("Preparing gene symbols for SingleR...\n")

gene_meta <- seu[["RNA"]]@meta.data
if ("gene_symbol" %in% colnames(gene_meta)) {
  symbol_lookup <- setNames(gene_meta$gene_symbol, rownames(seu))
  # Keep Ensembl ID if symbol is NA or empty
  new_names <- ifelse(is.na(symbol_lookup) | symbol_lookup == "",
                      names(symbol_lookup), symbol_lookup)
  # Handle duplicated symbols by appending Ensembl ID
  dup_mask <- duplicated(new_names) | duplicated(new_names, fromLast = TRUE)
  new_names[dup_mask] <- paste0(new_names[dup_mask], "_", names(symbol_lookup)[dup_mask])
  cat("  Genes with symbols:", sum(!is.na(symbol_lookup) & symbol_lookup != ""), "/", length(symbol_lookup), "\n")
  cat("  Duplicated symbols resolved:", sum(dup_mask), "\n")
} else {
  stop("gene_symbol column not found in feature metadata. Cannot run SingleR.")
}

# ---- Convert to SCE with gene symbols ----
cat("Converting to SingleCellExperiment...\n")

sce <- as.SingleCellExperiment(seu)
rownames(sce) <- new_names

cat("SCE ready -", nrow(sce), "genes\n")

# ---- Load references ----
cat("\n=============================\n")
cat("Loading SingleR references...\n")

ref_hpca <- celldex::HumanPrimaryCellAtlasData()
cat("  HumanPrimaryCellAtlasData:", ncol(ref_hpca), "samples,", length(unique(ref_hpca$label.main)), "main labels\n")

ref_blueprint <- celldex::BlueprintEncodeData()
cat("  BlueprintEncodeData:", ncol(ref_blueprint), "samples,", length(unique(ref_blueprint$label.main)), "main labels\n")

# ---- Run SingleR - HPCA ----
cat("\n=============================\n")
cat("Running SingleR with HumanPrimaryCellAtlasData...\n")

pred_hpca <- SingleR(
  test = sce,
  ref = ref_hpca,
  labels = ref_hpca$label.main
)

seu$singler_hpca <- pred_hpca$labels
cat("HPCA annotation complete\n")
cat("Labels found:", length(unique(pred_hpca$labels)), "\n")
cat("\nHPCA label counts:\n")
print(sort(table(pred_hpca$labels), decreasing = TRUE))

# ---- Run SingleR - Blueprint ----
cat("\n=============================\n")
cat("Running SingleR with BlueprintEncodeData...\n")

pred_blueprint <- SingleR(
  test = sce,
  ref = ref_blueprint,
  labels = ref_blueprint$label.main
)

seu$singler_blueprint <- pred_blueprint$labels
cat("Blueprint annotation complete\n")
cat("Labels found:", length(unique(pred_blueprint$labels)), "\n")
cat("\nBlueprint label counts:\n")
print(sort(table(pred_blueprint$labels), decreasing = TRUE))

# ---- Output directory ----
dir.create(annotation_output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Save plots ----
cat("\n=============================\n")
cat("Generating plots...\n")

# UMAP: clusters vs HPCA vs Blueprint
pdf(file.path(annotation_output_dir, paste0(dataset_id, "_singler_umap.pdf")), width = 18, height = 6)
p1 <- DimPlot(seu, group.by = annotation_cluster_col, label = TRUE) +
  ggtitle(paste("Clusters (res", annotation_resolution, ")")) + NoLegend()
p2 <- DimPlot(seu, group.by = "singler_hpca", label = TRUE, repel = TRUE) +
  ggtitle("SingleR - HPCA") + NoLegend()
p3 <- DimPlot(seu, group.by = "singler_blueprint", label = TRUE, repel = TRUE) +
  ggtitle("SingleR - Blueprint") + NoLegend()
print(p1 + p2 + p3)
dev.off()
cat("UMAP plots saved\n")

# Confusion matrix: cluster vs HPCA
ct_hpca <- table(
  Cluster = seu@meta.data[[annotation_cluster_col]],
  HPCA = seu$singler_hpca
)
ct_hpca_prop <- prop.table(ct_hpca, margin = 1) * 100

pdf(file.path(annotation_output_dir, paste0(dataset_id, "_singler_hpca_confusion.pdf")), width = 14, height = 8)
pheatmap(
  ct_hpca_prop,
  fontsize = 7,
  display_numbers = TRUE,
  number_format = "%.0f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "steelblue"))(100),
  main = paste("Cluster vs SingleR-HPCA (%) -", dataset_id)
)
dev.off()
cat("HPCA confusion matrix saved\n")

# Confusion matrix: cluster vs Blueprint
ct_bp <- table(
  Cluster = seu@meta.data[[annotation_cluster_col]],
  Blueprint = seu$singler_blueprint
)
ct_bp_prop <- prop.table(ct_bp, margin = 1) * 100

pdf(file.path(annotation_output_dir, paste0(dataset_id, "_singler_blueprint_confusion.pdf")), width = 14, height = 8)
pheatmap(
  ct_bp_prop,
  fontsize = 7,
  display_numbers = TRUE,
  number_format = "%.0f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "steelblue"))(100),
  main = paste("Cluster vs SingleR-Blueprint (%) -", dataset_id)
)
dev.off()
cat("Blueprint confusion matrix saved\n")

# ---- Save annotated object ----
cat("\n=============================\n")
cat("Saving annotated object...\n")

saveRDS(seu, annotated_seurat_path)
cat("Saved:", annotated_seurat_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(annotated_seurat_path) / 1e6))

# ---- Save SingleR labels as CSV (for reference) ----
labels_df <- data.frame(
  cell_id = colnames(seu),
  cluster = seu@meta.data[[annotation_cluster_col]],
  singler_hpca = seu$singler_hpca,
  singler_blueprint = seu$singler_blueprint
)
labels_csv <- file.path(annotation_output_dir, paste0(dataset_id, "_singler_labels.csv"))
write.csv(labels_df, labels_csv, row.names = FALSE)
cat("Labels CSV saved:", labels_csv, "\n")

# ---- Convert to h5ad (for CellTypist) ----
# Export manually via reticulate to avoid sceasy/Seurat v5 slot compatibility issues
cat("\n=============================\n")
cat("Converting to h5ad for CellTypist...\n")

# Get normalized counts with gene symbols as rownames
counts_for_export <- GetAssayData(seu, layer = "data")

# Swap rownames to gene symbols (CellTypist needs symbols)
gene_meta <- seu[["RNA"]]@meta.data
if ("gene_symbol" %in% colnames(gene_meta)) {
  symbol_map <- setNames(gene_meta$gene_symbol, rownames(seu))
  export_names <- ifelse(is.na(symbol_map) | symbol_map == "",
                         names(symbol_map), symbol_map)
  # Deduplicate
  dup <- duplicated(export_names)
  export_names[dup] <- paste0(export_names[dup], "_", names(symbol_map)[dup])
  rownames(counts_for_export) <- export_names
}

# Export via reticulate + anndata
anndata <- reticulate::import("anndata")
scipy_sparse <- reticulate::import("scipy.sparse")
np <- reticulate::import("numpy")

# Transpose: R sparse matrix is genes x cells, anndata needs cells x genes
mat_t <- Matrix::t(counts_for_export)  # now cells x genes, still dgCMatrix (CSC format)

# dgCMatrix (CSC) -> scipy csc_matrix (same format, direct slot mapping)
sp_mat <- scipy_sparse$csc_matrix(
  reticulate::tuple(
    np$array(mat_t@x),
    np$array(as.integer(mat_t@i)),
    np$array(as.integer(mat_t@p))
  ),
  shape = reticulate::tuple(as.integer(nrow(mat_t)), as.integer(ncol(mat_t)))
)

# Create AnnData
obs_df <- seu@meta.data
obs_df$cell_id <- rownames(obs_df)  # preserve cell barcodes

adata <- anndata$AnnData(
  X = sp_mat,
  obs = obs_df,
  var = data.frame(gene = colnames(mat_t), row.names = colnames(mat_t))
)

# Add UMAP if available
if ("umap" %in% names(seu@reductions)) {
  umap_coords <- Embeddings(seu, "umap")
  adata$obsm$update(list("X_umap" = np$array(umap_coords)))
}

# Save
adata$write_h5ad(h5ad_path)

cat("h5ad saved:", h5ad_path, "\n")
cat(sprintf("File size: %.1f MB\n", file.size(h5ad_path) / 1e6))

# ---- Done ----
cat("\n=============================\n")
cat("Done! SingleR pre-annotation complete\n")
cat("  HPCA labels:", length(unique(seu$singler_hpca)), "cell types\n")
cat("  Blueprint labels:", length(unique(seu$singler_blueprint)), "cell types\n")
cat("Next: run CellTypist (09b), then import (09c), then explore in notebook\n")

rm(seu, sce, pred_hpca, pred_blueprint)
gc()
