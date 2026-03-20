#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 04_qc_filtering.R
# Description: Filter low-quality cells
# Input: Seurat objects with QC metrics
# Output: Filtered Seurat objects
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Provide config file path.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Create output directories ----
dir.create(seurat_qc_filtered_dir, recursive = TRUE, showWarnings = FALSE)

# ---- List QC Seurat objects ----
seurat_files <- list.files(
  seurat_qc_metrics_dir,
  pattern = "_seurat_qc_metrics.rds$",
  full.names = TRUE
)

# ---- Objects list ----
seus_filtered <- list()


# ---- Loop over samples ----
for (file in seurat_files) {
  
  samp <- sub("_seurat_qc_metrics.rds", "", basename(file))
  
  message("Filtering: ", samp)

  seu <- readRDS(file)

  seu <- subset(
    seu,
    subset =
      nFeature_RNA > qc_min_features &
      nFeature_RNA < qc_max_features &
      nCount_RNA > qc_min_counts &
      percent.mt < qc_max_percent_mt
  )
  
  if (remove_doublets) {
    seu <- subset(seu, subset = scDblFinder.class == "singlet")
  }
  
  seu_filtered <- seu

  saveRDS(
    seu_filtered,
    file = file.path(
      seurat_qc_filtered_dir,
      paste0(samp, "_seurat_filtered.rds")
    )
  )
  
  rm(seu, seu_filtered)
  gc()
}