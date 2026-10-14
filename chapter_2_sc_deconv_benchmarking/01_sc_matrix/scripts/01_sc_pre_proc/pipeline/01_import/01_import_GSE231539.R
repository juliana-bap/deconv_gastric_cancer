#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 01_import_GSE231539.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE231539
# Description: Import pre-built Seurat object (GC.NAC.Anti.PD1.rds, ~5.6 GB),
#              extract count matrices, convert gene IDs to Ensembl if needed,
#              split by sample, save individual sparse RDS files.
#              File contains GC samples (NAC and Anti-PD1 treatment groups).
#              ESCC and normal esophageal samples are in a separate file and
#              are out of scope — exclude them if present.
# Input:  data/sc_reference/raw/GSE231539/GSE231539_GC.NAC.Anti.PD1.rds
# Output: data/sc_reference/raw/GSE231539/GSE231539_FIXED/<sample>.rds
# Usage:  Rscript 01_import_GSE231539.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# Ensembl version: 109 (GRCh38)
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript 01_import_GSE231539.R <config_path>\n  No config file provided.")
}
if (!file.exists(args[1])) {
  stop("Config file not found: ", args[1])
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Load functions ----
source(utils_path)

# ---- Libraries ----
# NOTE: Do NOT load Seurat before readRDS.
# The RDS may have been saved with a different Seurat version,
# causing S4 class validation errors on load.
# Loading only Matrix allows R to read the object as a generic S4,
# and we extract counts/metadata directly via @ slots.
suppressPackageStartupMessages({
  library(Matrix)
})

# ---- Check input file ----
if (!file.exists(seurat_rds_path)) {
  # Try .rds.gz if .rds not found (gunzip may not have run yet)
  gz_path <- paste0(seurat_rds_path, ".gz")
  if (file.exists(gz_path)) {
    stop("Found .rds.gz but not .rds — gunzip the file first:\n  gunzip ", gz_path)
  }
  stop("RDS file not found: ", seurat_rds_path)
}

# ---- Load Seurat object (without Seurat library) ----
cat("\n=============================\n")
cat("Dataset:", dataset_id, "\n")
cat("Loading Seurat object:", basename(seurat_rds_path), "\n")
cat("File size:", round(file.size(seurat_rds_path) / 1e9, 1), "GB\n")

seu <- readRDS(seurat_rds_path)

cat("Object loaded successfully\n")

# ---- Extract count matrix (handle Seurat v3/v4 and v5 structures) ----
counts <- tryCatch(
  seu@assays$RNA@counts,
  error = function(e1) {
    tryCatch(
      seu@assays$RNA@layers$counts,
      error = function(e2) {
        stop("Could not extract count matrix. Tried @counts (v3/v4) and @layers$counts (v5).")
      }
    )
  }
)

cat("Count matrix dimensions:", nrow(counts), "x", ncol(counts), "\n")

# ---- Extract cell metadata ----
cell_metadata <- seu@meta.data

cat("Available metadata columns:", paste(colnames(cell_metadata), collapse = ", "), "\n")

# ---- Identify samples ----
if (!sample_column %in% colnames(cell_metadata)) {
  cat("\nERROR: sample column '", sample_column, "' not found.\n")
  cat("Available columns:", paste(colnames(cell_metadata), collapse = ", "), "\n")
  cat("Update 'sample_column' in the config and re-run.\n")
  stop("Sample column not found.")
}

samples <- unique(cell_metadata[[sample_column]])
cat("Number of samples:", length(samples), "\n")
cat("Samples:", paste(sort(samples), collapse = ", "), "\n")

# ---- Check for non-GC samples ----
# The GC.NAC.Anti.PD1 file should contain only GC samples.
# Expected IDs: G48T, G54T, G57T, G60T, G65T.
# If ESCC or normal esophageal samples are present, exclude them.
escc_pattern  <- "^E[0-9]+(SCC|N)$"
normal_esoph  <- grepl(escc_pattern, samples, perl = TRUE)

if (any(normal_esoph)) {
  cat("\nWARNING: Non-GC samples detected and will be excluded:\n")
  cat("  ", paste(samples[normal_esoph], collapse = ", "), "\n")
}

# ---- Free Seurat object early (keep only counts + metadata) ----
rm(seu)
gc()

# ---- Check gene ID type and convert if needed ----
if (is_ensembl(rownames(counts))) {
  cat("-> Ensembl IDs detected\n")
  counts <- clean_ensembl_ids(counts)
} else {
  cat("-> Gene symbols detected, converting to Ensembl\n")
  mapping_table <- readRDS(mapping_table_path)
  counts <- convert_to_ensembl(counts, mapping_table)
  rm(mapping_table)
}

cat("Count matrix after conversion:", dim(counts), "\n")

# ---- Output directory ----
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Determine samples to process ----
# Exclude: non-GC (ESCC/normal esophageal) detected above
# Exclude: samples in samples_to_exclude_bio (from config — e.g. ESCC if named)
all_exclude <- union(
  samples[normal_esoph],
  samples_to_exclude_bio
)

samples_to_process <- setdiff(samples, all_exclude)

if (length(all_exclude) > 0) {
  cat("\nExcluded samples (bio):", paste(all_exclude, collapse = ", "), "\n")
}
cat("Samples to process:", paste(sort(samples_to_process), collapse = ", "), "\n")

# ---- Split by sample and save ----
for (s in samples_to_process) {

  cat("\n=============================\n")
  cat("Sample:", s, "\n")

  cells <- rownames(cell_metadata)[cell_metadata[[sample_column]] == s]
  cells <- intersect(cells, colnames(counts))
  mat <- counts[, cells]

  cat("Dimension:", dim(mat), "\n")

  outfile <- file.path(input_dir, paste0(s, ".rds"))
  saveRDS(mat, outfile)
  cat("Saved:", outfile, "\n")

  rm(mat)
  gc()
}

rm(counts, cell_metadata)
gc()

cat("\n=============================\n")
cat("Done! All GC samples saved to:", input_dir, "\n")
cat("Next step: gunzip done, run script 02 with the config.\n")
