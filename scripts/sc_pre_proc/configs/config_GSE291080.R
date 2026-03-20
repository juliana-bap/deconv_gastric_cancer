# ==============================
# CONFIG - GSE291080
# Project: Gastric Cancer Deconvolution
# Source: GEO
# Technology: 10x Genomics
# Reference genome: GRCh38
# Format: 10x (barcodes/features/matrix per sample dir)
# Gene IDs: Ensembl (from features.tsv column 1)
# Samples: 44 (1 excluded: GSM8828843_a187_es - barcode/matrix mismatch)
# ==============================

# ---- Dataset identity ----
dataset_id <- "GSE291080"

# ---- Project root ----
project_root <- "/Users/julianapinto/doutorado/deconv_gastric_cancer"

# ---- Base directories ----
raw_base_dir <- file.path(project_root, "data", "sc_reference", "raw", dataset_id)

processed_base_dir <- file.path(project_root, "data", "sc_reference", "processed", dataset_id)

# ---- Input directory (fixed RDS from script 01) ----
input_dir <- file.path(raw_base_dir, paste0(dataset_id, "_FIXED"))

# ---- Output directory (Seurat objects) ----
seurat_raw_dir <- file.path(processed_base_dir, "seurat_raw")

# ---- File list ----
rds_files <- list.files(
  input_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

# ---- Annotation ----
annotation_dir <- file.path(project_root, "annotation")
mapping_table_path <- file.path(annotation_dir, "ensembl109_full_mapping.rds")

# ---- Utils / Functions ----
utils_path <- file.path(project_root, "scripts", "sc_pre_proc", "pipeline", "00_utils_general.R")

# ---- Samples to exclude at import (data integrity issues) ----
samples_to_exclude_import <- c("GSM8828843_a187_es")  # barcode/matrix mismatch in GEO

# ---- Metadata ----
metadata_dir <- file.path(project_root, "data", "sc_reference", "metadata")
metadata_path <- file.path(metadata_dir, paste0("metadata_", dataset_id, ".rds"))

# ---- Seurat creation parameters ----
min_cells <- 3
min_features <- 200

####Script 03_QC_metrics

# ---- QC directories ----
seurat_qc_metrics_dir <- file.path(processed_base_dir, "seurat_qc_metrics")
qc_reports_dir <- file.path(processed_base_dir, "qc_reports")

# ---- Doublet detection ----
run_doublet_detection <- TRUE

# ---- QC gene patterns ----
mito_pattern <- "^MT-"
ribo_pattern <- "^RP[SL]"
hb_pattern   <- "^HB[ABDEG]"

# ---- QC filtering output ----
seurat_qc_filtered_dir <- file.path(processed_base_dir, "seurat_qc_filtered")

# ---- QC thresholds ----
qc_min_features <- 200
qc_max_features <- 7000
qc_min_counts   <- 1000
qc_max_percent_mt <- 20
remove_doublets <- TRUE
