# ==============================
# CONFIG - GSE246662
# Project: Gastric Cancer Deconvolution
# Source: GEO
# Technology: 10x
# Reference genome: GRCh38
# Notes: Some matrices were originally transposed
# ==============================

# ---- Dataset identity ----
dataset_id <- "GSE246662"

# ---- Project root ----
project_root <- "/Users/julianapinto/doutorado/deconv_gastric_cancer"

# ---- Base directories ----
raw_base_dir <- file.path(project_root, "data", "sc_reference", "raw", dataset_id)

processed_base_dir <- file.path(project_root, "data", "sc_reference", "processed", dataset_id)

# ---- Input directory (fixed CSVs from script 01) ----
input_dir <- file.path(raw_base_dir, paste0(dataset_id, "_FIXED"))

# ---- Output directory (Seurat objects) ----
seurat_raw_dir <- file.path(processed_base_dir, "seurat_raw")

# ---- File list ----
csvs <- list.files(
  input_dir,
  pattern = "\\.csv\\.gz$",
  full.names = TRUE
)

# ---- Annotation ----
annotation_dir <- file.path(project_root, "annotation")
mapping_table_path <- file.path(annotation_dir, "ensembl109_full_mapping.rds")

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