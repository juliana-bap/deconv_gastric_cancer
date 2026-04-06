# ==============================
# CONFIG - GSE163558
# Project: Gastric Cancer Deconvolution
# Source: GEO
# Technology: 10x Genomics
# Reference genome: GRCh38
# Format: 10x (barcodes/features/matrix per sample dir)
# Gene IDs: Ensembl (from features.tsv column 1)
# Samples: 10
# ==============================

# ---- Dataset identity ----
dataset_id <- "GSE163558"

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

# ---- Samples to exclude (biological criteria) ----
samples_to_exclude_bio <- c("GSM5004184_LN1", "GSM5004185_LN2", "GSM5004187_P1", "GSM5004188_Li1", "GSM5004189_Li2")  # metástases

# ---- Samples to exclude (after QC) ----
samples_to_exclude_qc <- c("GSM5004186_O1") #%MT

#### Script 05_merge

# ---- Merged object output ----
merged_seurat_path <- file.path(processed_base_dir, "merged", paste0(dataset_id, "_merged.rds"))

#### Script 06_normalize_hvg

# ---- Normalization parameters ----
norm_method   <- "RC"       # RC = Relative Counts (CPM without log)
scale_factor  <- 10000

# ---- Highly Variable Genes ----
hvg_method    <- "vst"
n_hvg         <- 2000

# ---- Normalized object output ----
norm_seurat_path <- file.path(processed_base_dir, "normalized", paste0(dataset_id, "_norm_hvg.rds"))

#### Script 07a_dimred

# ---- PCA ----
n_pcs <- 50

# ---- PCA object output ----
pca_seurat_path <- file.path(processed_base_dir, "pca", paste0(dataset_id, "_pca.rds"))

#### Script 07b_clustering

# ---- Neighbors / Clustering / UMAP ----
dims_to_use         <- 1:18
cluster_algorithm   <- 4           # 4 = Leiden
cluster_resolutions <- c(0.25, 0.5, 1.0)

# ---- Clustered object output ----
clustered_seurat_path <- file.path(processed_base_dir, "clustered", paste0(dataset_id, "_clustered.rds"))

#### Script 08_rogue

# ---- ROGUE parameters ----
rogue_resolution  <- 0.5          # resolution to evaluate
rogue_downsample  <- NULL         # NULL = use all cells; set integer (e.g. 1000) for large datasets

# ---- ROGUE output ----
rogue_output_dir <- file.path(processed_base_dir, "rogue")

