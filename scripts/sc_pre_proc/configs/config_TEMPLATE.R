# ============================================================================
# CONFIG TEMPLATE — sc_pre_proc pipeline
# Project: Gastric Cancer Deconvolution
#
# HOW TO USE
#   1. Copy this file to: config_<DATASET_ID>.R  (e.g. config_GSE246662.R)
#   2. Replace every <PLACEHOLDER> with the value for your dataset.
#   3. Adjust parameters per script as you progress through the pipeline.
#      You don't need to set everything at once — fill in the variables for
#      each script just before running it. The pipeline is sequential.
#   4. Source check (optional): `source("config_<DATASET>.R")` in R should
#      run without errors before you start the pipeline.
#
# CONVENTIONS
#   - All paths are absolute, derived from `project_root`.
#   - Each script reads only the variables it needs from this file.
#   - Values that differ from the GSE163558 pilot defaults should be
#     documented in the dataset's report (docs/reports/<DATASET>_report.md).
# ============================================================================


# ============================================================================
# DATASET METADATA (header — fill before anything else)
# ============================================================================
# Update the comment block below with what is known about the dataset.
# This is not parsed by R; it exists so anyone reading the config knows
# the basics at a glance.
#
# GEO accession    : <GSExxxxxx>
# Source           : GEO / ArrayExpress / other
# Technology       : 10x Genomics / Smart-seq2 / DropSeq / etc.
# Reference genome : GRCh38 / GRCh37 / other
# Format           : 10x (mtx) / CSV / H5 / Seurat RDS
# Gene IDs         : Ensembl / gene symbol (and source column if applicable)
# Original samples : <N>
# Notes            : <anything unusual: transposed matrices, suffix splits,
#                    metastatic samples to drop, etc.>
# ============================================================================


# ---- Dataset identity ------------------------------------------------------
dataset_id <- "<GSExxxxxx>"


# ---- Project root ----------------------------------------------------------
# Absolute path to the repository root. Update for your machine if needed.
project_root <- "/Users/julianapinto/doutorado/deconv_gastric_cancer"


# ---- Base directories (do not edit unless project layout changes) ----------
raw_base_dir       <- file.path(project_root, "data", "sc_reference", "raw",       dataset_id)
processed_base_dir <- file.path(project_root, "data", "sc_reference", "processed", dataset_id)


# ============================================================================
# Script 01 — IMPORT
# ----------------------------------------------------------------------------
# Each dataset has its own import script under pipeline/01_import/.
# That script's job is to read the raw data, fix any orientation/format
# issues, and write one sparse RDS per sample into `input_dir`.
# ============================================================================

# Where the cleaned, sparse-matrix RDS files (one per sample) will be written
# by the dataset-specific 01_import_<DATASET>.R script.
input_dir <- file.path(raw_base_dir, paste0(dataset_id, "_FIXED"))

# Annotation files (shared across datasets — usually no need to change)
annotation_dir     <- file.path(project_root, "annotation")
mapping_table_path <- file.path(annotation_dir, "ensembl109_full_mapping.rds")

# Utility functions (shared)
utils_path <- file.path(project_root, "scripts", "sc_pre_proc", "pipeline", "00_utils_general.R")

# GEO metadata RDS (downloaded once via 01_download_metadata script when
# needed). For most datasets it's downloaded automatically inside 01_import.
metadata_dir  <- file.path(project_root, "data", "sc_reference", "metadata")
metadata_path <- file.path(metadata_dir, paste0("metadata_", dataset_id, ".rds"))


# ============================================================================
# Script 02 — CREATE SEURAT OBJECT
# ----------------------------------------------------------------------------
# Reads the per-sample RDS files in `input_dir` and builds Seurat objects
# with gene_symbol, dataset_id, sample_id, and embedded GEO metadata.
# ============================================================================

# Output directory for raw (unfiltered) Seurat objects
seurat_raw_dir <- file.path(processed_base_dir, "seurat_raw")

# List of input RDS files (auto-populated from input_dir)
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

# Seurat::CreateSeuratObject parameters
# These are very loose initial filters — real QC happens in scripts 03–04.
min_cells    <- 3      # gene must be detected in at least N cells
min_features <- 200    # cell must have at least N detected genes


# ============================================================================
# Script 03 — QC METRICS
# ----------------------------------------------------------------------------
# Computes %mt, %ribo, %hb, runs scDblFinder, and produces violin/histogram
# plots saved to `qc_reports_dir`. No filtering happens here.
#
# After running this, render the QC exploration notebook
# (notebooks/qc_exploration/00_qc_template.qmd) to inspect the metrics
# and decide on the cutoffs for script 04.
# ============================================================================

# Output directories
seurat_qc_metrics_dir <- file.path(processed_base_dir, "seurat_qc_metrics")
qc_reports_dir        <- file.path(processed_base_dir, "qc_reports")

# Doublet detection (scDblFinder). Set to FALSE only if the dataset is
# extremely small or already deduplicated upstream.
run_doublet_detection <- TRUE

# Gene-name patterns used to compute %mt / %ribo / %hb metrics.
# Defaults work for human gene symbols (HGNC). For datasets where rownames
# are Ensembl IDs, the pipeline maps via gene_symbol metadata internally —
# you should not need to change these.
mito_pattern <- "^MT-"
ribo_pattern <- "^RP[SL]"
hb_pattern   <- "^HB[ABDEG]"


# ============================================================================
# Script 04 — QC FILTERING
# ----------------------------------------------------------------------------
# Applies fixed cutoffs (and optionally removes doublets and excluded
# samples). There are two variants of this script:
#   - 04_qc_filtering.R       (fixed cutoffs only)
#   - 04_qc_filtering_mad.R   (fixed cutoffs + MAD-based intra-sample)
#
# Set the cutoffs below AFTER inspecting the QC notebook.
# ============================================================================

# Output directory
seurat_qc_filtered_dir <- file.path(processed_base_dir, "seurat_qc_filtered")

# QC thresholds — review these per dataset using the QC notebook.
# The defaults below are the GSE163558 pilot values. Lower bounds protect
# against empty droplets / low-quality cells; upper bounds catch likely
# doublets or aggregates not caught by scDblFinder.
qc_min_features    <- 200      # min genes detected per cell
qc_max_features    <- 7000     # max genes (often empirical, dataset-specific)
qc_min_counts      <- 1000     # min total UMIs per cell
qc_max_percent_mt  <- 20       # max % mitochondrial reads (often 10–25%)
remove_doublets    <- TRUE     # remove cells flagged by scDblFinder

# ---- Sample exclusion ------------------------------------------------------
# Samples to drop for biological reasons (e.g. metastases when only the
# primary tumor reference is being built). Use the GSM identifiers as they
# appear in seurat_raw/.
samples_to_exclude_bio <- c(
  # "GSMxxxxxxx_LN1",
  # "GSMxxxxxxx_Li1"
)

# Samples to drop after QC inspection (e.g. very high MT or low cell count).
# Document the reason in an inline comment so it's auditable later.
samples_to_exclude_qc <- c(
  # "GSMxxxxxxx_O1"  # rejected: median %MT > 30
)


# ============================================================================
# Script 05 — MERGE
# ----------------------------------------------------------------------------
# Merges all retained samples into a single Seurat object, joins layers,
# rebuilds the gene_symbol feature metadata, and adds the immunoglobulin
# supplement (used downstream for B/plasma cell handling).
# ============================================================================

merged_seurat_path <- file.path(processed_base_dir, "merged",
                                paste0(dataset_id, "_merged.rds"))


# ============================================================================
# Script 06 — NORMALIZATION + HVG
# ----------------------------------------------------------------------------
# Normalizes counts and selects highly variable genes.
#
# We use Relative Counts (RC) without log because the downstream goal is a
# deconvolution reference matrix. Normalization will be reconsidered at the
# integration step. Don't change these defaults without discussing first.
# ============================================================================

# Normalization
norm_method  <- "RC"        # Relative Counts (CPM-like, no log)
scale_factor <- 10000       # standard 10K library size

# Highly Variable Genes
hvg_method <- "vst"         # Seurat default
n_hvg      <- 2000          # 2000 is the conventional starting point

# Output
norm_seurat_path <- file.path(processed_base_dir, "normalized",
                              paste0(dataset_id, "_norm_hvg.rds"))


# ============================================================================
# Script 07a — DIMRED (PCA only)
# ----------------------------------------------------------------------------
# Runs ScaleData + PCA + ElbowPlot. Output is reused by 07b — do not skip.
#
# After running this, render the dimred notebook
# (notebooks/clustering/00_dimred_template.qmd) to inspect the elbow plot
# and decide how many PCs to use in 07b.
# ============================================================================

n_pcs <- 50    # number of PCs to compute (we keep 50, then subset in 07b)

pca_seurat_path <- file.path(processed_base_dir, "pca",
                             paste0(dataset_id, "_pca.rds"))


# ============================================================================
# Script 07b — CLUSTERING
# ----------------------------------------------------------------------------
# FindNeighbors + Leiden clustering + UMAP. The PCs and resolutions chosen
# here directly affect the cluster structure that ROGUE and the annotation
# scripts will see.
#
# After running this, render the clustering notebook
# (notebooks/clustering/00_clustering_template.qmd).
# ============================================================================

# PCs to feed into FindNeighbors / RunUMAP. Set this from the dimred notebook.
dims_to_use <- 1:18

# Clustering algorithm: 1=Louvain, 2=Louvain (refined), 3=SLM, 4=Leiden.
# Leiden (4) is preferred — better community detection.
cluster_algorithm <- 4

# Resolutions to compute. Leave several so the clustering notebook can
# compare them. The "working" resolution is set below for downstream scripts.
cluster_resolutions <- c(0.25, 0.5, 1.0)

# Output
clustered_seurat_path <- file.path(processed_base_dir, "clustered",
                                   paste0(dataset_id, "_clustered.rds"))


# ============================================================================
# Script 08 — ROGUE (cluster purity)
# ----------------------------------------------------------------------------
# Computes ROGUE scores per cluster x sample. Use these to identify clusters
# that may be biologically heterogeneous.
# ============================================================================

# Resolution to evaluate (usually the working resolution chosen above)
rogue_resolution <- 0.5

# Downsampling: NULL = use all cells (recommended for moderate datasets).
# For very large datasets (>100k cells) set an integer like 1000 to keep
# the computation tractable.
rogue_downsample <- NULL

rogue_output_dir <- file.path(processed_base_dir, "rogue")


# ============================================================================
# Scripts 09a / 09b / 09c — PRE-ANNOTATION
# ----------------------------------------------------------------------------
# 09a: SingleR (HPCA + Blueprint) + h5ad export for CellTypist
# 09b: R wrapper that calls 09b_celltypist.py
# 09c: import CellTypist CSV back into Seurat
#
# After running all three, render the annotation notebook
# (notebooks/annotation/00_annotation_template.qmd).
# ============================================================================

# Working resolution for annotation (usually same as rogue_resolution)
annotation_resolution <- 0.5

# Output paths
annotation_output_dir <- file.path(processed_base_dir, "annotation")

annotated_seurat_path <- file.path(annotation_output_dir,
                                   paste0(dataset_id, "_annotated.rds"))

h5ad_path <- file.path(annotation_output_dir,
                       paste0(dataset_id, "_for_celltypist.h5ad"))


# ============================================================================
# END OF CONFIG
# ----------------------------------------------------------------------------
# Quick self-check: source this file in R. If you get errors about missing
# objects, you have not yet filled in a section that a downstream script
# needs. That's fine — fill them in as you progress through the pipeline.
# ============================================================================
