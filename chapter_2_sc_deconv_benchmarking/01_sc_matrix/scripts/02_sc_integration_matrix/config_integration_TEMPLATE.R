# ============================================================================
# CONFIG TEMPLATE — sc_integration_matrix pipeline
# Project: Gastric Cancer Deconvolution
#
# HOW TO USE
#   1. Copy this file to: config_integration.R  (gitignored — stays local)
#   2. Replace every <PLACEHOLDER> with the absolute path on your machine/server.
#   3. Adjust integration parameters below as needed.
#
# CONVENTIONS
#   - All paths are absolute.
#   - This config is sourced by scripts 10–16.
#   - Each script reads only the variables it needs.
# ============================================================================


# ============================================================================
# PROJECT ROOT
# ============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required.")
}
project_root <- here::here()


# ============================================================================
# INPUT: Annotated Seurat objects (output of script 09c per dataset)
# ============================================================================
# Replace <SERVER_PATH_TO_PROJECT> with the absolute project root on the server.
# These files are the final annotated objects from the pre-processing pipeline.

datasets <- c("GSE163558", "GSE275648", "GSE308231", "GSE201347", "GSE249874")

annotated_object_paths <- list(
  GSE163558 = file.path(project_root, "data", "sc_reference", "processed",
                        "GSE163558", "annotation", "GSE163558_annotated.rds"),
  GSE275648 = file.path(project_root, "data", "sc_reference", "processed",
                        "GSE275648", "annotation", "GSE275648_annotated.rds"),
  GSE308231 = file.path(project_root, "data", "sc_reference", "processed",
                        "GSE308231", "annotation", "GSE308231_annotated.rds"),
  GSE201347 = file.path(project_root, "data", "sc_reference", "processed",
                        "GSE201347", "annotation", "GSE201347_annotated.rds"),
  GSE249874 = file.path(project_root, "data", "sc_reference", "processed",
                        "GSE249874", "annotation", "GSE249874_annotated.rds")
)


# ============================================================================
# INPUT: Sample sheet
# ============================================================================

sample_sheet_path <- file.path(
  project_root,
  "chapter_2_sc_deconv_benchmarking", "01_sc_matrix",
  "metadata", "sample_sheet.csv"
)


# ============================================================================
# OUTPUT: Integration directories
# ============================================================================

integration_base_dir <- file.path(
  project_root, "data", "sc_reference", "processed", "integrated"
)

# Script 10 — merged + PCA (unintegrated baseline)
merged_integration_path <- file.path(
  integration_base_dir, "merged", "all_datasets_merged_pca.rds"
)

# Integration benchmark report directory (script 11)
benchmark_output_dir <- file.path(integration_base_dir, "benchmark")

# Script 12 — final integrated object
integrated_seurat_path <- file.path(
  integration_base_dir, "harmony", "all_datasets_harmony.rds"
)

# Script 14 — inferCNV output
infercnv_output_dir <- file.path(integration_base_dir, "infercnv")

# Script 15–16 — annotated + reference matrix
annotated_integrated_path <- file.path(
  integration_base_dir, "annotated", "all_datasets_annotated.rds"
)

reference_matrix_dir <- file.path(integration_base_dir, "reference_matrix")


# ============================================================================
# Script 10 — MERGE + PREPARE
# ============================================================================

# Samples to exclude by disease_group at integration time.
# Intestinal metaplasia is excluded: CDX2+/MUC2+ epithelium is transcriptionally
# similar to intestinal-type GC → would blur normal/malignant boundary.
disease_groups_to_exclude <- c("intestinal_metaplasia")

# Batch-aware HVG selection (SelectIntegrationFeatures)
# Recommended range: 2000–5000. Using 3000 as a compromise between
# capturing biology and limiting batch-specific noise.
n_hvg_integration <- 3000

# PCA: number of components to compute for the unintegrated baseline.
n_pcs_integration <- 50

# Dims to use for UMAP / downstream (set after inspecting elbow plot)
dims_integration <- 1:30


# ============================================================================
# Script 12 — HARMONY INTEGRATION
# ============================================================================

# Batch covariate: Harmony corrects on dataset_id (each dataset = one batch)
harmony_batch_var <- "dataset_id"

# Harmony convergence parameters (defaults are usually fine)
harmony_theta      <- 2    # diversity clustering penalty (higher = more mixing)
harmony_lambda     <- 1    # ridge regression penalty
harmony_max_iter   <- 10   # max correction iterations


# ============================================================================
# END OF CONFIG
# ============================================================================
