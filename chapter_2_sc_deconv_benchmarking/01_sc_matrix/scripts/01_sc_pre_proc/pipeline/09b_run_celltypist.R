#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 09b_run_celltypist.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Wrapper that reads the config and calls the CellTypist
#              Python script (09b_celltypist.py) with the correct arguments.
#              This keeps the pipeline interface consistent: all scripts
#              are called via Rscript + config.
# Input: Config file (same as all other scripts)
# Output: CellTypist annotations CSV (via Python script)
# Usage: Rscript 09b_run_celltypist.R <config_path>
# Author: Juliana Pinto
# Date: 2026
# ==============================

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript 09b_run_celltypist.R <config_path>\n  No config file provided.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Locate Python script ----
# The Python script lives in the same pipeline directory as this R wrapper.
python_script <- file.path(
  project_root,
  "chapter_2_sc_deconv_benchmarking", "01_sc_matrix",
  "scripts", "01_sc_pre_proc", "pipeline",
  "09b_celltypist.py"
)

if (!file.exists(python_script)) {
  stop("Python script not found: ", python_script)
}

# ---- Check h5ad exists ----
if (!file.exists(h5ad_path)) {
  stop("h5ad file not found: ", h5ad_path, "\n  Run 09a first.")
}

# ---- Run CellTypist ----
cat("\n=============================\n")
cat("Script: 09b_run_celltypist.R (wrapper)\n")
cat("Dataset:", dataset_id, "\n")
cat("Calling Python script:", basename(python_script), "\n")
cat("  h5ad:", h5ad_path, "\n")
cat("  output_dir:", annotation_output_dir, "\n")
cat("  dataset_id:", dataset_id, "\n")
cat("=============================\n\n")

cmd <- sprintf(
  'python3 "%s" "%s" "%s" "%s"',
  python_script,
  h5ad_path,
  annotation_output_dir,
  dataset_id
)

exit_code <- system(cmd)

if (exit_code != 0) {
  stop("CellTypist failed with exit code ", exit_code)
}

cat("\n=============================\n")
cat("CellTypist wrapper complete\n")
cat("Next: Rscript 09c_import_celltypist.R", basename(config_path), "\n")
