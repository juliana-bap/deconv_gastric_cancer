# =============================================================================
# SCRIPT: install_r_deps.R
# Project: Gastric Cancer Deconvolution
# Purpose: One-shot installer for all R packages required by the
#          sc_pre_proc pipeline (scripts 00–09).
#
# Usage:
#   Rscript scripts/setup/install_r_deps.R
#
# Notes:
#   - Idempotent: skips packages already installed.
#   - Run this once after cloning the repo, before executing the pipeline.
#   - For full reproducibility, after the first successful pipeline run,
#     initialize renv with `renv::init()` and commit the resulting renv.lock.
# =============================================================================

cran_packages <- c(
  "Seurat",       # >= 5.0
  "Matrix",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "dplyr",
  "tibble",       # required by ROGUE internals
  "here",         # portable project-root detection
  "reticulate",   # R <-> Python bridge for h5ad export
  "BiocManager",
  "remotes",
  "knitr",
  "rmarkdown"
)

bioc_packages <- c(
  "GEOquery",
  "biomaRt",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "SingleCellExperiment",
  "scDblFinder",
  "SingleR",
  "celldex"
)

github_packages <- list(
  "ROGUE" = "PaulingLiu/ROGUE"
)

cat("=============================================\n")
cat(" Installing R dependencies for sc_pre_proc\n")
cat("=============================================\n\n")

# ---- CRAN ----
cat("--- CRAN packages ---\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat("OK :", pkg, "\n")
  }
}

# ---- Bioconductor ----
cat("\n--- Bioconductor packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat("OK :", pkg, "\n")
  }
}

# ---- GitHub ----
cat("\n--- GitHub packages ---\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
for (pkg in names(github_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "from", github_packages[[pkg]], "\n")
    remotes::install_github(github_packages[[pkg]])
  } else {
    cat("OK :", pkg, "\n")
  }
}

cat("\n=============================================\n")
cat(" All dependencies installed\n")
cat("=============================================\n")
cat("R version:", R.version.string, "\n")
cat("Seurat   :", as.character(packageVersion("Seurat")), "\n")
