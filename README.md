# Deconvolution Project — Immune Signature of Gastric Cancer Patients from the Amazon Region

This repository contains the full pipeline for preprocessing, integration, and analysis of public scRNA-seq datasets from gastric cancer, with the final objective of performing bulk RNA-seq deconvolution and downstream clinical association analyses on patient samples from the Amazon region.

The project is structured in modular analytical phases to support reproducibility and future publications.

---

## Project Structure

```
deconv_gastric_cancer/
├── data/
│   └── sc_reference/
│       ├── raw/            # Raw scRNA-seq datasets (not versioned)
│       ├── processed/      # Processed Seurat objects (not versioned)
│       └── metadata/       # GEO metadata RDS files (not versioned)
│
├── scripts/
│   ├── sc_pre_proc/        # Single-cell preprocessing 
│   │   ├── pipeline/       # Numbered pipeline scripts (00–09)
│   │   │   └── 01_import/  # Dataset-specific import scripts
│   │   ├── configs/        # Per-dataset config files
│   │   └── pbs/            # PBS scripts for HPC server
│   ├── integration/        # Multi-dataset integration 
│   ├── deconvolution/      # Bulk RNA-seq deconvolution
│   └── clinical_analysis/  # Clinical association analyses
│
├── notebooks/              # Quarto exploration notebooks
│   ├── datasets/
│   ├── qc_exploration/
│   ├── normalization/
│   ├── clustering/
│   ├── rogue/
│   └── annotation/
│
├── results/                # Figures, tables (not versioned)
├── docs/                   # Documentation and notes
├── README.md
├── .gitignore
└── deconv_gastric_cancer.Rproj
```

---

## Datasets

Five public gastric cancer scRNA-seq datasets from GEO:

| GEO ID | Format | Notes |
|--------|--------|-------|
| GSE163558 | 10X (Read10X) | **Pilot dataset** — pipeline development |
| GSE246662 | CSV | Requires orientation fix |
| GSE264203 | H5 | Split by barcode suffix |
| GSE291080 | 10X (Read10X) | Excludes sample GSM8828843_a187_es |
| GSE201347 | Seurat RDS | Processed on HPC server (large) |

---

## Preprocessing Pipeline (sc_pre_proc)

Each script is **modular** and reads its inputs/outputs from a per-dataset config file. The standard call is always:

```bash
Rscript <script>.R <config_file>.R
```

### Pipeline steps

| Step | Script | Purpose |
|------|--------|---------|
| 00 | `00_utils_general.R` | Utility functions (orientation fix, Ensembl ID conversion) |
| 00 | `00_download_annotation_ensembl109.R` | Download Ensembl 109 annotation (run once) |
| 00 | `00_downloadannotation_org_hs_eg_db.R` | Download org.Hs.eg.db (run once) |
| 01 | `01_import/01_import_<DATASET>.R` | Dataset-specific raw import → sparse matrix |
| 02 | `02_create_seurat_object.R` | Build Seurat object with gene_symbol, dataset_id, sample_id, GEO metadata |
| 03 | `03_qc_metrics.R` | Compute %mt, %ribo, %hb, run scDblFinder, generate violin plots |
| 04 | `04_qc_filtering.R` *or* `04_qc_filtering_mad.R` | Apply fixed cutoffs (or MAD-based intra-sample) |
| 05 | `05_merge.R` | Merge samples + JoinLayers + rebuild gene_symbol |
| 06 | `06_normalize_hvg.R` | RC normalization + VST HVG selection |
| 07a | `07a_dimred.R` | ScaleData + PCA + ElbowPlot |
| 07b | `07b_clustering.R` | FindNeighbors + Leiden clustering + UMAP |
| 08 | `08_rogue.R` | ROGUE cluster purity scores per cluster × sample |
| 09a | `09a_pre_annotation_singler.R` | SingleR (HPCA + Blueprint) + h5ad export |
| 09b | `09b_run_celltypist.R` | R wrapper that calls `09b_celltypist.py` (CellTypist annotation) |
| 09c | `09c_import_celltypist.R` | Import CellTypist CSV back into Seurat object |

### Running the full pipeline for one dataset

Example for the pilot dataset GSE163558:

```bash
cd scripts/sc_pre_proc/pipeline

# Import + Seurat object
Rscript 01_import/01_import_GSE163558.R ../configs/config_GSE163558.R
Rscript 02_create_seurat_object.R       ../configs/config_GSE163558.R

# QC
Rscript 03_qc_metrics.R                 ../configs/config_GSE163558.R
Rscript 04_qc_filtering.R               ../configs/config_GSE163558.R

# Merge + normalization
Rscript 05_merge.R                      ../configs/config_GSE163558.R
Rscript 06_normalize_hvg.R              ../configs/config_GSE163558.R

# Dimred + clustering
Rscript 07a_dimred.R                    ../configs/config_GSE163558.R
Rscript 07b_clustering.R                ../configs/config_GSE163558.R

# Cluster purity
Rscript 08_rogue.R                      ../configs/config_GSE163558.R

# Pre-annotation (SingleR + CellTypist)
Rscript 09a_pre_annotation_singler.R    ../configs/config_GSE163558.R
Rscript 09b_run_celltypist.R            ../configs/config_GSE163558.R
Rscript 09c_import_celltypist.R         ../configs/config_GSE163558.R
```

### Exploration notebooks

After running the pipeline, render the corresponding Quarto notebook to inspect the results:

```bash
cd notebooks/<topic>
quarto render <notebook>_GSE163558.qmd
```

Available notebook topics: `qc_exploration`, `normalization`, `clustering`, `rogue`, `annotation`.

---

## Setup

A one-shot setup script and a conda environment file are provided:

```bash
# 1. R packages (CRAN + Bioconductor + GitHub)
Rscript scripts/setup/install_r_deps.R

# 2. Python environment (CellTypist)
conda env create -f environment.yml
conda activate deconv_gc_py
```

See [`scripts/setup/README.md`](scripts/setup/README.md) for full setup instructions, including how to point reticulate at the conda env and how to lock R package versions with `renv`.

---

## Dependencies

### R (≥ 4.4)

CRAN packages:

```r
install.packages(c(
  "Seurat",         # ≥ 5.0  — scRNA-seq core
  "Matrix",         # sparse matrix ops
  "ggplot2",
  "patchwork",
  "pheatmap",
  "dplyr",
  "tibble",         # required by ROGUE
  "reticulate",     # R ↔ Python bridge for h5ad export
  "remotes"
))
```

Bioconductor packages:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery",
  "biomaRt",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "SingleCellExperiment",
  "scDblFinder",
  "SingleR",
  "celldex"
))
```

GitHub-only packages:

```r
remotes::install_github("PaulingLiu/ROGUE")
```

### Python (≥ 3.10) — for CellTypist (script 09b)

```bash
pip install scanpy celltypist anndata pandas matplotlib scipy numpy
```

Or via conda:

```bash
conda create -n celltypist_env python=3.10
conda activate celltypist_env
pip install scanpy celltypist anndata pandas matplotlib
```

The R script `09a_pre_annotation_singler.R` uses **reticulate** to write h5ad directly (avoiding sceasy, which is incompatible with Seurat v5). Make sure reticulate sees the right Python environment — set `RETICULATE_PYTHON` in `~/.Renviron` or in the dataset config if needed.

### Quarto (≥ 1.4)

For rendering the exploration notebooks. Install from <https://quarto.org/docs/get-started/>.

### Recommended tested versions (pilot run)

| Tool | Version |
|------|---------|
| R | 4.4.x |
| Seurat | 5.x |
| Python | 3.10 |
| CellTypist | 1.6+ |
| Quarto | 1.7.26 |

---

## Working with the GitHub repository

### Clone

```bash
git clone https://github.com/juliana-bap/deconv_gastric_cancer.git
cd deconv_gastric_cancer
```

### Branch workflow

The default branch is `main`. For new work:

```bash
git checkout -b feature/<short-description>
# ... make changes ...
git add scripts/<file>.R
git commit -m "Short message describing the change"
git push -u origin feature/<short-description>
```

Then open a Pull Request on GitHub against `main`.

### What is **not** versioned (see `.gitignore`)

- All raw and processed data (`data/raw/`, `data/sc_reference/raw/`, `data/sc_reference/processed/`)
- All `.rds`/`.rda` Seurat objects
- All Quarto rendered outputs (`*.html`, `*.pdf`, `*.png`, `*_files/`)
- The `results/` folder
- R session files (`.Rhistory`, `.RData`, `.Rproj.user/`)
- macOS metadata (`.DS_Store`)
- `.claude/` (local environment)

### What **is** versioned

- All scripts (`scripts/`)
- All notebook source `.qmd` files
- Per-dataset config files
- Documentation (`README.md`, `docs/`)

---

## Conventions

- **Raw data is never modified.** Each dataset has its own import script.
- **Modular scripts:** every script reads inputs and writes outputs based on a config file (`scripts/sc_pre_proc/configs/config_<DATASET>.R`).
- **Consistent interface:** all scripts are called as `Rscript <script>.R <config>.R`.
- QC thresholds are defined per-dataset in the config, not hard-coded.
- Processed objects are saved as `.rds` (excluded from git).
- All paths are absolute via `project_root` defined in each config.
- Each dataset has its own exploration notebooks under `notebooks/<topic>/`.

---

## Project Status

- [x] Repository structure defined
- [x] Pilot dataset (GSE163558) preprocessed end-to-end (steps 01–09)
- [ ] All scRNA-seq datasets pre-processed
- [ ] Multi-dataset integration completed
- [ ] Reference matrix finalized
- [ ] Deconvolution applied to patient bulk data
- [ ] Clinical association analyses completed

---

## Author

Juliana Barreto Albuquerque Pinto — PhD Project, Gastric Cancer Deconvolution
