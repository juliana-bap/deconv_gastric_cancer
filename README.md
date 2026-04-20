# Deconvolution Project — Immune Signature of Gastric Cancer Patients from the Amazon Region

> **⚠️ Work in progress.** The single-cell preprocessing pipeline is actively running on multiple datasets. Integration, benchmarking, and clinical analyses are planned for subsequent phases. Scripts and notebooks are shared for transparency and reproducibility, not as a finished product.

This repository contains the code and notebooks for a PhD project that builds a gastric cancer immune reference matrix from public scRNA-seq data and uses it to deconvolve bulk RNA-seq samples from Amazon region patients.

---

## PhD Structure

The project is organized in three chapters, each corresponding to a publication:

| Chapter | Goal | Status |
|---------|------|--------|
| `chapter_1_review/` | Review of the immune landscape of gastric cancer | ✍️ Writing |
| `chapter_2_sc_deconv_benchmarking/` | scRNA-seq reference matrix + deconvolution benchmarking | 🔬 In progress |
| `chapter_3_clinical/` | Clinical deconvolution of Amazon patient samples | 🗓️ Planned |

---

## Repository Structure

```
deconv_gastric_cancer/
│
├── chapter_1_review/                    # Review article (manuscript not versioned)
│
├── chapter_2_sc_deconv_benchmarking/
│   ├── 01_sc_matrix/                    # Etapa 1: scRNA-seq preprocessing → reference matrix
│   │   ├── annotation/                  # Ensembl 109 mapping tables (not versioned)
│   │   ├── scripts/
│   │   │   ├── 01_sc_pre_proc/          # Pipeline scripts 01–09
│   │   │   │   ├── pipeline/            # Numbered scripts
│   │   │   │   │   └── 01_import/       # Dataset-specific import scripts
│   │   │   │   └── configs/             # config_TEMPLATE.R (dataset configs are local only)
│   │   │   └── 02_sc_integration_matrix/ # Multi-dataset integration (planned)
│   │   └── notebooks/                   # Quarto exploration notebooks
│   │       ├── datasets/
│   │       ├── qc_exploration/
│   │       ├── normalization/
│   │       ├── clustering/
│   │       ├── rogue/
│   │       └── annotation/
│   └── 02_deconv_benchmarking/          # Etapa 2: tool benchmarking (planned)
│
├── chapter_3_clinical/                  # Clinical application (planned)
│
├── data/                                # NOT versioned — all raw and processed data
│   └── sc_reference/
│       ├── raw/                         # Raw scRNA-seq datasets
│       ├── processed/                   # Processed Seurat objects (.rds)
│       └── metadata/                    # GEO metadata
│
├── results/                             # NOT versioned — figures and tables
├── docs/                                # Project documentation
├── setup/                               # Environment setup scripts
├── environment.yml                      # Conda environment (Python/CellTypist)
└── deconv_gastric_cancer.Rproj
```

---

## Datasets

Five public gastric cancer scRNA-seq datasets from GEO (subset may change):

| GEO ID | Format | Notes | Status |
|--------|--------|-------|--------|
| GSE163558 | 10X (Read10X) | **Pilot** — full pipeline completed | ✅ Done |
| GSE246662 | CSV | Requires matrix orientation fix | 🔬 In progress |
| GSE264203 | H5 | Split by barcode suffix | 🗓️ Pending |
| GSE291080 | 10X (Read10X) | Excludes sample GSM8828843_a187_es | 🗓️ Pending |
| GSE201347 | Seurat RDS | Large dataset, processed on HPC (PBS) | 🔬 In progress |

---

## Preprocessing Pipeline (`01_sc_pre_proc`)

Each script is **modular** and reads all inputs/outputs from a per-dataset config file.

### Running a script

```bash
# From the repo root (required for here::here() to work)
cd ~/doutorado/deconv_gastric_cancer

Rscript chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/pipeline/<script>.R \
        chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_<DATASET>.R
```

> **Note:** always run from within the repository directory so that `here::here()` correctly detects the project root.

### Pipeline steps

| Step | Script | Purpose | Pause for notebook? |
|------|--------|---------|---------------------|
| 00 | `00_utils_general.R` | Shared utility functions | — |
| 00 | `00_download_annotation_ensembl109.R` | Download Ensembl 109 mapping (run once) | — |
| 00 | `00_downloadannotation_org_hs_eg_db.R` | Download org.Hs.eg.db mapping (run once) | — |
| 01 | `01_import/01_import_<DATASET>.R` | Dataset-specific import → sparse matrix RDS | — |
| 02 | `02_create_seurat_object.R` | Build Seurat object with metadata | — |
| 03 | `03_qc_metrics.R` | Compute %mt, %ribo, %hb; run scDblFinder | ✋ QC notebook |
| 04 | `04_qc_filtering.R` / `04_qc_filtering_mad.R` | Apply QC cutoffs | — |
| 05 | `05_merge.R` | Merge samples + JoinLayers | — |
| 06 | `06_normalize_hvg.R` | RC normalization + VST HVG selection | ✋ Normalization notebook |
| 07a | `07a_dimred.R` | ScaleData + PCA | ✋ Dimred notebook |
| 07b | `07b_clustering.R` | FindNeighbors + Leiden + UMAP | ✋ Clustering notebook |
| 08 | `08_rogue.R` | ROGUE cluster purity scores | ✋ ROGUE notebook |
| 09a | `09a_pre_annotation_singler.R` | SingleR (HPCA + Blueprint) + h5ad export | — |
| 09b | `09b_run_celltypist.R` | CellTypist annotation (calls `09b_celltypist.py`) | — |
| 09c | `09c_import_celltypist.R` | Import CellTypist results into Seurat | ✋ Annotation notebook |

### Full example — pilot dataset GSE163558

```bash
PIPELINE="chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/pipeline"
CONFIG="chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_GSE163558.R"

# Import + Seurat object
Rscript $PIPELINE/01_import/01_import_GSE163558.R $CONFIG
Rscript $PIPELINE/02_create_seurat_object.R       $CONFIG

# QC → inspect QC notebook → set thresholds in config → filter
Rscript $PIPELINE/03_qc_metrics.R                 $CONFIG
Rscript $PIPELINE/04_qc_filtering.R               $CONFIG

# Merge + normalization → inspect normalization notebook
Rscript $PIPELINE/05_merge.R                      $CONFIG
Rscript $PIPELINE/06_normalize_hvg.R              $CONFIG

# Dimred → inspect elbow plot → set dims_to_use in config
Rscript $PIPELINE/07a_dimred.R                    $CONFIG
Rscript $PIPELINE/07b_clustering.R                $CONFIG

# Cluster purity → inspect ROGUE notebook
Rscript $PIPELINE/08_rogue.R                      $CONFIG

# Pre-annotation (SingleR + CellTypist)
Rscript $PIPELINE/09a_pre_annotation_singler.R    $CONFIG
Rscript $PIPELINE/09b_run_celltypist.R            $CONFIG
Rscript $PIPELINE/09c_import_celltypist.R         $CONFIG
```

### Adding a new dataset

1. Copy the config template and fill in the dataset-specific fields:
   ```bash
   cp chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_TEMPLATE.R \
      chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_<DATASET>.R
   ```
2. Write a dataset-specific import script under `pipeline/01_import/`.
3. Run the pipeline step by step, pausing at the notebook checkpoints (steps 03, 06, 07a, 07b, 08, 09c).

> Dataset-specific configs are **not versioned** — they contain local paths and dataset-specific thresholds that vary between machines and analysts.

### Exploration notebooks

After each pipeline checkpoint, render the corresponding Quarto notebook to inspect results visually:

```bash
# From the repo root
quarto render chapter_2_sc_deconv_benchmarking/01_sc_matrix/notebooks/qc_exploration/00_qc_template.qmd \
  -P dataset_id:GSE163558
```

Available notebooks: `qc_exploration`, `normalization`, `clustering`, `rogue`, `annotation`, `datasets`.

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/juliana-bap/deconv_gastric_cancer.git
cd deconv_gastric_cancer
```

### 2. R packages

```bash
Rscript setup/install_r_deps.R
```

Installs all CRAN, Bioconductor, and GitHub dependencies idempotently.
See [`setup/README.md`](setup/README.md) for detailed instructions including `renv` lockfile and HPC server setup.

### 3. Python environment (for CellTypist — script 09b)

```bash
conda env create -f environment.yml
conda activate deconv_gc_py
```

Then point reticulate to this environment by adding the following line to `~/.Renviron`:

```
RETICULATE_PYTHON=/path/to/miniconda3/envs/deconv_gc_py/bin/python
```

### 4. Quarto (for notebooks)

Install from <https://quarto.org/docs/get-started/> (version ≥ 1.4).

---

## Dependencies

### R (≥ 4.4)

| Package | Source | Role |
|---------|--------|------|
| Seurat ≥ 5.0 | CRAN | scRNA-seq core |
| Matrix | CRAN | Sparse matrix ops |
| ggplot2, patchwork, pheatmap | CRAN | Visualization |
| dplyr, tibble | CRAN | Data wrangling |
| here | CRAN | Portable path resolution |
| reticulate | CRAN | R ↔ Python bridge |
| GEOquery, biomaRt | Bioconductor | GEO + annotation |
| SingleCellExperiment, scDblFinder | Bioconductor | QC / doublet detection |
| SingleR, celldex | Bioconductor | Automated annotation |
| ROGUE | GitHub (PaulingLiu) | Cluster purity |

### Python (≥ 3.10)

| Package | Role |
|---------|------|
| celltypist ≥ 1.6 | Cell type annotation |
| scanpy, anndata | h5ad I/O |
| scipy, numpy, pandas | Data processing |

### Tested versions (pilot run — GSE163558)

| Tool | Version |
|------|---------|
| R | 4.4.x |
| Seurat | 5.x |
| Python | 3.10 |
| CellTypist | 1.6+ |
| Quarto | 1.7.26 |

---

## Key Design Decisions

- **RC normalization (no log):** downstream goal is a deconvolution reference matrix, not differential expression. Log normalization will be reconsidered at the integration step.
- **Leiden clustering (algorithm 4):** preferred over Louvain for better community detection.
- **Ensembl IDs as rownames:** gene symbols are held in feature metadata; the mapping comes from `01_sc_matrix/annotation/ensembl109_full_mapping.rds`.
- **Portable paths via `here::here()`:** no hardcoded absolute paths in any versioned script. Run R from inside the repo directory.
- **Modular configs:** every decision (thresholds, resolutions, sample exclusions) is documented in a per-dataset config file, not scattered across scripts.

---

## What Is and Is Not Versioned

**Versioned (public):**
- All pipeline scripts (`pipeline/`)
- Config template (`config_TEMPLATE.R`)
- Notebook templates (`00_*_template.qmd`)
- Setup scripts (`setup/`)
- Documentation (`docs/`, `README.md`)

**Not versioned (local only):**
- Raw and processed data (`data/`) — download from GEO
- Seurat objects (`.rds`, `.rda`)
- Rendered notebook outputs (`*.html`, `*.pdf`, `*_files/`)
- Results (`results/`)
- Per-dataset configs (`config_GSE*.R`) — machine and analyst-specific
- Per-dataset notebooks — created from templates for each run
- HPC job scripts (`pbs/`) — cluster-specific

---

## Project Status

- [x] Repository structure defined
- [x] Setup scripts (`install_r_deps.R`, `environment.yml`)
- [x] Full pipeline implemented (scripts 01–09)
- [x] Pilot dataset GSE163558 preprocessed end-to-end
- [ ] Remaining datasets preprocessed (GSE246662, GSE264203, GSE291080, GSE201347)
- [ ] Multi-dataset integration
- [ ] Reference matrix finalized
- [ ] Deconvolution tool benchmarking
- [ ] Deconvolution applied to clinical bulk RNA-seq
- [ ] Clinical association analyses

---

## Author

Juliana Barreto Albuquerque Pinto
PhD candidate — Universidade Federal do Pará (UFPA)
Contact: fcmoreira@ufpa.br
