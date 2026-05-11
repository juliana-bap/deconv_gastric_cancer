# Pipeline Guide — Gastric Cancer scRNA-seq Pre-processing

This guide walks you through running the scRNA-seq pre-processing pipeline
(Chapter 2, Stage 1) end-to-end on a new dataset. It assumes you have cloned
the repository and completed the [setup](../setup/README.md).

> **Audience:** collaborators and future users of the pipeline. If you just
> want an overview of the project, read the main [README](../README.md) first.

---

## 1. Prerequisites

- **R ≥ 4.4** with all packages installed (`Rscript setup/install_r_deps.R`)
- **Python ≥ 3.10** with the `deconv_gc_py` conda environment (used in step 09b)
- **Quarto ≥ 1.4** (used to render exploration notebooks)
- **RStudio** (optional, but recommended if you want to run notebooks interactively)

See [`setup/README.md`](../setup/README.md) for detailed setup instructions.

---

## 2. Mental model of the pipeline

The pipeline is a sequence of **11 numbered R scripts** (steps 01 → 09c). Every
script reads a single CLI argument: **the path to a dataset-specific config
file**. All inputs, outputs, and parameters live in that config — scripts
never hardcode paths.

```
config_<DATASET>.R   ──┐
                       ▼
01_import  →  02_create_seurat  →  03_qc_metrics  →  04_qc_filtering
                                        ✋ QC notebook
   →  05_merge  →  06_normalize_hvg  →  07a_dimred  →  07b_clustering
          ✋ normalization notebook   ✋ dimred/clustering notebook
   →  08_rogue  →  09a_singler  →  09b_celltypist  →  09c_import_celltypist
      ✋ ROGUE notebook                              ✋ annotation notebook
```

After each ✋ you pause, render (or run) the corresponding Quarto notebook,
inspect the output, and **update the config** with the decisions you made
(thresholds, dimensions, resolutions) before running the next script.

---

## 3. Adding a new dataset (full workflow)

### 3.1. Choose a GEO accession and download the raw data

Download raw files (matrices, barcodes, features, or H5) from GEO and place
them under:

```
data/sc_reference/raw/<DATASET_ID>/
```

This directory is **not versioned**. Keep the folder structure that GEO
provides — your import script will know how to navigate it.

### 3.2. Create the dataset config

Copy the template and rename it:

```bash
cp chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_TEMPLATE.R \
   chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_<DATASET>.R
```

Open the new config and:

1. Fill in the **dataset metadata header** (GEO ID, technology, format, notes).
2. Set `dataset_id <- "<DATASET>"`.
3. Leave the rest of the defaults untouched for now — you will tune them as
   you progress through the pipeline.

> **Note:** `config_<DATASET>.R` is gitignored. Each analyst keeps their own
> local copy. Only `config_TEMPLATE.R` is versioned.

### 3.3. Write the dataset-specific import script

Every dataset on GEO has its own quirks (transposed matrices, mixed samples
in one file, unusual naming, etc.), so step 01 is always dataset-specific.

Create a new file:

```
chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/pipeline/01_import/01_import_<DATASET>.R
```

Use one of the existing import scripts as a starting point:

| Reference script | Good for |
|---|---|
| `01_import_GSE163558.R` | Standard 10X (matrix.mtx + barcodes.tsv + features.tsv) |
| `01_import_GSE246662.R` | CSV files (detects transposition automatically) |
| `01_import_GSE264203.R` | H5 files with multiple samples split by barcode suffix |
| `01_import_GSE291080.R` | 10X with per-sample exclusions |
| `01_import_GSE201347.R` | Pre-built Seurat RDS (large, HPC-ready) |

The contract of any `01_import_*.R` script is:

- Takes **one CLI argument**: the config file path.
- Reads raw data from `raw_base_dir` (defined in the config).
- Writes **one sparse-matrix RDS per sample** into `input_dir`, with clean
  rownames (Ensembl IDs or gene symbols, consistent across samples) and
  barcodes as column names.
- Does **any fixes** needed (transposition, duplicate gene collapse, suffix
  splitting) before writing.

### 3.4. Run the pipeline step by step

From the repo root (so `here::here()` resolves correctly):

```bash
cd ~/doutorado/deconv_gastric_cancer

PIPELINE="chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/pipeline"
CONFIG="chapter_2_sc_deconv_benchmarking/01_sc_matrix/scripts/01_sc_pre_proc/configs/config_<DATASET>.R"
```

Now run each script in order. Below, each block ends with the expected output
and any decision you need to make before moving on.

---

#### Step 01 — Import raw → sparse RDS

```bash
Rscript $PIPELINE/01_import/01_import_<DATASET>.R $CONFIG
```

**Output:** one `.rds` file per sample in `data/sc_reference/raw/<DATASET>/<DATASET>_FIXED/`.

**Check:** open one of the RDS files in R and confirm dimensions are
`genes × cells` (not transposed) and rownames are consistent gene identifiers.

---

#### Step 02 — Build Seurat objects + attach GEO metadata

```bash
Rscript $PIPELINE/02_create_seurat_object.R $CONFIG
```

**Output:** one Seurat object per sample in `data/sc_reference/processed/<DATASET>/seurat_raw/`.

**Check:** `seu$sample_id`, `seu$dataset_id`, and GEO-derived metadata columns
are populated.

---

#### Step 03 — QC metrics (no filtering yet)

```bash
Rscript $PIPELINE/03_qc_metrics.R $CONFIG
```

**Output:** objects with `percent.mt`, `percent.ribo`, `percent.hb`, and
`scDblFinder.class` columns in `seurat_qc_metrics/`.

**✋ Checkpoint:** render the **QC notebook** to inspect per-sample
distributions and decide on thresholds.

See [Section 4](#4-running-the-exploration-notebooks) for how to run notebooks.

**Decisions to record in the config** (before step 04):
- `qc_min_features`, `qc_max_features` — gene count bounds per cell
- `qc_min_counts` — minimum UMIs per cell
- `qc_max_percent_mt` — % mitochondrial cutoff
- `samples_to_exclude_qc` — any sample you decide to drop completely

---

#### Step 04 — Apply QC thresholds

```bash
Rscript $PIPELINE/04_qc_filtering.R $CONFIG
```

**Output:** filtered objects in `seurat_qc_filtered/`.

**Check:** compare pre- vs post-filter cell counts — a report is printed to
the terminal. Excessive loss (>50%) usually means a threshold is too strict.

---

#### Step 05 — Merge all samples

```bash
Rscript $PIPELINE/05_merge.R $CONFIG
```

**Output:** single merged Seurat object at `merged/<DATASET>_merged.rds`.

---

#### Step 06 — Normalize + HVG

```bash
Rscript $PIPELINE/06_normalize_hvg.R $CONFIG
```

**Output:** normalized object at `normalized/<DATASET>_norm_hvg.rds`.

> **Why RC (Relative Counts) without log?** Our downstream goal is a
> deconvolution reference matrix, not differential expression. Log
> transformation will be re-evaluated at the integration step.
> **Do not change `norm_method` without discussing first.**

**✋ Checkpoint:** render the **normalization notebook** to confirm that
library-size normalization looks sensible and HVGs are reasonable.

---

#### Step 07a — PCA

```bash
Rscript $PIPELINE/07a_dimred.R $CONFIG
```

**Output:** object with PCA reduction at `pca/<DATASET>_pca.rds`.

**✋ Checkpoint:** render the **dimred notebook** to inspect the elbow plot.

**Decision to record in the config:**
- `dims_to_use` — which PCs to carry into clustering (usually the elbow)

---

#### Step 07b — Clustering

```bash
Rscript $PIPELINE/07b_clustering.R $CONFIG
```

**Output:** object with clusters at multiple resolutions + UMAP at
`clustered/<DATASET>_clustered.rds`.

**✋ Checkpoint:** render the **clustering notebook** to compare resolutions
side-by-side and choose the working resolution.

**Decision:** set `rogue_resolution` and `annotation_resolution` in the
config (usually the same value).

---

#### Step 08 — ROGUE cluster purity

```bash
Rscript $PIPELINE/08_rogue.R $CONFIG
```

**Output:** ROGUE score matrix (cluster × sample) in `rogue/`.

**✋ Checkpoint:** render the **ROGUE notebook** to flag impure clusters.

---

#### Step 09 — Automated annotation (3 sub-steps)

```bash
Rscript $PIPELINE/09a_pre_annotation_singler.R $CONFIG
Rscript $PIPELINE/09b_run_celltypist.R         $CONFIG   # calls Python via reticulate
Rscript $PIPELINE/09c_import_celltypist.R      $CONFIG
```

**Output:** fully annotated object at `annotation/<DATASET>_annotated.rds`.

**✋ Checkpoint:** render the **annotation notebook** to compare SingleR and
CellTypist predictions per cluster and pick the final label per cluster.

---

## 4. Running the exploration notebooks

Every ✋ checkpoint has a corresponding Quarto notebook template under
`chapter_2_sc_deconv_benchmarking/01_sc_matrix/notebooks/<topic>/`. There are
**two equally valid ways** to run them — pick the one that fits your workflow.

### 4.1. Option A — Render from the terminal (batch mode)

Useful when the dataset is already fully processed and you just want to
generate an HTML snapshot.

```bash
quarto render chapter_2_sc_deconv_benchmarking/01_sc_matrix/notebooks/qc_exploration/00_qc_template.qmd \
  -P dataset_id:<DATASET> \
  --output QC_exploration_<DATASET>.html
```

The `-P dataset_id:<DATASET>` argument overrides the `params$dataset_id`
field in the notebook's YAML, so the notebook knows which config to source.

### 4.2. Option B — Run interactively in RStudio (recommended for exploration)

This is the most comfortable way when you are **still deciding thresholds**.

1. Open the project by double-clicking `deconv_gastric_cancer.Rproj`.
   This sets the working directory to the repo root (required for
   `here::here()` and `source(config_path)` to work).
2. Copy the template to a dataset-specific notebook (this file is gitignored):
   ```bash
   cp chapter_2_sc_deconv_benchmarking/01_sc_matrix/notebooks/qc_exploration/00_qc_template.qmd \
      chapter_2_sc_deconv_benchmarking/01_sc_matrix/notebooks/qc_exploration/QC_exploration_<DATASET>.qmd
   ```
3. Open `QC_exploration_<DATASET>.qmd` in RStudio.
4. Edit **one line** in the YAML header at the top:
   ```yaml
   params:
     dataset_id: "<DATASET>"     # replace GSExxxxxx with your dataset
   ```
5. Run chunk-by-chunk with `Ctrl/Cmd + Enter` or the green play buttons.

> **Does this hurt reproducibility?** No. Reproducibility means that anyone
> re-running the same code with the same config reaches the same result.
> The config (`config_<DATASET>.R`) holds all decisions; the notebook just
> visualises them. Whether you ran it in RStudio or via `quarto render`
> does not change the output.

### 4.3. Available notebooks and when to run each

| Run after script | Notebook template | Decision to make |
|---|---|---|
| 03 | `qc_exploration/00_qc_template.qmd` | QC thresholds, bad samples |
| 06 | `normalization/00_normalization_template.qmd` | Sanity check only |
| 07a | `clustering/00_dimred_template.qmd` | `dims_to_use` |
| 07b | `clustering/00_clustering_template.qmd` | working resolution |
| 08 | `rogue/00_rogue_template.qmd` | Impure clusters to revisit |
| 09c | `annotation/00_annotation_template.qmd` | Final label per cluster |

There is also an exploratory `datasets/` notebook for raw-data inspection
before step 01 (useful for diagnosing transposition, duplicate genes, etc.).

---

## 5. What is versioned vs. local-only

| Type | Path pattern | Versioned? |
|---|---|---|
| Pipeline scripts | `pipeline/*.R`, `pipeline/01_import/*.R` | ✅ yes |
| Config template | `configs/config_TEMPLATE.R` | ✅ yes |
| Notebook templates | `notebooks/**/00_*_template.qmd` | ✅ yes |
| Dataset configs | `configs/config_GSE*.R` | ❌ local |
| Dataset notebooks | `notebooks/**/*_GSE*.qmd` | ❌ local |
| Raw + processed data | `data/` | ❌ local |
| Rendered outputs | `*.html`, `*_files/` | ❌ local |

---

## 6. Troubleshooting

### `Error: config_<DATASET>.R not found`
You are not running from the repo root. `cd` into the repository directory
(or open the `.Rproj` file) before running the scripts or notebooks.

### `Error in rbind(): number of columns does not match`
The per-sample Seurat objects have different metadata columns (e.g. scDblFinder
ran on some samples but not others, or GEO metadata fields differ between
samples). The current QC template uses `dplyr::bind_rows()` which tolerates
this by filling missing columns with `NA`.

If you kept an older dataset-specific copy of the template that still uses
`do.call(rbind, ...)`, replace that chunk with:

```r
combined_meta <- dplyr::bind_rows(lapply(names(seu_list), function(s) {
  df <- seu_list[[s]]@meta.data
  df$sample <- s
  df
}))
```

**Before moving on, it is worth understanding why the columns differ** —
run this diagnostic and see if a sample is missing `scDblFinder.class` or
similar, which would indicate step 03 did not complete for that sample:

```r
sapply(seu_list, function(x) ncol(x@meta.data))
lapply(seu_list, function(x) colnames(x@meta.data))
```

### CellTypist step fails with "Python not found"
You have not set `RETICULATE_PYTHON` in `~/.Renviron`. See
[`setup/README.md`](../setup/README.md#step-3-python-environment).

### `here::here()` returns the wrong directory
You launched R from outside the repo. Either `cd` into the repo and launch
`R`/`Rscript` from there, or open `deconv_gastric_cancer.Rproj` in RStudio.

### Notebook renders but plots are empty
The RDS files expected by the notebook do not exist yet. Re-run the
corresponding pipeline script and verify the output path in the config.

---

## 7. Where to get help

- **Pilot reference report:** `docs/reports/GSE163558_pilot_report.md`
  (the expected format for each dataset report — shows what a "done" dataset
  looks like end-to-end).
- **Config reference:** read `config_TEMPLATE.R` inline comments — every
  parameter is documented there.
- **Author:** Juliana Pinto (juliana.pinto@icb.ufpa.br, jualbup@gmail.com)
