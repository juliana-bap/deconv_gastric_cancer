# Deconvolution project for Immune Signature of gastric cancer patients from the Amazon region

This repository contains the full pipeline for preprocessing, integration and analysis of public scRNA-seq datasets from gastric cancer, with the final objective of performing bulk RNA-seq deconvolution and downstream clinical association analyses.

The project is structured in modular analytical phases to support reproducibility and future publications.

---

## Project Structure

data/
  sc_reference/
    raw/            # Raw scRNA-seq datasets (not versioned)
    processed/      # Processed Seurat objects (not versioned)

scripts/
  sc_pre_proc/      # Import, preprocessing and QC of scRNA-seq datasets
  integration/      # Multi-dataset integration and batch correction
  deconvolution/    # Bulk RNA-seq deconvolution analyses
  clinical_analysis/# Association with clinical variables and outcomes

results/
  figures/
  tables/

docs/               # Documentation and development notes

---

## Analytical Phases

### 1. scRNA-seq Preprocessing (sc_pre_proc)

- Dataset-specific import scripts
- Creation of Seurat objects
- Quality control (QC)
- Normalization
- Cell type annotation
- Generation of reference matrices

Raw data example:
data/sc_reference/raw/GSE246662/

Processed outputs:
data/sc_reference/processed/

---

### 2. Integration (integration)

- Integration of multiple scRNA-seq datasets
- Batch correction
- Construction of unified immune reference
- Export of integrated reference object

---

### 3. Deconvolution (deconvolution)

- Preparation of bulk RNA-seq data
- Application of deconvolution method
- Comparison of estimated immune proportions
- Sensitivity analyses

---

### 4. Clinical Analysis (clinical_analysis)

- Association of deconvolution outputs with:
  - Prognosis
  - Treatment response
  - Molecular subtypes
- Survival analysis
- Statistical modeling

---

## How to Run

### Step 1 – Preprocessing

Run dataset-specific import script located in:

scripts/sc_pre_proc/

Then execute the preprocessing pipeline steps sequentially.

### Step 2 – Integration

Run scripts inside:

scripts/integration/

### Step 3 – Deconvolution

Run scripts inside:

scripts/deconvolution/

### Step 4 – Clinical Analysis

Run scripts inside:

scripts/clinical_analysis/

---

## Dependencies

- R (>= 4.2)
- Seurat
- dplyr
- ggplot2
- survival
- tidyverse

---

## Conventions

- Raw data is never modified.
- Each dataset has its own import script.
- QC thresholds are consistent across datasets.
- Processed objects are saved as `.rds`.
- Data folders are excluded from version control.
- All paths are relative to the project root.

---

## Project Status

- [x] Repository structure defined
- [ ] All scRNA-seq datasets pre-processed
- [ ] Integration completed
- [ ] Reference matrix finalized
- [ ] Deconvolution applied to patient bulk data
- [ ] Clinical association analyses completed