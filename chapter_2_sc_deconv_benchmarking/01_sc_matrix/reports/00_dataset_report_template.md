# Dataset Report — `<DATASET_ID>`

**Author:** Juliana Barreto Albuquerque Pinto
**Date:** YYYY-MM-DD
**Pipeline phase:** sc_pre_proc (pre-integration, pre-annotation)
**Status:** _(in progress / complete)_

> **How to use this template**
> 1. Duplicate this file as `docs/reports/<DATASET_ID>_report.md`
> 2. Fill each section using the values from the dataset's config file and the rendered exploration notebooks
> 3. Keep the same headings — the integration script will read these reports as input for cross-dataset summary tables

---

## 1. Dataset overview

| Field | Value |
|---|---|
| GEO accession | `<GSExxxxxx>` |
| Technology | _(10x Genomics / Smart-seq2 / etc.)_ |
| Reference genome | _(GRCh38 / GRCh37 / other)_ |
| Gene IDs | _(Ensembl / gene symbol — and the source column)_ |
| Original samples | _(N)_ |
| Tissue types | _(primary tumor, normal adjacent, metastases, blood, ...)_ |
| Role in project | _(reference matrix contributor / sub-analysis only)_ |

### Sample exclusion rationale

- **Excluded biologically (N):** _(list samples and reason — e.g. metastases out of scope)_
- **Excluded after QC (N):** _(list samples and the QC metric that disqualified them)_
- **Samples retained (N):** _(final number that goes into the merged object)_

---

## 2. Pipeline parameters

> Copy these from `scripts/sc_pre_proc/configs/config_<DATASET_ID>.R`. Highlight in **bold** any value that differs from the pilot (GSE163558) defaults.

| Step | Parameter | Value |
|---|---|---|
| 02 — Seurat creation | `min_cells` |  |
| 02 — Seurat creation | `min_features` |  |
| 03 — QC metrics | Mito pattern |  |
| 03 — QC metrics | Ribo pattern |  |
| 03 — QC metrics | Hemoglobin pattern |  |
| 03 — QC metrics | Doublet detection |  |
| 04 — QC filtering | `min_features` |  |
| 04 — QC filtering | `max_features` |  |
| 04 — QC filtering | `min_counts` |  |
| 04 — QC filtering | `max_percent_mt` |  |
| 04 — QC filtering | Doublet removal |  |
| 06 — Normalization | Method |  |
| 06 — Normalization | Scale factor |  |
| 06 — HVG | Method |  |
| 06 — HVG | N HVGs |  |
| 07a — PCA | N PCs computed |  |
| 07b — Clustering | PCs used |  |
| 07b — Clustering | Algorithm |  |
| 07b — Clustering | Resolutions tested |  |
| 07b — Clustering | Working resolution |  |
| 09 — Annotation | Resolution |  |
| 09 — Annotation | Methods |  |

### Why these choices

- _(Document any parameter that differs from the pilot defaults and the reason. If everything matches the pilot, write "All parameters follow the GSE163558 pilot defaults.")_

---

## 3. Final dataset metrics

| Metric | Value |
|---|---|
| Cells (post-QC, post-merge) |  |
| Genes |  |
| Samples retained |  |
| Clusters at working resolution |  |
| HVGs |  |
| Annotation columns |  |
| SingleR HPCA labels |  |
| SingleR Blueprint labels |  |
| CellTypist labels |  |
| CellTypist mean confidence |  |
| CellTypist cells with confidence < 0.5 |  |

---

## 4. Cluster-level results (working resolution)

> Fill these from the annotation exploration notebook's "cluster-level consensus" table.

### High-confidence clusters (3 methods agree)

| Cluster | Cells | Identity | Notes |
|---|---|---|---|
|  |  |  |  |

### T cell clusters (subtype granularity differs)

| Cluster | Cells | HPCA | Blueprint | CellTypist | Interpretation |
|---|---|---|---|---|---|
|  |  |  |  |  |  |

### Discordant / problem clusters

| Cluster | Cells | Discordance | Action |
|---|---|---|---|
|  |  |  |  |

---

## 5. Cell composition summary

> Approximate composition as % of total cells.

| Compartment | Clusters | % cells |
|---|---|---|
| Epithelial (tumor + normal) |  |  |
| T cells (CD4 + CD8) |  |  |
| NK cells |  |  |
| B cells / Plasma cells |  |  |
| Myeloid (neutrophils + monocytes/macrophages) |  |  |
| Stroma (fibroblasts + endothelial) |  |  |
| Mast cells |  |  |
| _(other dataset-specific compartments)_ |  |  |

_(2–3 sentences describing whether the composition matches expectations for this tissue/condition.)_

---

## 6. ROGUE — cluster purity

- **Pure clusters** (mean ROGUE > 0.8): _(list)_
- **Moderate** (0.5–0.8): _(list)_
- **Lower purity** (< 0.6): _(list)_

**Decision:** _(re-cluster needed? merge clusters? proceed as-is to integration?)_

---

## 7. Pipeline issues encountered

> Document any technical issue not already in the pilot report. If everything ran cleanly, write "No new issues — pipeline ran without modification."

| Issue | Step | Fix |
|---|---|---|
|  |  |  |

---

## 8. Differences from the pilot dataset

> Brief notes on what is different about this dataset.

- _(e.g. uses CSV format instead of 10x and required orientation fix in script 01)_
- _(e.g. higher mitochondrial content in 2 samples; threshold raised to 25%)_
- _(e.g. only one sample type — no normal adjacent — composition expected to skew tumor)_

---

## 9. Open questions and pending actions for integration

1. _(clusters or annotations to revisit post-integration)_
2. _(samples that may need batch-specific handling)_
3. _(CellTypist confidence concerns)_
4. _(inferCNV decisions specific to this dataset)_

---

## 10. Outputs produced

| Output | Path |
|---|---|
| Raw Seurat objects (per sample) | `data/sc_reference/processed/<DATASET_ID>/seurat_raw/` |
| QC metrics objects | `data/sc_reference/processed/<DATASET_ID>/seurat_qc_metrics/` |
| QC reports (plots) | `data/sc_reference/processed/<DATASET_ID>/qc_reports/` |
| Filtered objects | `data/sc_reference/processed/<DATASET_ID>/seurat_qc_filtered/` |
| Merged object | `data/sc_reference/processed/<DATASET_ID>/merged/<DATASET_ID>_merged.rds` |
| Normalized object | `data/sc_reference/processed/<DATASET_ID>/normalized/<DATASET_ID>_norm_hvg.rds` |
| PCA object | `data/sc_reference/processed/<DATASET_ID>/pca/<DATASET_ID>_pca.rds` |
| Clustered object | `data/sc_reference/processed/<DATASET_ID>/clustered/<DATASET_ID>_clustered.rds` |
| ROGUE scores | `data/sc_reference/processed/<DATASET_ID>/rogue/` |
| Annotated object | `data/sc_reference/processed/<DATASET_ID>/annotation/<DATASET_ID>_annotated.rds` |
| h5ad for CellTypist | `data/sc_reference/processed/<DATASET_ID>/annotation/<DATASET_ID>_for_celltypist.h5ad` |

### Exploration notebooks (rendered HTML)

- `notebooks/datasets/<DATASET_ID>_exploration.html`
- `notebooks/qc_exploration/QC_exploration_<DATASET_ID>.html`
- `notebooks/normalization/norm_exploration_<DATASET_ID>.html`
- `notebooks/clustering/dimred_exploration_<DATASET_ID>.html`
- `notebooks/clustering/clustering_exploration_<DATASET_ID>.html`
- `notebooks/rogue/rogue_exploration_<DATASET_ID>.html`
- `notebooks/annotation/annotation_exploration_<DATASET_ID>.html`

---

## 11. Conclusion

_(2–4 sentences. Did the dataset go through the pipeline cleanly? Is the output ready for integration? What's the main caveat to remember when integrating?)_
