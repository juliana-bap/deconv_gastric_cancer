#!/usr/bin/env python3

# ==============================
# SCRIPT: 09b_convert_and_celltypist.py
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Description: Convert annotated Seurat object (h5ad) and run CellTypist.
#              Exports predictions as CSV for import back into R (script 09c).
# Input: h5ad file converted from Seurat (09a_convert_to_h5ad.R)
# Output: CellTypist annotations CSV + UMAP plots
# Usage: python 09b_convert_and_celltypist.py <h5ad_path> <output_dir> <dataset_id>
# Author: Juliana Pinto
# Date: 2026
# ==============================

import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"

import scanpy as sc
import celltypist
from celltypist import models
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for scripts
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 4:
        print("Usage: python 09b_convert_and_celltypist.py <h5ad_path> <output_dir> <dataset_id>")
        sys.exit(1)

    h5ad_path = sys.argv[1]
    output_dir = sys.argv[2]
    dataset_id = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    # ---- Load AnnData ----
    print("=" * 30)
    print(f"Loading: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    # ---- Normalize if not already (CellTypist expects log-normalized) ----
    # Check if data is already normalized by looking at max values
    max_val = adata.X.max()
    if max_val > 50:
        print("Data appears to be raw counts or CPM. Log-normalizing for CellTypist...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(f"Data max value: {max_val:.1f} - appears already log-normalized")

    # ---- Download / load CellTypist model ----
    print("=" * 30)
    print("Loading CellTypist model: Immune_All_Low.pkl")

    model = models.Model.load(model="Immune_All_Low.pkl")
    print(f"Model loaded: {model.description}")

    # ---- Run CellTypist ----
    print("=" * 30)
    print("Running CellTypist annotation...")

    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True
    )

    # Extract predictions
    pred_adata = predictions.to_adata()

    # Add to original adata
    adata.obs["celltypist_label"] = pred_adata.obs.loc[adata.obs.index, "majority_voting"]
    adata.obs["celltypist_conf_score"] = pred_adata.obs.loc[adata.obs.index, "conf_score"]

    print("CellTypist annotation complete")
    print(f"Labels found: {adata.obs['celltypist_label'].nunique()}")
    print("\nLabel counts:")
    print(adata.obs["celltypist_label"].value_counts())

    # ---- Export CSV ----
    print("=" * 30)
    print("Exporting annotations...")

    annotations = adata.obs[["celltypist_label", "celltypist_conf_score"]]
    csv_path = os.path.join(output_dir, f"{dataset_id}_celltypist_annotations.csv")
    annotations.to_csv(csv_path, index=True, index_label="cell_id")
    print(f"CSV saved: {csv_path}")

    # ---- UMAP plots ----
    print("=" * 30)
    print("Generating UMAP plots...")

    # Check if UMAP coordinates exist
    if "X_umap" in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))

        sc.pl.umap(adata, color="celltypist_label", frameon=False,
                    title=f"CellTypist - {dataset_id}", ax=axes[0], show=False)
        sc.pl.umap(adata, color="celltypist_conf_score", frameon=False,
                    title=f"Confidence Score - {dataset_id}", ax=axes[1], show=False)

        umap_path = os.path.join(output_dir, f"{dataset_id}_celltypist_umap.pdf")
        fig.savefig(umap_path, bbox_inches="tight")
        plt.close(fig)
        print(f"UMAP saved: {umap_path}")
    else:
        print("No UMAP coordinates found in h5ad. Skipping UMAP plots.")

    # ---- Done ----
    print("=" * 30)
    print("Done! CellTypist pre-annotation complete")
    print(f"  Labels: {adata.obs['celltypist_label'].nunique()} cell types")
    print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")
    print(f"  CSV: {csv_path}")
    print("Next: run 09c_compare_annotations.R to import and compare")


if __name__ == "__main__":
    main()
