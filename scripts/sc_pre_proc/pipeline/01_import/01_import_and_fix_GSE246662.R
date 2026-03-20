# ==============================
# SCRIPT: import_GSE246662.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE246662
# Description: Import, orientation fix and gene symbol to Ensembl conversion
# Input: data/sc_reference/raw/GSE246662/
# Output: data/sc_reference/raw/GSE246662/GSE246662_FIXED/
# Author: Juliana Pinto
# Date: 23-02-2026
# Ensembl version: 109 (GRCh38)
# ==============================

# Packages

library(biomaRt)

# Files

csvs <- c(
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874169_HL1.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874170_HL2.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874171_HL3.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874172_GC1.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874173_GC2.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874174_GC3.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874175_LM1.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874176_LM2.csv.gz",
  "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSM7874177_LM3.csv.gz"
)

# Functions

source("/Users/julianapinto/doutorado/deconv_gastric_cancer/scripts/sc_pre_proc/pipeline/00_utils_general.R")

mapping_table <- readRDS("annotation/ensembl109_full_mapping.rds")

# Output directory

outdir <- "/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/raw/GSE246662/GSE246662_FIXED"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Main loop

for (path in csvs) {
  
  cat("\n=============================\n")
  cat("Sample:", basename(path), "\n")
  
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  
  # Fix orientation
  df <- fix_orientation(df)
  cat("Dimension after orientation:", dim(df), "\n")
  
  # Convert to Ensembl
  df <- convert_to_ensembl(df, mapping_table)
  
  cat("Final dimension (Ensembl):", dim(df), "\n")
  
  # Save counts
  outfile <- file.path(outdir, basename(path))
  write.csv(df, gzfile(outfile))
  
  rm(df)
  gc()
}

