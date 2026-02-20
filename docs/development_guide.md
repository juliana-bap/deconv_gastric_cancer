# Development Guide

This document describes repository usage, update workflow and script header standards.

---

## Repository Usage

### Clone repository

git clone <repository_url>

### Daily workflow

Before starting work:

git pull

After finishing changes:

git add .
git commit -m "Clear and descriptive message"
git push

---

## Commit Guidelines

Use short and descriptive commit messages.

Examples:

Add QC thresholds for sc_reference
Fix gene name formatting in import script
Refactor integration workflow
Update README

Avoid:

update
changes
test

---

## Directory Rules

- Raw data is never modified.
- Raw and processed data are not versioned.
- All paths must be relative to project root.
- No absolute local paths are allowed.
- Intermediate objects must be saved as .rds files.

---

## Script Organization

Scripts are organized by analytical phase:

scripts/sc_pre_proc/
scripts/integration/
scripts/deconvolution/
scripts/clinical_analysis/

Dataset-specific scripts are allowed only in preprocessing phase.

---

## Script Header Standard

All scripts must include a standardized header.

Required fields:

- Script name
- Project name
- Author
- Description
- Input
- Output
- Date
- Dataset (if applicable)

---

## Header Templates

# ==============================
# CONFIG - GSE246662
# Project: Gastric Cancer Deconvolution
# Dataset: GSE246662
# Source: GEO
# Technology: 10x
# Reference genome: GRCh38
# Notes: Some matrices provided transposed
# Created: YYYY-MM-DD
# ==============================


# ==============================
# SCRIPT: import_GSE246662.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE246662
# Description: Import and structural correction of raw matrices
# Input: data/sc_reference/raw/GSE246662/
# Output: data/sc_reference/processed/GSE246662_raw_counts.rds
# Author: Juliana Pinto
# Date: YYYY-MM-DD
# ==============================

# ==============================
# SCRIPT: import_GSE246662.R
# Project: Gastric Cancer Deconvolution
# Phase: sc_pre_proc
# Dataset: GSE246662
# Description: Import and structural correction of raw matrices
# Input: data/sc_reference/raw/GSE246662/
# Output: data/sc_reference/processed/GSE246662_raw_counts.rds
# Author: Juliana Pinto
# Date: YYYY-MM-DD
# ==============================