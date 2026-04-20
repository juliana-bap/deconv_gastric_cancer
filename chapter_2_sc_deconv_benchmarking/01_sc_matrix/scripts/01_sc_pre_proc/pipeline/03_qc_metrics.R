#!/usr/bin/env Rscript

# ==============================
# SCRIPT: 03_qc_metrics.R
# Description: Calculate QC metrics and doublet detection
# Input: Seurat raw objects
# Output: Seurat objects with QC metadata
# ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
})

# ---- Load config ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Provide config file path.")
}

config_path <- normalizePath(args[1])
source(config_path)

# ---- Create output directories ----
dir.create(seurat_qc_metrics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_reports_dir, recursive = TRUE, showWarnings = FALSE)

# ---- List raw Seurat objects ----
seurat_files <- list.files(
  seurat_raw_dir,
  pattern = "_seurat_raw.rds$",
  full.names = TRUE
)

message("Number of Seurat objects found: ", length(seurat_files))

total_cells_pre_qc <- 0

# ---- Storage for summary ----
all_pre_qc_stats <- list()


# ---- Loop over samples ----
for (file in seurat_files) {
  
  samp <- sub("_seurat_raw.rds", "", basename(file))
  message("Processing: ", samp)
  
  seu <- readRDS(file)
  
  # ---- Calculate QC percentages using gene symbols ----
  
  gene_symbols <- seu[["RNA"]]@meta.data$gene_symbol
  
  mito_genes <- rownames(seu)[grepl(mito_pattern, gene_symbols)]
  ribo_genes <- rownames(seu)[grepl(ribo_pattern, gene_symbols)]
  hb_genes   <- rownames(seu)[grepl(hb_pattern, gene_symbols)]
  
  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, features = mito_genes)
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = ribo_genes)
  seu[["percent.hb"]]   <- PercentageFeatureSet(seu, features = hb_genes)
  
  # ---- Doublet detection (optional) ----
  if (run_doublet_detection) {
    
    message("Running scDblFinder...")
    
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce)
    
    seu$scDblFinder.class <- sce$scDblFinder.class
    seu$doublet_score <- sce$scDblFinder.score
    
    message("Doublet classification:")
    print(table(seu$scDblFinder.class))
  }
  
  # ---- Pre-QC statistics ----
  pre_qc_stats <- data.frame(
    Sample = samp,
    number_of_cells = ncol(seu),
    
    mean_counts_per_cell = mean(seu$nCount_RNA),
    median_counts_per_cell = median(seu$nCount_RNA),
    
    mean_features_per_cell = mean(seu$nFeature_RNA),
    median_features_per_cell = median(seu$nFeature_RNA),
    
    mean_percent_mt = mean(seu$percent.mt),
    median_percent_mt = median(seu$percent.mt),
    
    mean_percent_ribo = mean(seu$percent.ribo),
    median_percent_ribo = median(seu$percent.ribo),
    
    mean_percent_hb = mean(seu$percent.hb),
    median_percent_hb = median(seu$percent.hb)
  )
  
  all_pre_qc_stats[[samp]] <- pre_qc_stats
  
  total_cells_pre_qc <- total_cells_pre_qc + ncol(seu)
  
  print(pre_qc_stats)
  
  # ---- Save QC plots automatically ----
  pdf(file.path(qc_reports_dir, paste0(samp, "_QC_metrics.pdf")),
      width = 10, height = 6)
  
  print(VlnPlot(
    seu,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.ribo",
      "percent.hb"
    ),
    ncol = 3
  ))
  
  if (run_doublet_detection) {
    print(VlnPlot(
      seu,
      features = "doublet_score"
    ))
  }
  
  dev.off()
  
  # ---- Save updated Seurat object ----
  saveRDS(
    seu,
    file = file.path(
      seurat_qc_metrics_dir,
      paste0(samp, "_seurat_qc_metrics.rds")
    )
  )
  
  rm(seu)
  gc()
}

# ==============================
# Global summary
# ==============================

summary_df <- do.call(rbind, all_pre_qc_stats)

total_summary <- summary_df[1, ]  # copia estrutura
total_summary[1, ] <- NA          # limpa valores
total_summary$Sample <- "TOTAL"
total_summary$number_of_cells <- sum(summary_df$number_of_cells)

summary_with_total <- rbind(summary_df, total_summary)

# ---- Save CSV summary ----
write.csv(
  summary_with_total,
  file = file.path(seurat_qc_metrics_dir, "QC_summary_pre_filtering.csv"),
  row.names = FALSE
)

# ---- Save log file ----
log_file <- file.path(seurat_qc_metrics_dir, "QC_log.txt")

sink(log_file)

cat("QC METRICS SUMMARY\n")
cat("=========================\n\n")

for (i in 1:nrow(summary_df)) {
  cat("Sample:", summary_df$Sample[i], "\n")
  cat("Number of cells:", summary_df$number_of_cells[i], "\n")
  cat("Mean % MT:", round(summary_df$mean_percent_mt[i], 2), "\n")
  cat("Mean % Ribo:", round(summary_df$mean_percent_ribo[i], 2), "\n")
  cat("Mean % HB:", round(summary_df$mean_percent_hb[i], 2), "\n\n")
}

cat("TOTAL CELLS (ALL SAMPLES):", total_cells_pre_qc, "\n")

sink()

message("Summary and log files saved.")

message("Total cells before QC: ", total_cells_pre_qc)

message("QC metrics calculation completed.")