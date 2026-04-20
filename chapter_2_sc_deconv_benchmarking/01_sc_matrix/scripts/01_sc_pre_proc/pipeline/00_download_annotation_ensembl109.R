# ==============================
# SCRIPT: 00_download_ensembl109_full_mapping.R
# Description: Download full HGNC ↔ Ensembl mapping
# Ensembl version: 109 (GRCh38)
# ==============================

library(biomaRt)

options(timeout = 300)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 109
)

# Download human gene (HGNC symbol)
full_mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  mart = ensembl
)

# Remove empt symbols
full_mapping <- full_mapping[full_mapping$hgnc_symbol != "", ]

# Remover duplicados por símbolo
full_mapping <- full_mapping[!duplicated(full_mapping$hgnc_symbol), ]

dir.create("annotation", showWarnings = FALSE)

saveRDS(full_mapping, "annotation/ensembl109_full_mapping.rds")

cat("Total genes com HGNC:", nrow(full_mapping), "\n")
