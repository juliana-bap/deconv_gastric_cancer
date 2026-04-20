# ---- Fix orientation from cell x gene to gene x cell ----

fix_orientation <- function(df) {

  genes_in_lines <- any(grepl("^MT-", rownames(df))) |
    any(grepl("^RP[SL]", rownames(df)))

  if (!genes_in_lines) {
    cat("-> Transposed matrix detected\n")
    df <- t(df)
  }

  return(df)

}


# ---- Detect if rownames are Ensembl IDs ----

is_ensembl <- function(ids, threshold = 0.8) {
  mean(grepl("^ENSG[0-9]", ids)) >= threshold
}


# ---- Merge duplicated rownames in sparse matrix (sum counts) ----
# Avoids as.matrix() conversion to prevent memory issues with large matrices

merge_duplicate_rows_sparse <- function(mat) {

  rn <- rownames(mat)
  if (!any(duplicated(rn))) return(mat)

  n_dup <- sum(duplicated(rn))
  cat("Merging", n_dup, "duplicated IDs\n")

  # Group matrix: creates a sparse indicator matrix for aggregation
  groups <- factor(rn, levels = unique(rn))
  group_mat <- Matrix::sparseMatrix(
    i = as.integer(groups),
    j = seq_along(groups),
    x = 1,
    dims = c(nlevels(groups), length(groups)),
    dimnames = list(levels(groups), NULL)
  )

  result <- group_mat %*% mat
  return(result)
}


# ---- Convert gene symbols to Ensembl IDs (bioMart) ----

convert_to_ensembl <- function(df, mapping_table) {

  # Ensure sparse
  if (!inherits(df, "dgCMatrix")) {
    df <- Matrix::Matrix(as.matrix(df), sparse = TRUE)
  }

  # Clean mapping
  mapping_table <- mapping_table[
    mapping_table$hgnc_symbol != "" &
      !duplicated(mapping_table$hgnc_symbol),
  ]

  mapping_table$hgnc_symbol <- toupper(mapping_table$hgnc_symbol)

  # Clean gene symbols
  gene_symbols <- trimws(rownames(df))
  gene_symbols <- toupper(gene_symbols)

  symbol_to_ens <- setNames(
    mapping_table$ensembl_gene_id,
    mapping_table$hgnc_symbol
  )

  ens_ids <- symbol_to_ens[gene_symbols]
  valid <- !is.na(ens_ids)

  cat("Genes before mapping:", length(gene_symbols), "\n")
  cat("Genes mapped:", sum(valid), "\n")

  df <- df[valid, ]
  rownames(df) <- ens_ids[valid]

  # Remove Ensembl version
  rownames(df) <- sub("\\..*", "", rownames(df))

  # Merge duplicated Ensembl IDs (sparse-safe)
  df <- merge_duplicate_rows_sparse(df)

  cat("Genes after deduplication:", nrow(df), "\n")

  return(df)
}


# ---- Clean Ensembl IDs (strip version, deduplicate) ----

clean_ensembl_ids <- function(counts) {

  # Strip version numbers (ENSG00000123456.5 -> ENSG00000123456)
  rownames(counts) <- sub("\\..*", "", rownames(counts))

  # Merge duplicated Ensembl IDs (sparse-safe)
  counts <- merge_duplicate_rows_sparse(counts)

  cat("Genes after deduplication:", nrow(counts), "\n")

  return(counts)
}