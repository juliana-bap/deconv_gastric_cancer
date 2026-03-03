#Fix orientation from cell x gene to gene x cell 

fix_orientation <- function(df) { 
  
  genes_in_lines <- any(grepl("^MT-", rownames(df))) |
    any(grepl("^RP[SL]", rownames(df))) 
  
  if (!genes_in_lines) {
    message(" → Transposed Matrix detected")
    df <- t(df) 
  } 
  
  return(df) 
  
} 


#Convert to ensembl (bioMart)

convert_to_ensembl <- function(df, mapping_table) {
  
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
  
  message("Genes before mapping: ", length(gene_symbols))
  message("Genes mapped: ", sum(valid))
  
  df <- df[valid, ]
  rownames(df) <- ens_ids[valid]
  
  # Remove Ensembl version
  rownames(df) <- sub("\\..*", "", rownames(df))
  
  df <- as.matrix(df)
  mode(df) <- "numeric"
  
  # Merge duplicated Ensembl IDs
  df <- rowsum(df, group = rownames(df))
  
  message("Genes after deduplication: ", nrow(df))
  
  return(df)
}