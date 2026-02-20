#Função que corrige matriz célula X gene para gene X célula

fix_orientation <- function(df) {
  
  genes_nas_linhas <- any(grepl("^MT-", rownames(df))) |
    any(grepl("^RP[SL]", rownames(df)))
  
  if (!genes_nas_linhas) {
    message("  → Matriz transposta detectada")
    df <- t(df)
  }
  
  return(df)
}