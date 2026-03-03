library(org.Hs.eg.db)
library(AnnotationDbi)

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "SYMBOL"
)

mapping <- mapping[!is.na(mapping$ENSEMBL), ]

dir.create("annotation", showWarnings = FALSE)
saveRDS(mapping, "annotation/hgnc_to_ensembl_bioc.rds")
