#Verificando genes não mapeados (Acima de ~15%)

merged = readRDS("/Users/julianapinto/doutorado/deconv_gastric_cancer/data/sc_reference/processed/GSE308231/merged/GSE308231_merged.rds")

# Após o script de merge, no objeto Seurat:
# Ver alguns IDs não mapeados
unmapped <- rownames(merged)[grepl("^ENSG", rownames(merged))]
head(merged, 20)

# Checar se algum marcador clássico ficou de fora
markers <- c("CD3D", "CD8A", "CD4", "CD14", "PTPRC", "EPCAM",
             "PECAM1", "COL1A1", "MS4A1", "NCAM1")
# Esses devem estar presentes como gene symbol, não como ENSG...
markers[!markers %in% rownames(merged)]

#Checar quais não foram mapeados
unmapped <- rownames(merged)[grepl("^ENSG", rownames(merged))]
length(unmapped)
head(unmapped, 30)

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = unmapped,
  mart       = mart
)
table(info$gene_biotype)

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
      filters = "ensembl_gene_id",
      values = c("ENSG00000100604", "ENSG00000184956"),
      mart = mart)