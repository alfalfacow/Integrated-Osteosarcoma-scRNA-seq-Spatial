##Hallmark pathway analysis using GSVA
##Setup: Average normalized RNA expression within each malignant scVI cluster
##Gene sets: MSigDB Hallmark
##Cluster identities: malignant.scvi.clusters
##Note: GSVA scores are descriptive pathway activity scores; no significance test is performed here

library(Seurat)
library(msigdbr)
library(GSVA)
library(pheatmap)


#A: Data entry and preparation

validation <- readRDS("validation_subset.rds")
DefaultAssay(validation) <- "RNA"

#Join RNA layers
validation <- JoinLayers(
  validation,
  assay = "RNA"
)

#Use scVI-derived malignant clusters
Idents(validation) <- "malignant.scvi.clusters"

table(
  Idents(validation),
  useNA = "ifany"
)


#Load Hallmark gene sets
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark_list <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)


#B: Average normalized expression by malignant cluster

avg_expr_validation <- AverageExpression(
  validation,
  assays = "RNA",
  group.by = "malignant.scvi.clusters",
  layer = "data",
  verbose = FALSE
)$RNA

dim(avg_expr_validation)

#Check gene overlap with Hallmark gene sets
length(
  intersect(
    rownames(avg_expr_validation),
    unique(hallmark_df$gene_symbol)
  )
)


#C: GSVA

gsva_par <- gsvaParam(
  avg_expr_validation,
  hallmark_list,
  minSize = 10,
  maxSize = 500
)

validation_gsva <- gsva(
  gsva_par
)

dim(validation_gsva)
range(validation_gsva)


#D: Save GSVA scores

saveRDS(
  validation_gsva,
  file = "validation_GSVA_scores.rds"
)

write.csv(
  as.data.frame(validation_gsva),
  file = "validation_GSVA_scores.csv",
  row.names = TRUE
)


#E: GSVA heatmap
#Rows are scaled to show relative pathway activity across clusters

pheatmap(
  validation_gsva,
  scale = "row",
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Hallmark GSVA - Validation Malignant States"
)


#Save heatmap

pheatmap(
  validation_gsva,
  scale = "row",
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Hallmark GSVA - Validation Malignant States",
  filename = "validation_GSVA_heatmap.pdf",
  width = 9,
  height = 11
)
