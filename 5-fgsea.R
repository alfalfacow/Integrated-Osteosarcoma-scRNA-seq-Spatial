##Hallmark pathway analysis using fgsea
##Setup: Run separately for each malignant scVI cluster in the validation cohort
##Uses cluster-vs-rest differential expression and MSigDB Hallmark gene sets
##Cluster identities: malignant.scvi.clusters
##Note: umap.scvi.malignant is used for visualization only; fgsea uses RNA expression + cluster labels

#A: Data entry and preparation
library(Seurat)
library(dplyr)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(ggplot2)

validation <- readRDS("validation_subset.rds")
DefaultAssay(validation) <- "RNA"

#Join RNA layers for differential expression
validation <- JoinLayers(validation, assay = "RNA")

#Use scVI-derived malignant clusters
Idents(validation) <- "malignant.scvi.clusters"
table(Idents(validation))

#Load Hallmark gene sets
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark_list <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)


#B: Differential expression for fgsea rankings
#Cluster-vs-rest differential expression
validation_markers <- FindAllMarkers(
  validation,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = FALSE,
  min.pct = 0.10,
  logfc.threshold = 0
)

saveRDS(validation_markers, file = "validation_fgsea_rankings.rds")
write.csv(validation_markers, file = "validation_fgsea_rankings.csv", row.names = FALSE)


#C: fgsea
#Rank genes by avg_log2FC for each cluster
validation_fgsea_list <- list()

clusters <- sort(unique(validation_markers$cluster))

for (cl in clusters) {

  mk <- validation_markers[
    validation_markers$cluster == cl,
  ]

  ranks <- mk$avg_log2FC
  names(ranks) <- mk$gene

  #Remove missing/infinite values
  ranks <- ranks[
    is.finite(ranks) &
      !is.na(names(ranks)) &
      names(ranks) != ""
  ]

  #Remove duplicate gene names
  ranks <- ranks[
    !duplicated(names(ranks))
  ]

  #Rank genes from highest to lowest fold change
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgseaMultilevel(
    pathways = hallmark_list,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )

  fg$cluster <- cl

  validation_fgsea_list[[as.character(cl)]] <- fg
}

#Combine all clusters
validation_fgsea <- bind_rows(validation_fgsea_list)


#D: Significant fgsea pathways
#FDR-adjusted p-value < 0.05
validation_fgsea_sig <- validation_fgsea %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

table(validation_fgsea_sig$cluster)

#Save results
saveRDS(validation_fgsea, file = "validation_fgsea_results.rds")

validation_fgsea_csv <- validation_fgsea
validation_fgsea_csv$leadingEdge <- vapply(
  validation_fgsea_csv$leadingEdge,
  paste,
  collapse = ";",
  FUN.VALUE = character(1)
)

write.csv(
  validation_fgsea_csv,
  file = "validation_fgsea_results.csv",
  row.names = FALSE
)


#E: fgsea heatmap
#Only significant pathway-cluster combinations (padj < 0.05)
fg_heat <- validation_fgsea_sig %>%
  select(
    pathway,
    cluster,
    NES
  )

nes_mat <- xtabs(
  NES ~ pathway + cluster,
  data = fg_heat
)

pheatmap(
  nes_mat,
  scale = "none",
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Hallmark fgsea - Validation Malignant States"
)

#Save heatmap
pheatmap(
  nes_mat,
  scale = "none",
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Hallmark fgsea - Validation Malignant States",
  filename = "validation_fgsea_heatmap.pdf",
  width = 9,
  height = 11
)
