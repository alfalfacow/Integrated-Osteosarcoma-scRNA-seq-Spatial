##Cell-Cell interactions using CellChat
##Setup: One CellChat object for each condition (primary vs metastatic), then compare
##
##Tutorials:
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#multiple comparisons:
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html


library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(grid)
library(ggplot2)


#A: Load primary and metastatic CellChat objects

primary_CellChat <- readRDS(
  "/Users/elizabethkim/osteosarcoma/cell chat/primary_CellChat.rds"
)

metastatic_CellChat <- readRDS(
  "/Users/elizabethkim/osteosarcoma/cell chat/metastatic_CellChat.rds"
)


#B: Compute network centrality
#Needed for signaling-role analyses

primary_CellChat <- netAnalysis_computeCentrality(
  primary_CellChat,
  slot.name = "netP"
)

metastatic_CellChat <- netAnalysis_computeCentrality(
  metastatic_CellChat,
  slot.name = "netP"
)


#C: Merge primary and metastatic CellChat objects

object.list <- list(
  Primary = primary_CellChat,
  Metastatic = metastatic_CellChat
)

cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

cellchat


#D: Compare total number and strength of interactions

gg1 <- compareInteractions(
  cellchat,
  show.legend = FALSE,
  group = c(1, 2)
)

gg2 <- compareInteractions(
  cellchat,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

gg1 + gg2


#E: Differential interaction networks
#Compare primary vs metastatic interaction number and strength

pdf(
  "differential_interactions.pdf",
  width = 16,
  height = 8
)

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

#Difference in number of interactions
netVisual_diffInteraction(
  cellchat,
  weight.scale = TRUE,
  vertex.label.cex = 0.5
)

#Difference in interaction strength
netVisual_diffInteraction(
  cellchat,
  weight.scale = TRUE,
  measure = "weight",
  vertex.label.cex = 0.5
)

dev.off()


#F: Differential interaction heatmaps

gg1 <- netVisual_heatmap(
  cellchat
)

gg2 <- netVisual_heatmap(
  cellchat,
  measure = "weight"
)

gg1 + gg2


#G: Primary vs metastatic circle plots
#Use the same scale across conditions

weight.max <- getMaxWeight(
  object.list,
  attribute = c(
    "idents",
    "count"
  )
)

pdf(
  "number_of_interactions_primary_vs_metastatic.pdf",
  width = 16,
  height = 8
)

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

for (i in 1:length(object.list)) {

  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    vertex.weight = 8,
    vertex.label.cex = 0.5,
    title.name = paste0(
      "Number of interactions - ",
      names(object.list)[i]
    )
  )
}

dev.off()


#H: Overall outgoing vs incoming signaling roles

num.link <- sapply(
  object.list,
  function(x) {
    rowSums(x@net$count) +
      colSums(x@net$count) -
      diag(x@net$count)
  }
)

weight.MinMax <- c(
  min(num.link),
  max(num.link)
)

gg <- list()

for (i in 1:length(object.list)) {

  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}

patchwork::wrap_plots(
  plots = gg
)


#I: Functional signaling similarity

cellchat <- computeNetSimilarityPairwise(
  cellchat,
  type = "functional"
)

cellchat <- netEmbedding(
  cellchat,
  type = "functional"
)

cellchat <- netClustering(
  cellchat,
  type = "functional"
)

netVisual_embeddingPairwise(
  cellchat,
  type = "functional",
  label.size = 3.5
)


#J: Structural signaling similarity

cellchat <- computeNetSimilarityPairwise(
  cellchat,
  type = "structural"
)

cellchat <- netEmbedding(
  cellchat,
  type = "structural"
)

cellchat <- netClustering(
  cellchat,
  type = "structural"
)

netVisual_embeddingPairwise(
  cellchat,
  type = "structural",
  label.size = 3.5
)


#K: Compare signaling pathway information flow

gg1 <- rankNet(
  cellchat,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)

gg2 <- rankNet(
  cellchat,
  mode = "comparison",
  stacked = FALSE,
  do.stat = TRUE
)

combined <- gg1 + gg2

ggsave(
  "rankNet_primary_vs_metastatic.pdf",
  combined,
  width = 16,
  height = 18
)


#L: Signaling-role heatmaps
#Compare outgoing, incoming, and overall signaling

i <- 1

pathway.union <- union(
  object.list[[i]]@netP$pathways,
  object.list[[i + 1]]@netP$pathways
)


#Outgoing signaling

ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 5,
  height = 12,
  font.size = 6
)

ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i + 1]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[i + 1],
  width = 5,
  height = 12,
  font.size = 6
)

#Shrink pathway labels
ht1@row_names_param$gp <- grid::gpar(
  fontsize = 4
)

ht2@row_names_param$gp <- grid::gpar(
  fontsize = 4
)

#Shrink cell-type labels
ht1@column_names_param$gp <- grid::gpar(
  fontsize = 6
)

ht2@column_names_param$gp <- grid::gpar(
  fontsize = 6
)

pdf(
  "outgoing_signaling_role_heatmap.pdf",
  width = 16,
  height = 22
)

draw(
  ht1 + ht2,
  ht_gap = unit(
    0.5,
    "cm"
  )
)

dev.off()


#Incoming signaling

ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 5,
  height = 12,
  color.heatmap = "GnBu"
)

ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i + 1]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(object.list)[i + 1],
  width = 5,
  height = 12,
  color.heatmap = "GnBu"
)

#Shrink pathway labels
ht1@row_names_param$gp <- grid::gpar(
  fontsize = 4
)

ht2@row_names_param$gp <- grid::gpar(
  fontsize = 4
)

#Shrink cell-type labels
ht1@column_names_param$gp <- grid::gpar(
  fontsize = 6
)

ht2@column_names_param$gp <- grid::gpar(
  fontsize = 6
)

pdf(
  "incoming_signaling_role_heatmap.pdf",
  width = 16,
  height = 22
)

draw(
  ht1 + ht2,
  ht_gap = unit(
    0.5,
    "cm"
  )
)

dev.off()


#Overall signaling

ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "all",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 5,
  height = 12,
  color.heatmap = "OrRd"
)

ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i + 1]],
  pattern = "all",
  signaling = pathway.union,
  title = names(object.list)[i + 1],
  width = 5,
  height = 12,
  color.heatmap = "OrRd"
)

#Shrink pathway labels
ht1@row_names_param$gp <- grid::gpar(
  fontsize = 4
)

ht2@row_names_param$gp <- grid::gpar(
  fontsize = 4
)

#Shrink cell-type labels
ht1@column_names_param$gp <- grid::gpar(
  fontsize = 6
)

ht2@column_names_param$gp <- grid::gpar(
  fontsize = 6
)

pdf(
  "all_signaling_role_heatmap.pdf",
  width = 16,
  height = 22
)

draw(
  ht1 + ht2,
  ht_gap = unit(
    0.5,
    "cm"
  )
)

dev.off()


#M: Cell-type / ligand-receptor comparison
#Example comparison used during follow-up analysis

netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = c(5:11),
  comparison = c(1, 2),
  angle.x = 45
)


#N: MHC-II pathway-specific follow-up

pathways.show <- c("MHC-II")

weight.max <- getMaxWeight(
  object.list,
  slot.name = "netP",
  attribute = pathways.show
)

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

for (i in 1:length(object.list)) {

  netVisual_aggregate(
    object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    vertex.label.cex = 0.6,
    signaling.name = paste(
      pathways.show,
      names(object.list)[i]
    )
  )
}

par(
  mfrow = c(1, 1)
)


#O: Example cell-type-specific pathway follow-up
#CD99 signaling received by myoblasts

#Find myoblast index if needed
myoblast.idx <- which(
  levels(object.list[[1]]@idents) == "Myoblast"
)

gg_cd99_receiver <- netVisual_bubble(
  cellchat,
  sources.use = 1:length(
    levels(object.list[[1]]@idents)
  ),
  targets.use = myoblast.idx,
  signaling = "CD99",
  comparison = c(1, 2),
  angle.x = 45,
  remove.isolate = TRUE
)

ggsave(
  "Myoblast_CD99_receiver.pdf",
  gg_cd99_receiver,
  width = 10,
  height = 7
)


#P: Save merged CellChat comparison object

saveRDS(
  cellchat,
  file = "CellChat_Primary_vs_Metastatic_merged.rds"
)
