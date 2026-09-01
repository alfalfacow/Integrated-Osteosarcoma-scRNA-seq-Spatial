##Cell-Cell interactions using CellChat
##Setup: One CellChat object for each condition (primary vs metastatic), then compare
##Tutorials:
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#multiple comparisons:
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
#merging is ok:
#https://github.com/jinworks/CellChat/issues/62


library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)


#A: Data entry and preparation

#Split integrated Seurat object by tissue type
metastatic <- subset(Integrated, subset = Tissue_Type == "lung_metastasis")
primary <- subset(Integrated, subset = Tissue_Type == "primary")


#Create metastatic CellChat object
data.input <- GetAssayData(
  metastatic,
  assay = "RNA",
  layer = "data"
) #normalized data matrix

labels <- metastatic$cell2loc_labels
samples <- metastatic$sample

meta <- data.frame(
  labels = labels,
  samples = samples,
  row.names = names(labels)
)

metastatic_CellChat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "labels"
)


#Create primary CellChat object
data.input <- GetAssayData(
  primary,
  assay = "RNA",
  layer = "data"
) #normalized data matrix

labels <- primary$cell2loc_labels
samples <- primary$sample

meta <- data.frame(
  group = labels,
  samples = samples,
  row.names = names(labels)
)

primary_CellChat <- createCellChat(
  object = data.input,
  meta = meta,
  group.by = "group"
)


#B: Load CellChat database

#Load human CellChat database
CellChatDB <- CellChatDB.human

#Look at database categories/interactions
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

#Subset database to protein-mediated signaling categories
CellChatDB.use <- subsetDB(
  CellChatDB,
  search = c(
    "Secreted Signaling",
    "ECM-Receptor",
    "Cell-Cell Contact"
  ),
  key = "annotation"
)

#Assign filtered database to CellChat objects
metastatic_CellChat@DB <- CellChatDB.use
primary_CellChat@DB <- CellChatDB.use

#Subset expression matrices to signaling genes in database
metastatic_CellChat <- subsetData(metastatic_CellChat)
primary_CellChat <- subsetData(primary_CellChat)


#C: Identify overexpressed signaling genes/interactions

metastatic_CellChat <- identifyOverExpressedGenes(
  metastatic_CellChat
)

metastatic_CellChat <- identifyOverExpressedInteractions(
  metastatic_CellChat
)

primary_CellChat <- identifyOverExpressedGenes(
  primary_CellChat
)

primary_CellChat <- identifyOverExpressedInteractions(
  primary_CellChat
)


#D: Compute communication probabilities

#Main analysis uses triMean
metastatic_CellChat <- computeCommunProb(
  metastatic_CellChat,
  type = "triMean"
)

primary_CellChat <- computeCommunProb(
  primary_CellChat,
  type = "triMean"
)

#Compute pathway-level communication probabilities
metastatic_CellChat <- computeCommunProbPathway(
  metastatic_CellChat
)

primary_CellChat <- computeCommunProbPathway(
  primary_CellChat
)

#Filter communication involving cell populations with fewer than 10 cells
metastatic_CellChat <- filterCommunication(
  metastatic_CellChat,
  min.cells = 10
)

primary_CellChat <- filterCommunication(
  primary_CellChat,
  min.cells = 10
)


#E: Extract communication tables

#All inferred ligand-receptor interactions
met_paths <- subsetCommunication(
  metastatic_CellChat
)

prim_paths <- subsetCommunication(
  primary_CellChat
)

#Pathway-level communication tables
met_paths_pathway <- subsetCommunication(
  metastatic_CellChat,
  slot.name = "netP"
)

prim_paths_pathway <- subsetCommunication(
  primary_CellChat,
  slot.name = "netP"
)


#F: Aggregate networks

metastatic_CellChat <- aggregateNet(
  metastatic_CellChat
)

primary_CellChat <- aggregateNet(
  primary_CellChat
)


#G: Initial network plots for each condition

#Metastatic
groupSize <- as.numeric(
  table(metastatic_CellChat@idents)
)

netVisual_circle(
  metastatic_CellChat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength - Metastatic"
)

#Primary
groupSize <- as.numeric(
  table(primary_CellChat@idents)
)

netVisual_circle(
  primary_CellChat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength - Primary"
)


#H: Compute network centrality
#Needed for signaling-role analyses

primary_CellChat <- netAnalysis_computeCentrality(
  primary_CellChat,
  slot.name = "netP"
)

metastatic_CellChat <- netAnalysis_computeCentrality(
  metastatic_CellChat,
  slot.name = "netP"
)


#I: Merge primary and metastatic CellChat objects

object.list <- list(
  Primary = primary_CellChat,
  Metastatic = metastatic_CellChat
)

cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)


#J: Compare total number and strength of interactions

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


#K: Differential interaction networks

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

#Difference in number of interactions
netVisual_diffInteraction(
  cellchat,
  weight.scale = TRUE
)

#Difference in interaction strength
netVisual_diffInteraction(
  cellchat,
  weight.scale = TRUE,
  measure = "weight"
)

par(
  mfrow = c(1, 1)
)


#L: Differential interaction heatmaps

#Number of interactions
ht_count <- netVisual_heatmap(
  cellchat
)

#Interaction strength
ht_weight <- netVisual_heatmap(
  cellchat,
  measure = "weight"
)

ht_count + ht_weight


#M: Primary vs metastatic circle plots using same scale

#Number of interactions
weight.max <- getMaxWeight(
  object.list,
  attribute = c(
    "idents",
    "count"
  )
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
    title.name = paste0(
      "Number of interactions - ",
      names(object.list)[i]
    )
  )
}

par(
  mfrow = c(1, 1)
)


#Interaction strength
weight.max <- getMaxWeight(
  object.list,
  attribute = c(
    "idents",
    "weight"
  )
)

par(
  mfrow = c(1, 2),
  xpd = TRUE
)

for (i in 1:length(object.list)) {

  netVisual_circle(
    object.list[[i]]@net$weight,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    title.name = paste0(
      "Interaction strength - ",
      names(object.list)[i]
    )
  )
}

par(
  mfrow = c(1, 1)
)


#N: Overall outgoing vs incoming signaling roles
#This generated the Primary vs Metastatic scatter plot

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

gg_role <- list()

for (i in 1:length(object.list)) {

  gg_role[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}

patchwork::wrap_plots(
  plots = gg_role
)


#O: Cell-type-specific changes in signaling

#Myeloid
gg_myeloid <- netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "Myeloid Cell"
)

print(gg_myeloid)


#T/NK
gg_TNK <- netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "T/NK Cell"
)

print(gg_TNK)


#Alveolar epithelial cells
gg_lung <- netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "Alveolar epithelial cells (Lung)"
)

print(gg_lung)


#Myoblast
gg_myoblast <- netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "Myoblast"
)

print(gg_myoblast)


#Malignant subcluster 2
gg_sub2 <- netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "Osteoblast_subcluster_2"
)

print(gg_sub2)


#Malignant subcluster 5
gg_sub5 <- netAnalysis_signalingChanges_scatter(
  cellchat,
  idents.use = "Osteoblast_subcluster_5"
)

print(gg_sub5)


#P: Compare signaling pathway information flow

gg_rank <- rankNet(
  cellchat,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)

print(gg_rank)


#Q: Signaling-role heatmaps
#Outgoing, incoming, and overall signaling

pathway.union <- union(
  object.list[[1]]@netP$pathways,
  object.list[[2]]@netP$pathways
)


#Outgoing signaling
ht_out_primary <- netAnalysis_signalingRole_heatmap(
  object.list[[1]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = "Primary"
)

ht_out_metastatic <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = "Metastatic"
)

ComplexHeatmap::draw(
  ht_out_primary + ht_out_metastatic,
  ht_gap = grid::unit(
    0.5,
    "cm"
  )
)


#Incoming signaling
ht_in_primary <- netAnalysis_signalingRole_heatmap(
  object.list[[1]],
  pattern = "incoming",
  signaling = pathway.union,
  title = "Primary"
)

ht_in_metastatic <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "incoming",
  signaling = pathway.union,
  title = "Metastatic"
)

ComplexHeatmap::draw(
  ht_in_primary + ht_in_metastatic,
  ht_gap = grid::unit(
    0.5,
    "cm"
  )
)


#Overall signaling
ht_all_primary <- netAnalysis_signalingRole_heatmap(
  object.list[[1]],
  pattern = "all",
  signaling = pathway.union,
  title = "Primary"
)

ht_all_metastatic <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "all",
  signaling = pathway.union,
  title = "Metastatic"
)

ComplexHeatmap::draw(
  ht_all_primary + ht_all_metastatic,
  ht_gap = grid::unit(
    0.5,
    "cm"
  )
)


#R: Functional signaling similarity

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

rankSimilarity(
  cellchat,
  type = "functional"
)


#S: Pathway-specific follow-up
#Example: VEGF

pathways.show <- c("VEGF")


#VEGF circle plots with matched scale
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
    signaling.name = paste(
      pathways.show,
      names(object.list)[i]
    )
  )
}

par(
  mfrow = c(1, 1)
)


#VEGF heatmaps
ht_vegf <- list()

for (i in 1:length(object.list)) {

  ht_vegf[[i]] <- netVisual_heatmap(
    object.list[[i]],
    signaling = "VEGF",
    color.heatmap = "Reds",
    title.name = paste(
      "VEGF signaling",
      names(object.list)[i]
    )
  )
}

ComplexHeatmap::draw(
  ht_vegf[[1]] + ht_vegf[[2]],
  ht_gap = grid::unit(
    0.5,
    "cm"
  )
)


#T: Save CellChat objects

saveRDS(
  metastatic_CellChat,
  file = "metastatic_CellChat.rds"
)

saveRDS(
  primary_CellChat,
  file = "primary_CellChat.rds"
)

saveRDS(
  cellchat,
  file = "CellChat_Primary_vs_Metastatic_merged.rds"
)
