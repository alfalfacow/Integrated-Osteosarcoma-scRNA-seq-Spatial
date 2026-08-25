##Cell-Cell interactions using cellchat
##Setup: One cellchat object for each condition (primary vs metastatic), then compare
##Tutorials:
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
#multiple comparisons: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
#merging is ok: https://github.com/jinworks/CellChat/issues/62

library(Seurat)
library(CellChat)

#A: Data Entry and preparation
metastatic <- subset(Integrated, subset = Tissue_Type == "lung_metastasis")
primary <- subset(Integrated, subset = Tissue_Type == "primary")

data.input <- GetAssayData(metastatic, assay = "RNA", layer = "data") # normalized data matrix
labels <- metastatic$cell2loc_labels
samples <- metastatic$sample
meta <- data.frame(labels = labels, samples = samples, row.names = names(labels)) # create a dataframe of the cell labels
metastatic_CellChat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

data.input <- GetAssayData(primary, assay = "RNA", layer = "data") # normalized data matrix
labels <- primary$cell2loc_labels
samples <- primary$sample
meta <- data.frame(group = labels, samples = samples, row.names = names(labels)) # create a dataframe of the cell labels
primary_CellChat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#loading in human CellChat database
CellChatDB <- CellChatDB.human #loads the database
showDatabaseCategory(CellChatDB) #shows type of data in the database
glimpse(CellChatDB$interaction) #shows cool info
#?subsetDB subsetting the cell chat database to only look at the types of interactions we want (in this case we ignore non protein signaling), key is just column for types?
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), key = "annotation")

#Assign our "filtered" database to the cell chat objects
metastatic_CellChat@DB <- CellChatDB.use
primary_CellChat@DB <- CellChatDB.use

#Simplifies object to only include genes that are part of the ligand receptor dataset we set up to save computation cost
metastatic_CellChat <- subsetData(metastatic_CellChat)
primary_CellChat <- subsetData(primary_CellChat)

#B: Getting into the analysis
#Identify highly variable ligand-receptor pairs
metastatic_CellChat <- identifyOverExpressedGenes(metastatic_CellChat)
metastatic_CellChat <- identifyOverExpressedInteractions(metastatic_CellChat) #found 2057 highly variable LR pairs

primary_CellChat <- identifyOverExpressedGenes(primary_CellChat)
primary_CellChat <- identifyOverExpressedInteractions(primary_CellChat) #found 2021 highly variable LR pairs

##Compute communication probability between cell types
#This does the actual probability calculations
metastatic_CellChat <- computeCommunProb(metastatic_CellChat, type="triMean") #ligand-receptor
metastatic_CellChat <- computeCommunProbPathway(metastatic_CellChat) #signaling pathways

primary_CellChat <- computeCommunProb(primary_CellChat, type="triMean") #ligand-receptor
primary_CellChat <- computeCommunProbPathway(primary_CellChat) #signaling pathways

#Filter out communication for low cell count
metastatic_CellChat <- filterCommunication(metastatic_CellChat, min.cells=10)
primary_CellChat <- filterCommunication(primary_CellChat, min.cells=10)

#Extract cell-cell communication in dataframe
met_paths <- subsetCommunication(metastatic_CellChat) #all
prim_paths <- subsetCommunication(primary_CellChat) #all
met_paths <- subsetCommunication(metastatic_CellChat, sources.use = c(12), targets.use = c(1:20)) #1 sending to 5
prim_paths <- subsetCommunication(primary_CellChat, sources.use = c(14, 20), targets.use = c(14, 20)) #1 sending to 5

#Aggregate the object
metastatic_CellChat <- aggregateNet(metastatic_CellChat)
primary_CellChat <- aggregateNet(primary_CellChat)

#Plot circle plots!
groupSize <- as.numeric(table(metastatic_CellChat@idents))
netVisual_circle(metastatic_CellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

groupSize <- as.numeric(table(primary_CellChat@idents))
netVisual_circle(primary_CellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#individual cells circle plot
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#PATHWAY SPECIFIC
pathways.show <- c("TGFb") 
netVisual_aggregate(metastatic_CellChat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(primary_CellChat, signaling = pathways.show, layout = "circle")

#HEATMAP
par(mfrow=c(1,1))
netVisual_heatmap(metastatic_CellChat, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(primary_CellChat, signaling = pathways.show, color.heatmap = "Reds")

#C: redoing with truncatedMeans to check for interferon gamma signaling (only about ~5% T/NK express)
metastatic_CellChatTrunc <- metastatic_CellChat
primary_CellChatTrunc <- primary_CellChat

##Compute communication probability between cell types
#This does the actual probability calculations
metastatic_CellChatTrunc <- computeCommunProb(metastatic_CellChatTrunc, type="truncatedMean", trim = 0.05) #ligand-receptor
metastatic_CellChatTrunc <- computeCommunProbPathway(metastatic_CellChatTrunc) #signaling pathways

primary_CellChatTrunc <- computeCommunProb(primary_CellChatTrunc, type="truncatedMean", trim = 0.05) #ligand-receptor
primary_CellChatTrunc <- computeCommunProbPathway(primary_CellChatTrunc) #signaling pathways

#Filter out communication for low cell count
metastatic_CellChatTrunc <- filterCommunication(metastatic_CellChatTrunc, min.cells=10)
primary_CellChatTrunc <- filterCommunication(primary_CellChatTrunc, min.cells=10)

#aggregate data
metastatic_CellChatTrunc <- aggregateNet(metastatic_CellChatTrunc)
primary_CellChatTrunc <- aggregateNet(primary_CellChatTrunc)

#get pathways
met_paths <- subsetCommunication(metastatic_CellChatTrunc) #all
prim_paths <- subsetCommunication(primary_CellChatTrunc) #all

#heatmap
pathways.show <- c("IFN-II") 
netVisual_heatmap(metastatic_CellChatTrunc, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(primary_CellChatTrunc, signaling = pathways.show, color.heatmap = "Reds")

#C: saving RDS files
saveRDS(metastatic_CellChat, file = "metastatic_CellChat.rds")
saveRDS(primary_CellChat, file = "primary_CellChat.rds")
saveRDS(metastatic_CellChatTrunc, file = "metastatic_CellChat_trunc.rds")
saveRDS(primary_CellChatTrunc, file = "primary_CellChat_trunc.rds")
