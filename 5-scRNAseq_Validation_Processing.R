##scRNAseq validation preprocessing
Val1 <- readRDS("/Users/alfred/Downloads/Archive/GSM9504297_with_scevan.rds")
Val2 <- readRDS("/Users/alfred/Downloads/Archive/GSM9504298_with_scevan.rds")
Val3 <- readRDS("/Users/alfred/Downloads/Archive/GSM9504303_with_scevan.rds")
Val4 <- readRDS("/Users/alfred/Downloads/Archive/GSM9504304_with_scevan.rds")
validation <- merge(Val1, c(Val2, Val3, Val4))
validation <- NormalizeData(validation)
validation <- FindVariableFeatures(validation)
validation <- ScaleData(validation)
validation <- RunPCA(validation)

#visualization of unintegrated clusters
validation <- FindNeighbors(validation, dims = 1:20, reduction = "pca")
validation <- FindClusters(validation, resolution = 0.5, cluster.name = "unintegrated_clusters")
validation <- RunUMAP(validation, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(validation, reduction = "umap.unintegrated", group.by = "unintegrated_clusters") 
DimPlot(validation, reduction = "umap.unintegrated", group.by = "orig.ident")
DimPlot(validation, reduction = "umap.unintegrated", group.by = "sample") #not very well integrated...
DimPlot(validation, reduction = "umap.unintegrated", group.by = "Tissue_Type")
DimPlot(validation, reduction = "umap.unintegrated", group.by = "SCEVAN.Malignancy")
#integration, begins at train_loss = 710
validation <- IntegrateLayers(
  object          = validation,
  method          = scVIIntegration,
  new.reduction   = "integrated.scvi",
  conda_env       = "/Users/alfred/micromamba/envs/scvi_env",
  features        = Features(validation),   # all genes, matching their unfiltered adata, null = variablefeatures
  layers          = "counts",
  groups.name     = "sample",
  groups          = validation[[]],      
  ndims.out       = 30,      # matches their n_latent=30
  n_hidden        = 256,     # matches
  n_layers        = 2,       # matches
  dropout_rate    = 0.1,     # matches
  max_epochs      = 200,     # matches
  gene_likelihood = "zinb",  # matches their unset default (scvi-tools defaults to zinb)
  seed.use        = 0,
  early_stopping = TRUE,
  early_stopping_patience = 30,
  early_stopping_monitor = "elbo_validation",
  torch.intraop.threads = 1L,
  torch.interop.threads = 1L,
  verbose = TRUE
)

validation <- FindNeighbors(validation, reduction = "integrated.scvi", dims = 1:30)
validation <- FindClusters(validation, resolution = 0.8, cluster.name = "scvi_clusters") #harrison's script uses 0.8
validation <- RunUMAP(validation, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
DimPlot(validation, reduction = "umap.scvi", group.by = "sample", combine = FALSE, label.size = 2)
DimPlot(validation, reduction = "umap.scvi", group.by = "Tissue_Type", combine = FALSE, label.size = 2)
a<-DimPlot(validation, reduction = "umap.scvi", group.by = "scvi_clusters", combine = FALSE, label.size = 2, label = TRUE)
b<-DimPlot(validation, reduction = "umap.scvi", group.by = "SCEVAN.Malignancy", combine = FALSE, label.size = 2)

#saveRDS(validation, "validation_post_scvi.rds")

FeaturePlot(validation, features = c("COL1A1", "CDH11", "ITGA10", "RUNX2", "IBSP", "FGFR3"), reduction = "umap.scvi") #osteoblasts
FeaturePlot(validation, features = c("COL1A1", "CDH11", "RUNX2", "IBSP", "CLEC11A", "ALPL", "SATB2"), reduction = "umap.scvi") #osteoblasts
FeaturePlot(validation, features = c("CTSK", "MMP9", "ACP5"), reduction = "umap.scvi") #osteoclasts
FeaturePlot(validation, features = c("IL7R", "CD3D", "NKG7", "CD3E", "GZMK", "GZMA", "TRAC", "CD3G", "TRBC1"), reduction = "umap.scvi") #T/NK Cells
FeaturePlot(validation, features = c("CD74", "CD14", "FCGR3A", "LYZ", "CD163", "MS4A7", "S100A9", "APOE", "C1QA", "CD68"), reduction = "umap.scvi") #Myeloid Cells
FeaturePlot(validation, features = c("COL1A1", "LUM", "DCN", "FBLN1", "ACTA2", "TAGLN", "COL3A1", "COL6A1", "FAP"), reduction = "umap.scvi") #Fibroblast
FeaturePlot(validation, features = c("PECAM1", "VWF", "EGFL7", "CDH5", "CLDN5", "CD34", "CAV1", "CLEC14A", "PLVAP"), reduction = "umap.scvi") #endothelial cell
FeaturePlot(validation, features = c("RUNX2", "ACAN"), reduction = "umap.scvi") # chondrocytes
FeaturePlot(validation, features = c("JCHAIN", "IGHG1", "IGLC2", "IGHG4"), reduction = "umap.scvi") #plasma cell
FeaturePlot(validation, features = c("CD79A", "MS4A1", "IGHM", "CD19", "BANK1"), reduction = "umap.scvi") #B cell
FeaturePlot(validation, features = c("MS4A2", "TPSB2", "TPSAB1", "CPA3"), reduction = "umap.scvi") #mast cells
FeaturePlot(validation, features = c("MYLPF", "MYL1", "DES", "MYBPC2"), reduction = "umap.scvi") #myoblast
FeaturePlot(validation, features = c("ACTA2", "RGS5", "PDGFRB", "NOTCH3", "MYL9", "TAGLN", "RGS5"), reduction = "umap.scvi") #pericytes
FeaturePlot(validation, features = c("CXCL12", "SFRP2", "MME"), reduction = "umap.scvi") #MSC (mesenchymal stem cells)
FeaturePlot(validation, features = c("SFTA3", "SFTPB", "SFTPC", "ROS1"), reduction = "umap.scvi") #lung cells
FeaturePlot(validation, features = c("COL1A1", "LUM", "DCN", "FBLN1", "ACTA2", "TAGLN", "COL3A1", "COL6A1", "FAP"), reduction = "umap.scvi") #Fibroblast
FeaturePlot(validation, features = c("CLEC4C", "LILRA4", "IL3RA"), reduction = "umap.scvi") #lung cells

validation <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/validation_post_scvi.rds")
Idents(validation) <- validation$scvi_clusters
validation <- RenameIdents(validation,
                           "0" = "Myeloid Cell",
                           "1" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "2" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "3" = "T/NK Cell", 
                           "4" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "5" = "Endothelial Cell",
                           "6" = "Endothelial Cell",
                           "7" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "8" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "9" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "10" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "11" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "12" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "13" = "Alveolar epithelial cells (Lung)", 
                           "14" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "15" = "Myeloid Cell",
                           "16" = "Endothelial Cell",
                           "17" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "18" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "19" = "Myeloid Cell",
                           "20" = "Osteoclasts",
                           "21" = "Plasma Cell",
                           "22" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "23" = "B Cell", 
                           "24" = "Alveolar epithelial cells (Lung)",
                           "25" = "osteoblastic/fibroblastic (OB/Fib) population",
                           "26" = "Myeloid Cell",
                           "27" = "Endothelial Cell",
                           "28" = "Myeloid Cell",
                           "29" = "Myoblasts",
                           "30" = "Mast Cells",
                           "31" = "Plasmacytoid dendritic cells (pDCs)")
validation <- JoinLayers(validation)
validation$cell_type <- Idents(validation)

DimPlot(validation, reduction = "umap.scvi", group.by = "cell_type", combine = FALSE, label.size = 2, label = TRUE)
DimPlot(validation, reduction = "umap.scvi", group.by = "sample", combine = FALSE, label.size = 2, label = TRUE)
DimPlot(validation, reduction = "umap.scvi", group.by = "SCEVAN.Malignancy", combine = FALSE, label.size = 2, label = TRUE)

validation_malignant <- subset(validation, subset = subset = SCEVAN.Malignancy == "tumor" & cell_type_annotated == "Osteoblast (OS)")

marker_genes <- c(
  # MSCs
  "CXCL12", "SFRP2", "MME",
  
  # Myeloid
  "CD74", "CD14", "FCGR3A", "LYZ", "CD163", "MS4A7",
  "S100A9", "APOE", "C1QA", "CD68",
  
  # Osteoblasts
  "CDH11", "RUNX2", "IBSP", "CLEC11A", "ALPL", "SATB2",
  
  # T/NK cells
  "IL7R", "CD3D", "NKG7", "CD3E", "GZMK", "GZMA",
  "TRAC", "CD3G", "TRBC1",
  
  # Endothelial cells
  "PECAM1", "VWF", "EGFL7", "CDH5", "CLDN5", "CD34",
  "CAV1", "CLEC14A", "PLVAP",
  
  # Alveolar epithelial cells (lung)
  "SFTA3", "SFTPB", "SFTPC", "ROS1", "ITGB6", "GPRC5A",
  
  # Osteoclasts
  "CTSK", "MMP9", "ACP5",
  
  # Plasma cells
  "JCHAIN", "IGHG1", "IGLC2", "IGHG4",
  
  # B cells
  "CD79A", "MS4A1", "IGHM", "CD19", "BANK1",
  
  # Myoblasts
  "MYLPF", "MYL1", "DES", "MYBPC2",
  
  # Mast cells
  "MS4A2", "TPSB2", "TPSAB1", "CPA3",
  
  # pDCs
  "CLEC4C", "LILRA4", "IL3RA"
)
##Marker gene dotplot
markers_dotplot <- DotPlot(validation, features = marker_genes, group.by = "cell_type") + RotatedAxis() + ggtitle("Top cluster markers")+
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1))
ggsave("marker_dotplot_val.pdf", markers_dotplot, width = 12, height = 8)


validation_subset <- subset(validation, subset = cell_type == 'osteoblastic/fibroblastic (OB/Fib) population' & SCEVAN.Malignancy == 'tumor')
saveRDS(validation_subset, file = "validation_subset.rds")
validation_subset <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/validation_subset")

#normal workflow
validation_subset <- NormalizeData(validation_subset)
validation_subset <- FindVariableFeatures(validation_subset)
validation_subset <- ScaleData(validation_subset)
validation_subset <- RunPCA(validation_subset)
ElbowPlot(validation_subset)
validation_subset <- FindNeighbors(validation_subset, dims = 1:15, reduction = "pca")
validation_subset <- FindClusters(validation_subset, resolution = 0.1, cluster.name = "malignant.clusters")
validation_subset <- RunUMAP(validation_subset, dims = 1:15, reduction = "pca", reduction.name = "umap.malignant")

DimPlot(validation_subset, reduction = "umap.malignant", group.by = "malignant.clusters", label = TRUE)
DimPlot(validation_subset, reduction = "umap.malignant", group.by = "sample") 
DimPlot(validation_subset, reduction = "umap.malignant", group.by = "dataset")
FeaturePlot(validation_subset, features = c("CD74"), reduction = "umap.malignant") #osteoblasts

validation_subset <- IntegrateLayers(
  object          = validation_subset,
  method          = scVIIntegration,
  new.reduction   = "integrated.malignant.scvi",
  conda_env       = "/Users/alfred/micromamba/envs/scvi_env",
  features        = Features(validation_subset),   # all genes, matching their unfiltered adata, null = variablefeatures
  layers          = "counts",
  groups.name     = "sample",
  groups          = validation_subset@meta.data[colnames(validation_subset), , drop = FALSE],      
  ndims.out       = 30,      # matches their n_latent=30
  n_hidden        = 256,     # matches
  n_layers        = 2,       # matches
  dropout_rate    = 0.1,     # matches
  max_epochs      = 200,     # matches
  gene_likelihood = "zinb",  # matches their unset default (scvi-tools defaults to zinb)
  seed.use        = 0,
  early_stopping = TRUE,
  early_stopping_patience = 30,
  early_stopping_monitor = "elbo_validation",
  torch.intraop.threads = 1L,
  torch.interop.threads = 1L,
  verbose = TRUE,
  model.save.dir = "/Users/alfred/Downloads"
)  #did the early stopping
saveRDS(validation_subset, file = "validation_subset.rds")
validation_subset <- FindNeighbors(validation_subset, dims = 1:30, reduction = "integrated.malignant.scvi")
validation_subset <- FindClusters(validation_subset, resolution = 0.2, cluster.name = "malignant.scvi.clusters")
validation_subset <- RunUMAP(validation_subset, dims = 1:30, reduction = "integrated.malignant.scvi", reduction.name = "umap.scvi.malignant")
DimPlot(validation_subset, reduction = "umap.scvi.malignant", group.by = "malignant.scvi.clusters", label = TRUE)
DimPlot(validation_subset, reduction = "umap.scvi.malignant", group.by = "sample") 
DimPlot(validation_subset, reduction = "umap.scvi.malignant", group.by = "orig.ident")
DimPlot(validation_subset, reduction = "umap.scvi.malignant", group.by = "Tissue_Type")

valmark1<-FindMarkers(validation_subset, ident.1 = "2", min.pct = 0.25)

FeaturePlot(validation_subset, features = c("SPP1"), reduction = "umap.scvi.malignant") #osteoblasts
FeaturePlot(validation_subset, features = c("IFI44L", "IFI6", "IFI44L", "ISG15"), reduction = "umap.scvi.malignant") #interferon
FeaturePlot(validation_subset, features = c("CENPF", "TOP2A", "KIF20A", "CDC20"), reduction = "umap.scvi.malignant") #proliferating
FeaturePlot(validation_subset, features = c("COL2A1", "S100A1"), reduction = "umap.scvi.malignant") #chondroblastic
