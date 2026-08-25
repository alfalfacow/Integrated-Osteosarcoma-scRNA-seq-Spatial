## Step 0: Preparing for SeuratIntegrate
# Python / scVI initialization (should restart R session)
Sys.getenv("RETICULATE_PYTHON") #to test
library(reticulate)
py_discover_config() #to test
py_config()
torch <- import("torch")

# Prevent the R + PyTorch threading crash
torch$set_num_threads(1L)
torch$set_num_interop_threads(1L)

# Load scVI + settings
scvi <- import("scvi")
scvi$settings$num_threads <- 1L
scvi$settings$dl_num_workers <- 0L

# Forcing R to use correct C++ compiler (run once to edit Renviron)
mkdir -p ~/.R
echo "CXX14 = clang++ -std=gnu++14" >> ~/.R/Makevars
echo "CXX14STD = -std=gnu++14" >> ~/.R/Makevars
echo "CXX11STD = -std=gnu++14" >> ~/.R/Makevars
cat ~/.R/Makevars
#Then, install gfortran binary https://mac.r-project.org/tools/

#Loading packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(Seurat)
library(SeuratWrappers)
library(SeuratIntegrate)
library(SCEVAN)
library(DoubletFinder) #https://github.com/chris-mcginnis-ucsf/DoubletFinder
library(CellChat)
set.seed(123)

##Vignettes/Tutorials/Resources used:
#https://satijalab.org/seurat/articles/pbmc3k_tutorial
#https://satijalab.org/seurat/articles/seurat5_integration
#https://www.datanovia.com/learn/bioinformatics/single-cell/scrnaseq-quality-control
#http://datanovia.com/learn/bioinformatics/single-cell/scrnaseq-quality-control
#https://github.com/cbib/Seurat-Integrate
#https://biostatsquid.com/doubletfinder-tutorial/

## Step 1: Loading, QC/Filtering, and processing datasets
LoadQCFilterProcess <- function(path, projectID, sampleID, tissue_type){
  #loading data
  Seu_temp <- Read10X(data.dir = path)
  if(is.list(Seu_temp)){ #if multiple data types are present (ex: ATAC-seq)
    Seu_temp <- Seu_temp[["Gene Expression"]] #then only keep gene expression data
    print("Keeping only Gene Expression type")
  }
  Seu <- CreateSeuratObject(counts = Seu_temp, project = projectID, min.cells = 3, min.features = 200)
  Seu[["Tissue_Type"]] <- tissue_type
  Seu[["sample"]] <- sampleID
  
  #Filtering
  Seu[["percent.mt"]] <- PercentageFeatureSet(Seu, pattern = "^MT-")
  Seu <- subset(Seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20 & nCount_RNA < 30000)
  
  #Processing
  Seu <- NormalizeData(Seu)
  Seu <- FindVariableFeatures(Seu)
  Seu <- ScaleData(Seu)
  Seu <- RunPCA(Seu, nfeatures.print = 10)
  return(Seu)
}
Obj1 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE162454/OS_1", "GSE162454", "GSE162454_OS1", "primary")
Obj2 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE162454/OS_2", "GSE162454", "GSE162454_OS2", "primary")
Obj3 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE162454/OS_3", "GSE162454", "GSE162454_OS3", "primary")
Obj4 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE162454/OS_4", "GSE162454", "GSE162454_OS4", "primary")
Obj5 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE162454/OS_5", "GSE162454", "GSE162454_OS5", "primary")
Obj6 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE162454/OS_6", "GSE162454", "GSE162454_OS6", "primary")
Obj7 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/S0058/filtered_feature_bc_matrix", "GSE270231", "GSE270231_S0058", "lung_metastasis")
Obj8 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/S0059/filtered_feature_bc_matrix", "GSE270231", "GSE270231_S0059", "lung_metastasis")
Obj9 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/S0217/filtered_feature_bc_matrix", "GSE270231", "GSE270231_S0217", "lung_metastasis")
Obj10 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/S0218/filtered_feature_bc_matrix", "GSE270231", "GSE270231_S0218", "lung_metastasis")
Obj11 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/SC069/filtered_feature_bc_matrix", "GSE270231", "GSE270231_SC069", "lung_metastasis")
Obj12 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/SC072/filtered_feature_bc_matrix", "GSE270231", "GSE270231_SC072", "lung_metastasis")
Obj13 <- LoadQCFilterProcess("/Users/alfred/Downloads/Osteosarcoma/GSE270231/SC073/filtered_feature_bc_matrix", "GSE270231", "GSE270231_SC073", "lung_metastasis")

## Step 2: DoubletFinder
# Function from: https://biostatsquid.com/doubletfinder-tutorial/, slightly modified
# run_doubletfinder_custom runs Doublet_Finder()
run_doubletfinder_custom <- function(sample, multiplet_rate = NULL){
 
  # DoubletFinder requires 5 inputs: Seurat Object, PCs, pk, pN (default 25%), nExp.
  
  ## This step Calculates multiplet rate based on published proportions by 10x
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(sample@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  ## pK identification (no ground-truth) 
  # introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  #prints number of doublets vs singles predicted
  print(table(sample$doublet_finder))
  
  return(sample)
} #Returns full seurat object with metadata column designating singlet vs doublet
run_doubletfinder_custom_subset <- function(sample, multiplet_rate = NULL){
  
  # DoubletFinder requires 5 inputs: Seurat Object, PCs, pk, pN (default 25%), nExp.
  
  ## This step Calculates multiplet rate based on published proportions by 10x
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(sample@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  ## pK identification (no ground-truth) 
  # introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  #prints number of doublets vs singles predicted
  print(table(sample$doublet_finder))
  sample <- subset(sample, subset = doublet_finder == 'Singlet')
  return(sample)
} #Returns subset seurat object with only singlets
Obj1 <- run_doubletfinder_custom_subset(Obj1) #Doublets: 217, Singlets: 6389
Obj2 <- run_doubletfinder_custom_subset(Obj2) #Doublets: 227, Singlets: 6123
Obj3 <- run_doubletfinder_custom_subset(Obj3) #Doublets: 253, Singlets: 7634
Obj4 <- run_doubletfinder_custom_subset(Obj4) #Doublets: 58, Singlets: 3281
Obj5 <- run_doubletfinder_custom_subset(Obj5) #Doublets: 322, Singlets: 8605
Obj6 <- run_doubletfinder_custom_subset(Obj6) #Doublets: 333, Singlets: 7540
Obj7 <- run_doubletfinder_custom_subset(Obj7) #Doublets: 39, Singlets: 2918
Obj8 <- run_doubletfinder_custom_subset(Obj8) #Doublets: 220, Singlets: 6289
Obj9 <- run_doubletfinder_custom_subset(Obj9) #Doublets: 98, Singelts: 4757
Obj10 <- run_doubletfinder_custom_subset(Obj10) #Doublets: 144, Singlets: 4974
Obj11 <- run_doubletfinder_custom_subset(Obj11) #Doublets: 317, Singlets: 8171
Obj12 <- run_doubletfinder_custom_subset(Obj12) #Doublets: 3, Singlets: 992
Obj13 <- run_doubletfinder_custom_subset(Obj13) #Doublets: 26, Singlets: 2594

#After running, remove doublet score metadata (no longer needed)
Obj1$pANN_0.25_0.03_217 <- NULL
Obj2$pANN_0.25_0.28_227 <- NULL
Obj3$pANN_0.25_0.26_253 <- NULL
Obj4$pANN_0.25_0.005_58 <- NULL
Obj5$pANN_0.25_0.08_322 <- NULL
Obj6$pANN_0.25_0.005_333 <- NULL
Obj7$pANN_0.25_0.005_39 <- NULL
Obj8$pANN_0.25_0.3_220 <- NULL
Obj9$pANN_0.25_0.18_98 <- NULL
Obj10$pANN_0.25_0.005_144 <- NULL
Obj11$pANN_0.25_0.16_317 <- NULL
Obj12$pANN_0.25_0.3_3 <- NULL
Obj13$pANN_0.25_0.06_26 <- NULL

## Step 3: SCEVAN to distinguish malignant vs normal cells
runSCEVAN <- function(sample, sampleID){
  start_time <- Sys.time()
  count_mtx <- GetAssayData(
    sample,
    assay = "RNA",
    layer = "counts"
  )
  SCEVAN_results <- pipelineCNA(
    count_mtx = count_mtx,
    sample = sampleID,
    par_cores = 1,
    FIXED_NORMAL_CELLS = FALSE,
    SUBCLONES = FALSE,
    beta_vega = 0.5,
    plotTree = FALSE,
    SCEVANsignatures = TRUE,
    organism = "human"
  )
  print(head(SCEVAN_results)) #prints part of results to cross reference
  sample[["SCEVAN.Malignancy"]] <- SCEVAN_results[colnames(sample),"class"]
  end_time <- Sys.time()
  
  print(end_time - start_time)
  cat("SCEVAN took", end_time - start_time, "\n")
  
  return(sample)
}
Obj1 <- runSCEVAN(Obj1, "GSE162454_OS1")
Obj2 <- runSCEVAN(Obj2, "GSE162454_OS2")
Obj3 <- runSCEVAN(Obj3, "GSE162454_OS3")
Obj4 <- runSCEVAN(Obj4, "GSE162454_OS4")
Obj5 <- runSCEVAN(Obj5, "GSE162454_OS5")
Obj6 <- runSCEVAN(Obj6, "GSE162454_OS6")
Obj7 <- runSCEVAN(Obj7, "GSE270231_S0058")
Obj8 <- runSCEVAN(Obj8, "GSE270231_S0059")
Obj9 <- runSCEVAN(Obj9, "GSE270231_S0217")
Obj10 <- runSCEVAN(Obj10, "GSE270231_S0218")
Obj11 <- runSCEVAN(Obj11, "GSE270231_S0069")
Obj12 <- runSCEVAN(Obj12, "GSE270231_S0072")
Obj13 <- runSCEVAN(Obj13, "GSE270231_S0073")

#saving RDS files
saveRDS(Obj1, file = "GSE162454_OS1.rds")
saveRDS(Obj2, file = "GSE162454_OS2.rds")
saveRDS(Obj3, file = "GSE162454_OS3.rds")
saveRDS(Obj4, file = "GSE162454_OS4.rds")
saveRDS(Obj5, file = "GSE162454_OS5.rds")
saveRDS(Obj6, file = "GSE162454_OS6.rds")
saveRDS(Obj7, file = "GSE270231_S0058.rds")
saveRDS(Obj8, file = "GSE270231_S0059.rds")
saveRDS(Obj9, file = "GSE270231_S0217.rds")
saveRDS(Obj10, file = "GSE270231_S0218.rds")
saveRDS(Obj11, file = "GSE270231_SC069.rds")
saveRDS(Obj12, file = "GSE270231_SC072.rds")
saveRDS(Obj13, file = "GSE270231_SC073.rds")

Obj1 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE162454_OS1.rds")
Obj2 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE162454_OS2.rds")
Obj3 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE162454_OS3.rds")
Obj4 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE162454_OS4.rds")
Obj5 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE162454_OS5.rds")
Obj6 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE162454_OS6.rds")
Obj7 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_S0058.rds")
Obj8 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_S0059.rds")
Obj9 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_S0217.rds")
Obj10 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_S0218.rds")
Obj11 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_SC069.rds")
Obj12 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_SC072.rds")
Obj13 <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/GSE270231_SC073.rds")


#Step 4: Loading RDS objects from Elizabeth + merging
Obj14 <- readRDS("/Users/alfred/Downloads/Archive/BC2_with_scevan.rds")
Obj15 <- readRDS("/Users/alfred/Downloads/Archive/BC3_with_scevan.rds")
Obj16 <- readRDS("/Users/alfred/Downloads/Archive/BC5_with_scevan.rds")
Obj17 <- readRDS("/Users/alfred/Downloads/Archive/BC6_with_scevan.rds")
Obj18 <- readRDS("/Users/alfred/Downloads/Archive/BC16_with_scevan.rds")
Obj19 <- readRDS("/Users/alfred/Downloads/Archive/BC21_with_scevan.rds")
Obj20 <- readRDS("/Users/alfred/Downloads/Archive/BC22_with_scevan.rds")
Obj21 <- readRDS("/Users/alfred/Downloads/Archive/BC10_with_scevan.rds")
Obj22 <- readRDS("/Users/alfred/Downloads/Archive/BC17_with_scevan.rds")

Integrated <- merge(Obj1, c(Obj2, Obj3, Obj4, Obj5, Obj6, Obj7, Obj8, Obj9, Obj10, Obj11, Obj12, Obj13, Obj14, Obj15, Obj16, Obj17, Obj18, Obj19, Obj20, Obj21, Obj22), merge.data = TRUE, project = "Integrated_Data")
Integrated$pANN_0.25_0.2_49 <- NULL
saveRDS(Integrated, file = "Integrated.rds")
Integrated <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/Integrated.rds")
rm(Obj1, Obj2, Obj3, Obj4, Obj5, Obj6, Obj7, Obj8, Obj9, Obj10, Obj11, Obj12, Obj13, Obj14, Obj15, Obj16, Obj17, Obj18, Obj19, Obj20, Obj21, Obj22)

#Step 5: Standard Workflow on merged object
#verify QC and filtering
VlnPlot(Integrated, features = c("percent.mt"), group.by = "sample")
#"nCount_RNA, "nFeature_RNA", "percent.mt"

Integrated <- NormalizeData(Integrated)
Integrated <- FindVariableFeatures(Integrated)
Integrated <- ScaleData(Integrated)
Integrated <- RunPCA(Integrated)

#visualization of unintegrated clusters
Integrated <- FindNeighbors(Integrated, dims = 1:20, reduction = "pca")
Integrated <- FindClusters(Integrated, resolution = 0.5, cluster.name = "unintegrated_clusters")
Integrated <- RunUMAP(Integrated, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(Integrated, reduction = "umap.unintegrated", group.by = "unintegrated_clusters") 
DimPlot(Integrated, reduction = "umap.unintegrated", group.by = "orig.ident")
DimPlot(Integrated, reduction = "umap.unintegrated", group.by = "sample") #not very well integrated...
DimPlot(Integrated, reduction = "umap.unintegrated", group.by = "Tissue_Type")
DimPlot(Integrated, reduction = "umap.unintegrated", group.by = "SCEVAN.Malignancy")


#Step 6: scVI Integration
Integrated <- IntegrateLayers(
  object          = Integrated,
  method          = scVIIntegration,
  new.reduction   = "integrated.scvi",
  conda_env       = "/Users/alfred/micromamba/envs/scvi_env",
  features        = Features(Integrated),   # all genes, matching their unfiltered adata, null = variablefeatures
  layers          = "counts",
  groups.name     = "sample",
  groups          = Integrated[[]],      
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

#Get and save embeddings
scvi_embeddings <- Embeddings(Integrated, reduction = "integrated.scvi") #pulls embeddings out of object
scvi_df <- as.data.frame(scvi_embeddings)
scvi_df <- cbind(cell = rownames(scvi_df), scvi_df) #makes rownames (cell ID) into its own column
write.csv(scvi_df, file = "scvi_latent.csv", row.names = FALSE)
rm(scvi_embeddings)
rm(scvi_df)

Integrated <- FindNeighbors(Integrated, reduction = "integrated.scvi", dims = 1:30)
Integrated <- FindClusters(Integrated, resolution = 0.8, cluster.name = "scvi_clusters") #harrison's script uses 0.8
Integrated <- RunUMAP(Integrated, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
DimPlot(Integrated, reduction = "umap.scvi", group.by = "orig.ident", combine = FALSE, label.size = 2)
DimPlot(Integrated, reduction = "umap.scvi", group.by = "sample", combine = FALSE, label.size = 2)
DimPlot(Integrated, reduction = "umap.scvi", group.by = "Tissue_Type", combine = FALSE, label.size = 2)
DimPlot(Integrated, reduction = "umap.scvi", group.by = "scvi_clusters", combine = FALSE, label.size = 2, label = TRUE)
DimPlot(Integrated, reduction = "umap.scvi", group.by = "SCEVAN.Malignancy", combine = FALSE, label.size = 2)
?DimPlot

#Step 7: cell type annotation
FeaturePlot(Integrated, features = c("COL1A1", "CDH11", "RUNX2", "IBSP", "CLEC11A", "ALPL", "SATB2"), reduction = "umap.scvi") #osteoblasts
FeaturePlot(Integrated, features = c("CTSK", "MMP9", "ACP5"), reduction = "umap.scvi") #osteoclasts
FeaturePlot(Integrated, features = c("IL7R", "CD3D", "NKG7", "CD3E", "GZMK", "GZMA", "TRAC", "CD3G", "TRBC1"), reduction = "umap.scvi") #T/NK Cells
FeaturePlot(Integrated, features = c("CD74", "CD14", "FCGR3A", "LYZ", "CD163", "MS4A7", "S100A9", "APOE", "C1QA", "CD68"), reduction = "umap.scvi") #Myeloid Cells
FeaturePlot(Integrated, features = c("COL1A1", "LUM", "DCN", "FBLN1", "ACTA2", "TAGLN", "COL3A1", "COL6A1", "FAP"), reduction = "umap.scvi") #Fibroblast
FeaturePlot(Integrated, features = c("PECAM1", "VWF", "EGFL7", "CDH5", "CLDN5", "CD34", "CAV1", "CLEC14A", "PLVAP"), reduction = "umap.scvi") #endothelial cell
FeaturePlot(Integrated, features = c("RUNX2", "ACAN"), reduction = "umap.scvi") # chondrocytes
FeaturePlot(Integrated, features = c("JCHAIN", "IGHG1", "IGLC2", "IGHG4"), reduction = "umap.scvi") #plasma cell
FeaturePlot(Integrated, features = c("CD79A", "MS4A1", "IGHM", "CD19", "BANK1"), reduction = "umap.scvi") #B cell
FeaturePlot(Integrated, features = c("MS4A2", "TPSB2", "TPSAB1", "CPA3"), reduction = "umap.scvi") #mast cells
FeaturePlot(Integrated, features = c("MYLPF", "MYL1", "DES", "MYBPC2"), reduction = "umap.scvi") #myoblast
FeaturePlot(Integrated, features = c("ACTA2", "RGS5", "PDGFRB", "NOTCH3", "MYL9", "TAGLN", "RGS5"), reduction = "umap.scvi") #pericytes
FeaturePlot(Integrated, features = c("CXCL12", "SFRP2", "MME"), reduction = "umap.scvi") #MSC (mesenchymal stem cells)
FeaturePlot(Integrated, features = c(), reduction = "umap.scvi")



#Renaming cluster numbers with cell type
Idents(Integrated) <- "scvi_clusters" #set column values of active.ident to be changed to scvi clusters
Integrated <- RenameIdents(Integrated,
                                    "0" = "Myeloid Cell",
                                    "1" = "Osteoblast (OS)",
                                    "2" = "Osteoblast (OS)",
                                    "3" = "T/NK Cell", 
                                    "4" = "Myeloid Cell",
                                    "5" = "5",
                                    "6" = "Osteoblast (OS)",
                                    "7" = "Osteoclast",
                                    "8" = "Myeloid Cell",
                                    "9" = "Endothelial Cell",
                                    "10" = "Osteoblast (OS)",
                                    "11" = "Osteoblast (OS)",
                                    "12" = "T/NK Cell",
                                    "13" = "Osteoclast", 
                                    "14" = "14",
                                    "15" = "Myeloid Cell",
                                    "16" = "16",
                                    "17" = "Mesenchymal Stem Cell",
                                    "18" = "T/NK Cell",
                                    "19" = "Osteoblast (OS)",
                                    "20" = "Myeloid Cell",
                                    "21" = "21",
                                    "22" = "22",
                                    "23" = "T/NK Cell", 
                                    "24" = "Osteoclast",
                                    "25" = "Plasma Cell",
                                    "26" = "Mast Cell",
                                    "27" = "B Cell",
                                    "28" = "Osteoclast",
                                    "29" = "29",
                                    "30" = "Osteoblast (OS)",
                                    "31" = "31",
                                    "32" = "Myeloid Cell",
                                    "33" = "33", 
                                    "34" = "Myeloid Cell",
                                    "35" = "Myoblast",
                                    "36" = "36",
                                    "37" = "Osteoblast (OS)",
                                    "38" = "38",
                                    "39" = "Osteoblast (OS)")
Integrated[["cell_type_annotated"]] <- Idents(Integrated) #create new metadata column
DimPlot(Integrated, reduction = "umap.scvi", group.by = "cell_type_annotated", combine = FALSE, label.size = 2, label = TRUE)

#resolving non annotated clusters
Integrated <- JoinLayers(Integrated) #needed to run find markers
markers5 <- FindMarkers(Integrated, ident.1 = "5", min.pct = 0.25)
markers14 <- FindMarkers(Integrated, ident.1 = "14", min.pct = 0.25)
markers16 <- FindMarkers(Integrated, ident.1 = "16", min.pct = 0.25)
markers21 <- FindMarkers(Integrated, ident.1 = "21", min.pct = 0.25)
markers22 <- FindMarkers(Integrated, ident.1 = "22", min.pct = 0.25) #alveolar epithelial cells (from lung metastasis)
markers29 <- FindMarkers(Integrated, ident.1 = "29", min.pct = 0.25) 
markers31 <- FindMarkers(Integrated, ident.1 = "31", min.pct = 0.25) #plasmacytoid dendritic cells
markers33 <- FindMarkers(Integrated, ident.1 = "33", min.pct = 0.25)
markers36 <- FindMarkers(Integrated, ident.1 = "36", min.pct = 0.25) #erythrocyte/blood cell
markers38 <- FindMarkers(Integrated, ident.1 = "38", min.pct = 0.25)
Integrated[["RNA"]] <- split(Integrated[["RNA"]], f = Integrated$sample) #resplit the layers after markers

FeaturePlot(Integrated, features = c("SFTPB", "ITGB6", "SFTA3"), reduction = "umap.scvi") #22
FeaturePlot(Integrated, features = c("CLEC4C", "LRRC26", "LILRA4", "TCL1A"), reduction = "umap.scvi") #31
FeaturePlot(Integrated, features = c("HBD", "HBB", "HBA1", "HBA2"), reduction = "umap.scvi") #36
FeaturePlot(Integrated, features = c("PTPRC"), reduction = "umap.scvi") #36

Integrated <- RenameIdents(Integrated,
                           "5" = "Non-malignant mesenchymal/stromal population", 
                           "14" = "Non-malignant mesenchymal/stromal population",
                           "16" = "Non-malignant mesenchymal/stromal population",
                           "21" = "Non-malignant mesenchymal/stromal population",
                           "22" = "Alveolar epithelial cells (Lung)",
                           "29" = "Non-malignant mesenchymal/stromal population",
                           "31" = "Plasmacytoid dendritic cells (pDCs)",
                           "33" = "Non-malignant mesenchymal/stromal population",
                           "36" = "Erythrocytes",
                           "38" = "Non-malignant mesenchymal/stromal population")
Integrated[["cell_type_annotated"]] <- Idents(Integrated) #update cell_type
cell_types <- as.data.frame(Integrated$cell_type_annotated)
OS_subtypes <- as.data.frame(malignant$malignant.scvi.clusters)

#Step 7: Compile cell-type markers and generate heatmap + dot plot
marker_genes <- c("SFTA3", "SFTPB", "SFTPC", "ROS1", "ITGB6", "GPRC5A", #Alveolar epithelial cells (Lung)
                  "CLEC4C", "LILRA4", "IL3RA", #pDCs
                  "HBD", "HBB", "HBA1", "HBA2", "AHSP", "CA1", #erythrocytes
                  "CD74", "CD14", "FCGR3A", "LYZ", "CD163", "MS4A7", "S100A9", "APOE", "C1QA", "CD68", #Myeloid
                  "COL1A1", "CDH11", "RUNX2", "IBSP", "CLEC11A", "ALPL", "SATB2", #osteoblasts
                  "IL7R", "CD3D", "NKG7", "CD3E", "GZMK", "GZMA", "TRAC", "CD3G", "TRBC1", #T/NK
                  "CTSK", "MMP9", "ACP5", #osteoclasts
                  "PECAM1", "VWF", "EGFL7", "CDH5", "CLDN5", "CD34", "CAV1", "CLEC14A", "PLVAP", #endothelial
                  "CXCL12", "SFRP2", "MME", #MSCs
                  "JCHAIN", "IGHG1", "IGLC2", "IGHG4", #plasma
                  "MS4A2", "TPSB2", "TPSAB1", "CPA3", #mast cells
                  "CD79A", "MS4A1", "IGHM", "CD19", "BANK1", #B cells
                  "MYLPF", "MYL1", "DES", "MYBPC2") #myoblasts
  
#scrapped, won't work due to some populations being too small
markers_heatmap <- DoHeatmap(Integrated, features = marker_genes, group.by = "cell_type_annotated", slot = "data", size = 3) +
  ggtitle("Cell-type gene markers")+theme(axis.text.y = element_text(size = 3))
ggsave("marker_heatmap.pdf", markers_heatmap, width = 12, height = 10)

markers_dotplot <- DotPlot(Integrated, features = marker_genes, group.by = "cell_type_annotated") + RotatedAxis() + ggtitle("Top cluster markers")+
  theme(axis.text.x = element_text(size = 4, angle = 90, hjust = 1))
ggsave("marker_dotplot.pdf", markers_dotplot, width = 12, height = 8)

stromal_mark <- FindMarkers(Integrated, ident.1 = "Non-malignant mesenchymal/stromal population", min.pct = 0.25)
epithelial_mark <- FindMarkers(Integrated, ident.1 = "Alveolar epithelial cells (Lung)", min.pct = 0.25)
pDC_mark <- FindMarkers(Integrated, ident.1 = "Plasmacytoid dendritic cells (pDCs)", min.pct = 0.25)
blood_mark <- FindMarkers(Integrated, ident.1 = "Erythrocytes", min.pct = 0.25)
myeloid_mark <- FindMarkers(Integrated, ident.1 = "Myeloid Cell", min.pct = 0.25)
OS_mark <- FindMarkers(Integrated, ident.1 = "Osteoblast (OS)", min.pct = 0.25)
TNK_mark<- FindMarkers(Integrated, ident.1 = "T/NK Cell", min.pct = 0.25)
Osteoclast_mark <- FindMarkers(Integrated, ident.1 = "Osteoclast", min.pct = 0.25)
endo_mark <- FindMarkers(Integrated, ident.1 = "Endothelial Cell", min.pct = 0.25)
MSC_mark <- FindMarkers(Integrated, ident.1 = "Mesenchymal Stem Cell", min.pct = 0.25)
plasma_mark <- FindMarkers(Integrated, ident.1 = "Plasma Cell", min.pct = 0.25)
mast_mark <- FindMarkers(Integrated, ident.1 = "Mast Cell", min.pct = 0.25)
B_mark <- FindMarkers(Integrated, ident.1 = "B Cell", min.pct = 0.25)
myoblast_mark <- FindMarkers(Integrated, ident.1 = "Myoblast", min.pct = 0.25)

markers <- FindAllMarkers(
  Integrated,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  group.by = "cell_type_annotated"
)
write.csv(markers, "cell_type_markers.csv", row.names = FALSE)

#Step 8: Subclustering malignant osteoblast population
malignant <- subset(Integrated, subset = SCEVAN.Malignancy == "tumor" & cell_type_annotated == "Osteoblast (OS)")
saveRDS(malignant, file = "malignant.rds")
malignant <- readRDS("/Users/alfred/Desktop/Lab/Osteosarcoma/malignant.rds")

#normal workflow
malignant <- NormalizeData(malignant)
malignant <- FindVariableFeatures(malignant)
malignant <- ScaleData(malignant)
malignant <- RunPCA(malignant)
ElbowPlot(malignant)
malignant <- FindNeighbors(malignant, dims = 1:15, reduction = "pca")
malignant <- FindClusters(malignant, resolution = 0.1, cluster.name = "malignant.clusters")
malignant <- RunUMAP(malignant, dims = 1:15, reduction = "pca", reduction.name = "umap.malignant")

table(malignant$cell_type_annotated)

DimPlot(malignant, reduction = "umap.malignant", group.by = "malignant.clusters", label = TRUE)
DimPlot(malignant, reduction = "umap.malignant", group.by = "sample") 
DimPlot(malignant, reduction = "umap.malignant", group.by = "dataset")
FeaturePlot(malignant, features = c("CD74"), reduction = "umap.malignant") #osteoblasts

#significant batch effects, will have to rerun scVI
malignant <- IntegrateLayers(
  object          = malignant,
  method          = scVIIntegration,
  new.reduction   = "integrated.malignant.scvi",
  conda_env       = "/Users/alfred/micromamba/envs/scvi_env",
  features        = Features(Integrated),   # all genes, matching their unfiltered adata, null = variablefeatures
  layers          = "counts",
  groups.name     = "sample",
  groups          = Integrated@meta.data[colnames(malignant), , drop = FALSE],      
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
)
malignant <- FindNeighbors(malignant, dims = 1:30, reduction = "integrated.malignant.scvi")
malignant <- FindClusters(malignant, resolution = 0.5, cluster.name = "malignant.scvi.clusters")
malignant <- RunUMAP(malignant, dims = 1:30, reduction = "integrated.malignant.scvi", reduction.name = "umap.scvi.malignant")
DimPlot(malignant, reduction = "umap.scvi.malignant", group.by = "malignant.scvi.clusters", label = TRUE)
DimPlot(malignant, reduction = "umap.scvi.malignant", group.by = "sample") 
DimPlot(malignant, reduction = "umap.scvi.malignant", group.by = "orig.ident")
DimPlot(malignant, reduction = "umap.scvi.malignant", group.by = "Tissue_Type")

FeaturePlot(malignant, features = c("CDK3", "MYC", "UHRF2", "STC2", "COL5A2", "MMD", "EHMT2"), reduction = "umap.scvi.malignant")

#checking for biological differences
markers_malignant <- FindAllMarkers(
  malignant,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  group.by = "malignant.scvi.clusters"
)
markers_malignant <- markers_malignant[
  markers_malignant$p_val_adj < 0.05,
] #filter by p value
write.csv(markers_malignant, "malignant_subcluster_markers.csv")

saveRDS(malignant, file = "malignant.rds")
FeaturePlot(malignant, features = c("RELA"), reduction = "umap.scvi.malignant")

#Plotting dotplot for markers
marker_genes <- c("SPP1", "IBSP", "SERPINF1", #0
                  "CENPF", "TOP2A", "BIRC5", "PCLAF", "TACC3", #1
                  "S100A1", "VEGFA", "S100P", "SNORC", "COL2A1", "COL9A1", "SOX9", "ACAN", #2
                  "OASL", "RSAD2", "IFIT2", "IFIT1", "IFIT3", #3
                  "CD74", "C1QA", "C1QB", "C1QC", #4
                  "LUZP4", "C5orf58", "LINC02387", "LINC02087", "MAGEC2", "CCNE12" #5
                  )
DotPlot(malignant, features = marker_genes) + RotatedAxis() + ggtitle("Selected cluster markers")+
  theme(axis.text.x = element_text(size = 3, angle = 90, hjust = 1))
DoHeatmap(malignant, features = marker_genes, size = 3) +
  ggtitle("Top cluster markers")+theme(axis.text.y = element_text(size = 3))

# D. Plot
library(RColorBrewer)
n_patients <- length(unique(df_pct$patient))
pal <- colorRampPalette(brewer.pal(12, "Paired"))(n_patients)

d<-ggplot(df_pct, aes(x = cluster, y = pct, fill = patient)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Percentage", fill = "Patient") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8)
  )

#Step 9: Preparing scRNA-seq data as reference .h5ad for Cell2Location
cell_types <- Integrated$cell_type_annotated
OS_subtypes <- malignant$malignant.scvi.clusters
OS_subtypes$cellID <- rownames(OS_subtypes)
cell_types$cellID <- rownames(cell_types)
cell_types_merged <- left_join(cell_types, OS_subtypes, by = "cellID")
colnames(cell_types_merged) <- c("Cell_type", "CellID", "Osteoblast_subcluster")
cell_types_merged$Cell_type <- as.character(cell_types_merged$Cell_type)
cell_types_merged$Osteoblast_subcluster <- as.character(cell_types_merged$Osteoblast_subcluster)
for(i in 1:162802){
  if(cell_types_merged$Cell_type[i] == "Osteoblast (OS)"){
    if(is.na(cell_types_merged$Osteoblast_subcluster[i])){
      cell_types_merged$Cell_type[i] <- "Osteoblast (non malignant)"
    }
    else{
      cell_types_merged$Cell_type[i] <- paste("Osteoblast_subcluster", cell_types_merged$Osteoblast_subcluster[i], sep = "_")
    }
  }
}

cell_types_merged$Cell_type <- as.factor(cell_types_merged$Cell_type)
Integrated[["cell2loc_labels"]] <- cell_types_merged$Cell_type

library(anndata)
ad <- reticulate::import("anndata")
counts_t <- t(GetAssayData(Integrated, layer = "counts"))

adata <- ad$AnnData(
  X = counts_t,
  obs = Integrated@meta.data,
  var = data.frame(row.names = colnames(counts_t))
)

write_h5ad(adata, "osteosarcoma_reference.h5ad")