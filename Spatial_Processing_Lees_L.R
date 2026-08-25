library(Seurat)
library(UCell)
library(Voyager) #BiocManager::install("Voyager")
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(terra)
set.seed(123)

#Examples of loading in ST datasets and preprocessing with custom function, https://github.com/satijalab/seurat/issues/9072
load_spatial <- function(spatial_path, cell2loc_path){
  data.dir <- spatial_path
  seurat <- Load10X_Spatial(data.dir, filename="filtered_feature_bc_matrix.h5")
  seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
  seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")
  seurat <- subset(
    seurat, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 &
      nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)
  seurat <- SCTransform(seurat, assay = "Spatial")
  cell_ab <- read.csv(cell2loc_path)
  rownames <- cell_ab[,1]
  cell_ab <- cell_ab[,-1]
  rownames(cell_ab) <- rownames
  seurat<-AddMetaData(seurat, cell_ab)
}
MR15 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/MR15", "/Users/alfred/Downloads/STOsteosarcoma/MR15/abundance.csv")
MR10A <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/MR10A", "/Users/alfred/Downloads/STOsteosarcoma/MR10A/abundance.csv")
PT5 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/PT5", "/Users/alfred/Downloads/STOsteosarcoma/PT5/abundance.csv")
PT6 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/PT6", "/Users/alfred/Downloads/STOsteosarcoma/PT6/abundance.csv")
MR12 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/MR12", "/Users/alfred/Downloads/STOsteosarcoma/MR12/abundance.csv")
MR13 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/MR13", "/Users/alfred/Downloads/STOsteosarcoma/MR13/abundance.csv")
MR3 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/MR3", "/Users/alfred/Downloads/STOsteosarcoma/MR3/abundance.csv")
MR4 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/MR4", "/Users/alfred/Downloads/STOsteosarcoma/MR4/abundance.csv")
PT3 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/PT3", "/Users/alfred/Downloads/STOsteosarcoma/PT3/abundance.csv")
PT6 <- load_spatial("/Users/alfred/Downloads/Osteosarcoma/PT6", "/Users/alfred/Downloads/STOsteosarcoma/PT6/abundance.csv")

#Adding signature scores for hypoxia, EMT, and angiogenesis via UCell package
signatures <- list(
  Angiogenesis<-c(
    "APOH", "APP", "CCND2", "COL3A1", "COL5A2", "CXCL6", "FGFR1", "FSTL1",
    "ITGAV", "JAG1", "JAG2", "KCNJ8", "LPL", "LRPAP1", "LUM", "MSX1", "NRP1",
    "OLR1", "PDGFA", "PF4", "PGLYRP1", "POSTN", "PRG2", "PTK2", "S100A4",
    "SERPINA5", "SLCO2A1", "SPP1", "STC1", "THBD", "TIMP1", "TNFRSF21", "VAV2",
    "VCAN", "VEGFA", "VTN"
  ),
  Hypoxia<-c(
    "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3",
    "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2",
    "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12",
    "CASP6", "CAV1", "CCNG2", "NOCT", "CDKN1A", "CDKN1B", "CDKN1C", "CHST2",
    "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CCN2", "CXCR4", "ACKR3",
    "CCN1", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1", "EDN2",
    "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1A", "ERRFI1",
    "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", "GAA",
    "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3",
    "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2",
    "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3",
    "IL6", "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A",
    "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE1", "LDHA", "LDHC", "LOX",
    "LXN", "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK",
    "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NR3C1", "P4HA1",
    "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP",
    "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1",
    "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A",
    "PPP1R3C", "PRDX5", "PRKCA", "CAVIN3", "CAVIN1", "PYGM", "RBPJ", "RORA",
    "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1",
    "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1", "SLC2A3", "SLC2A5", "SLC37A4",
    "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3",
    "TGFBI", "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPBG", "TPD52",
    "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "CCN5", "WSB1", "XPNPEP1",
    "ZFP36", "ZNF292"
  ),
  EMT<-c(
    "ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1", "BDNF",
    "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG", "CD44", "CD59",
    "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", "COL1A1",
    "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL5A3",
    "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COMP", "COPA", "CRLF1", "CCN2",
    "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CCN1", "DAB2", "DCN", "DKK1",
    "DPYSL3", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", "EMP3",
    "ENO2", "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2", "FERMT2",
    "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1", "FSTL3", "FUCA1", "FZD8",
    "GADD45A", "GADD45B", "GAS1", "GEM", "GJA1", "GLIPR1", "COLGALT1", "GPC1",
    "GPX7", "GREM1", "HTRA1", "ID2", "IGFBP2", "IGFBP3", "IGFBP4", "IL15",
    "IL32", "IL6", "CXCL8", "INHBA", "ITGA2", "ITGA5", "ITGAV", "ITGB1",
    "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMC1", "LAMC2",
    "P3H1", "LGALS1", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC15", "LUM",
    "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5", "MGP", "MMP1",
    "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "MYL9", "MYLK", "NID2", "NNMT",
    "NOTCH2", "NT5E", "NTM", "OXTR", "PCOLCE", "PCOLCE2", "PDGFRB", "PDLIM4",
    "PFN2", "PLAUR", "PLOD1", "PLOD2", "PLOD3", "PMEPA1", "PMP22", "POSTN",
    "PPIB", "PRRX1", "PRSS2", "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB",
    "SAT1", "SCG2", "SDC1", "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1",
    "SFRP4", "SGCB", "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3", "SNAI2",
    "SNTB1", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TFPI2", "TGFB1", "TGFBI",
    "TGFBR3", "TGM2", "THBS1", "THBS2", "THY1", "TIMP1", "TIMP3", "TNC",
    "TNFAIP3", "TNFRSF11B", "TNFRSF12A", "TPM1", "TPM2", "TPM4", "VCAM1",
    "VCAN", "VEGFA", "VEGFC", "VIM", "WIPF1", "WNT5A"
  )
)

#examples of adding this signature to Seurat objects
MR15 <- AddModuleScore_UCell(MR15, features = signatures)
MR10A <- AddModuleScore_UCell(MR10A, features = signatures)
PT5 <- AddModuleScore_UCell(PT5, features = signatures)
PT6 <- AddModuleScore_UCell(PT6, features = signatures)
PT3 <- AddModuleScore_UCell(PT3, features = signatures)

#Examples of plotting features of interest (malignant osteoblast subcluster abundance, gene expression, pathway expression, etc.)
SpatialFeaturePlot(MR15, features=c("Endothelial.Cell", "Osteoblast_subcluster_2", "VEGFA", "FGFBP2", "S100A1", "Hypoxia_UCell"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(MR10A, features=c("Endothelial.Cell", "Osteoblast_subcluster_2", "VEGFA", "FGFBP2", "S100A1", "Hypoxia_UCell"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(PT5, features=c("Endothelial.Cell", "Osteoblast_subcluster_2", "VEGFA", "FGFBP2", "S100A1", "Hypoxia_UCell"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(PT6, features=c("Endothelial.Cell", "Osteoblast_subcluster_2", "VEGFA", "FGFBP2", "S100A1", "Hypoxia_UCell"), pt.size.factor = 3, alpha = 1)

SpatialFeaturePlot(MR15, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(MR10A, features=c("Osteoblast_subcluster_0", "SPP1", "EMT_UCell"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(PT3, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(PT5, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(PT6, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(MR12, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(MR13, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(MR3, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(MR4, features=c("Osteoblast_subcluster_0", "SPP1"), pt.size.factor = 3, alpha = 1)

SpatialFeaturePlot(PT3, features=c("Osteoblast_subcluster_3", "T.NK.Cell"), pt.size.factor = 3, alpha = 1)
SpatialFeaturePlot(PT6, features=c("Osteoblast_subcluster_3", "T.NK.Cell"), pt.size.factor = 3, alpha = 1)

##Lee's L statistic
##Custom function to convert Seurat to SpatialExperiment format: https://github.com/drighelli/SpatialExperiment/issues/115
seurat_to_spe <- function(seu, sample_id, img_id) {
  library(dplyr)
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(GetTissueCoordinates(seu)[, c("x", "y")])
  #spatialCoords <- as.matrix(
  #seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")]) #this is visium v1 but need different for v2
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = NULL #set to null for now!!
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
} #set imgData null due to spatraster issues
MR15_Voyager <- seurat_to_spe(MR15, "MR15", "slice1")

#Converting SpatialExperiment to SpatialFeatureExperiment
MR15_Voyager<-toSpatialFeatureExperiment( #https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-coercion.html
  MR15_Voyager,
  colGeometries = NULL,
  rowGeometries = NULL,
  annotGeometries = NULL,
  spatialCoordsNames = c("x", "y"),
  annotGeometryType = "POLYGON",
  spatialGraphs = NULL,
  spotDiameter = NA,
  unit = NULL
)

#Example using MR15, following vignette for Voyager:https://pachterlab.github.io/voyager/articles/bivariate.html
#create neighborhood graph
colGraph(MR15_Voyager, "visium") <- findVisiumGraph(MR15_Voyager, zero.policy = TRUE)

##Calculating Lee's L
#gene + gene
res <- calculateBivariate(
  MR15_Voyager,
  type = "lee.mc",
  feature1 = "TP53",
  feature2 = "MDM2",
  colGraphName = "visium",
  nsim = 999,
  zero.policy = TRUE
)
res #display values

#metadata + metadata
set.seed(123)
res <- calculateBivariate(
  x = MR15_Voyager$Osteoblast_subcluster_2,
  y = MR15_Voyager$Hypoxia_UCell,
  type = "lee.mc",
  listw = colGraph(MR15_Voyager, "visium"),
  nsim = 999,
  zero.policy = TRUE
)

#gene + metadata
ob <- MR15_Voyager$Osteoblast_subcluster_2
gene <- logcounts(MR15_Voyager)["VEGFA", ]
set.seed(123)
res <- calculateBivariate(
  x = ob,
  y = gene,
  type = "lee.mc",
  listw = colGraph(MR15_Voyager, "visium"),
  nsim = 999,
  zero.policy = TRUE
)
res
