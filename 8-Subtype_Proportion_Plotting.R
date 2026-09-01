#Plotting cell type proportions from Cell2Location output
spatial_celltype_proportions_df <- function(path, patient){
  abundance <- read.csv(path)
  abundance <- abundance[,-1]
  colnames(abundance) <- c("Alveolar epithelial cells (Lung)", "B Cell", "Endothelial Cell",
                           "Erythrocytes", "Mast Cell", "Mesenchymal Stem Cell", "Myeloid Cell", 
                           "Myoblast", "Non-malignant mesenchymal/stromal population", "Osteoblast (non malignant)",
                           "Osteoblast_subcluster_0", "Osteoblast_subcluster_1", "Osteoblast_subcluster_2",
                           "Osteoblast_subcluster_3", "Osteoblast_subcluster_4", "Osteoblast_subcluster_5",
                           "Osteoclast", "Plasma Cell", "Plasmacytoid dendritic cells (pDCs)", "T/NK Cell")
  Sum <- 0
  for(i in 1:nrow(abundance)){
    for(j in 1:ncol(abundance)){
      Sum <- Sum + abundance[i,j]
    }
  }
  proportions <- data.frame(
    'Patient' = patient,
    'Alveolar epithelial cells (Lung)' = sum(abundance[,1])/Sum,
    'B Cell' = sum(abundance[,2])/Sum,
    'Endothelial Cell' = sum(abundance[,3])/Sum,
    'Erythrocytes' = sum(abundance[,4])/Sum,
    'Mast Cell' = sum(abundance[,5])/Sum,
    'Mesenchymal Stem Cell' = sum(abundance[,6])/Sum,
    'Myeloid Cell' = sum(abundance[,7])/Sum,
    'Myoblast' = sum(abundance[,8])/Sum,
    'Non-malignant mesenchymal/stromal population' = sum(abundance[,9])/Sum,
    'Osteoblast (non malignant)' = sum(abundance[,10])/Sum,
    'Osteoblast_subcluster_0' = sum(abundance[,11])/Sum,
    'Osteoblast_subcluster_1' = sum(abundance[,12])/Sum,
    'Osteoblast_subcluster_2' = sum(abundance[,13])/Sum,
    'Osteoblast_subcluster_3' = sum(abundance[,14])/Sum,
    'Osteoblast_subcluster_4' = sum(abundance[,15])/Sum,
    'Osteoblast_subcluster_5' = sum(abundance[,16])/Sum,
    'Osteoclast' = sum(abundance[,17])/Sum,
    'Plasma Cell' = sum(abundance[,18])/Sum,
    'Plasmacytoid dendritic cells (pDCs)' = sum(abundance[,19])/Sum,
    'T/NK Cell' = sum(abundance[,20])/Sum
  )
  return(proportions) 
}
spatial_celltype_proportions_OS_df <- function(path, patient){
  abundance <- read.csv(path)
  abundance <- abundance[,12:17]
  Sum <- 0
  for(i in 1:nrow(abundance)){
    for(j in 1:ncol(abundance)){
      Sum <- Sum + abundance[i,j]
    }
  }
  proportions <- data.frame(
    'Patient' = patient,
    'Osteoblast_subcluster_0' = sum(abundance[,1])/Sum,
    'Osteoblast_subcluster_1' = sum(abundance[,2])/Sum,
    'Osteoblast_subcluster_2' = sum(abundance[,3])/Sum,
    'Osteoblast_subcluster_3' = sum(abundance[,4])/Sum,
    'Osteoblast_subcluster_4' = sum(abundance[,5])/Sum,
    'Osteoblast_subcluster_5' = sum(abundance[,6])/Sum
  )
  return(proportions)
}
#samples: PT1-9 (primary), MR1-4 (lung, relapse1), MR 10-13, MR15, MR18

#primary samples data (n=9)
PT1 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT1/abundance.csv", "PT1")
PT2 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT2/abundance.csv", "PT2")
PT3 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT3/abundance.csv", "PT3")
PT4 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT4/abundance.csv", "PT4")
PT5 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT5/abundance.csv", "PT5")
PT6 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT6/abundance.csv", "PT6")
PT7 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT7/abundance.csv", "PT7")
PT8 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT8/abundance.csv", "PT8")
PT9 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/PT9/abundance.csv", "PT9")

#metastasis
MR1 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR1/abundance.csv", "MR1")
MR2 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR2/abundance.csv", "MR2")
MR3 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR3/abundance.csv", "MR3")
MR4 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR4/abundance.csv", "MR4")
MR10A <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR10A/abundance.csv", "MR10A")
MR10B <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR10B/abundance.csv", "MR10B")
MR11A <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR11A/abundance.csv", "MR11A")
MR11B <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR11B/abundance.csv", "MR11B")
MR12 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR12/abundance.csv", "MR12")
MR13 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR13/abundance.csv", "MR13")
MR15 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR15/abundance.csv", "MR15")
#MR18 <- spatial_celltype_proportions_OS_df("/Users/alfred/Downloads/STOsteosarcoma/MR18/abundance.csv", "MR18")

#merging into single dataframe
proportions <- bind_rows(PT1, PT2, PT3, PT4, PT5, PT6, PT7, PT8, PT9,
                         MR1, MR2, MR3, MR4, MR10A, MR10B, MR11A,
                         MR11B, MR12, MR13, MR15)
rm(PT1, PT2, PT3, PT4, PT5, PT6, PT7, PT8, PT9,
   MR1, MR2, MR3, MR4, MR10A, MR10B, MR11A,
   MR11B, MR12, MR13, MR15)
proportions <- pivot_longer(proportions, 2:7, names_to = "subcluster", values_to = "proportion")
ggplot(proportions, aes(x = Patient, y = proportion, fill = subcluster))+
  geom_bar(stat = "identity", position = "stack")+
  theme_minimal()+
  scale_fill_brewer(palette = "YlGnBu")


#for all cell types
#primary samples data (n=9)
PT1 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT1/abundance.csv", "PT1")
PT2 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT2/abundance.csv", "PT2")
PT3 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT3/abundance.csv", "PT3")
PT4 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT4/abundance.csv", "PT4")
PT5 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT5/abundance.csv", "PT5")
PT6 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT6/abundance.csv", "PT6")
PT7 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT7/abundance.csv", "PT7")
PT8 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT8/abundance.csv", "PT8")
PT9 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/PT9/abundance.csv", "PT9")

#metastasis
MR1 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR1/abundance.csv", "MR1")
MR2 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR2/abundance.csv", "MR2")
MR3 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR3/abundance.csv", "MR3")
MR4 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR4/abundance.csv", "MR4")
MR10A <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR10A/abundance.csv", "MR10A")
MR10B <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR10B/abundance.csv", "MR10B")
MR11A <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR11A/abundance.csv", "MR11A")
MR11B <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR11B/abundance.csv", "MR11B")
MR12 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR12/abundance.csv", "MR12")
MR13 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR13/abundance.csv", "MR13")
MR15 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR15/abundance.csv", "MR15")
MR18 <- spatial_celltype_proportions_df("/Users/alfred/Downloads/STOsteosarcoma/MR18/abundance.csv", "MR18")

#merging into single dataframe
proportions <- bind_rows(PT1, PT2, PT3, PT4, PT5, PT6, PT7, PT8, PT9,
                         MR1, MR2, MR3, MR4, MR10A, MR10B, MR11A,
                         MR11B, MR12, MR13, MR15)
rm(PT1, PT2, PT3, PT4, PT5, PT6, PT7, PT8, PT9,
   MR1, MR2, MR3, MR4, MR10A, MR10B, MR11A,
   MR11B, MR12, MR13, MR15, MR18)
proportions <- pivot_longer(proportions, 2:21, names_to = "subcluster", values_to = "proportion")
ggplot(proportions, aes(x = Patient, y = proportion, fill = subcluster))+
  geom_bar(stat = "identity", position = "stack")+
  theme_minimal()
