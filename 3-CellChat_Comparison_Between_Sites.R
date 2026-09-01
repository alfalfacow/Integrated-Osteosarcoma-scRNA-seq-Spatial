##Cell-cell communication comparison using CellChat
##Primary and metastatic CellChat objects were compared following the
##official CellChat multiple-dataset comparison tutorial:
#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html

library(CellChat)

primary_CellChat <- readRDS("primary_CellChat.rds")
metastatic_CellChat <- readRDS("metastatic_CellChat.rds")

object.list <- list(
  Primary = primary_CellChat,
  Metastatic = metastatic_CellChat
)

##Subsequent comparison of interaction number, interaction strength,
##signaling roles, pathway information flow, and pathway-specific
##networks followed the CellChat multiple-dataset comparison tutorial.
