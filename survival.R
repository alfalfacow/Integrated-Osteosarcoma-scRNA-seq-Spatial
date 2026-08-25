#Survival (Bulk RNAseq) Validation
library(survminer)
library(survival)
library(data.table)
library(tidyverse)

RNAseq<-fread("/Users/alfred/Downloads/TARGET-OS.star_counts.tsv.gz") #Accessed via UCSC Xena Browser: https://xenabrowser.net/datapages/?cohort=GDC%20TARGET-OS&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

#Converting ensembl IDs into gene names using biomaRt
library(biomaRt)
#strip the suffix
RNAseq$Ensembl_ID <- sub("\\..*", "", RNAseq$Ensembl_ID)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                   host = "https://useast.ensembl.org")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = RNAseq$Ensembl_ID,   # your column of Ensembl IDs
  mart = ensembl
)
RNAseq <- merge(RNAseq, gene_map, by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE) #+1 row because 2 forms of LINC mapped to same ID, 29241
genes <- as.data.frame(t(RNAseq))
gene_names <- genes[90,]
genes <- genes[-c(1,90),]
colnames(genes) <- gene_names
genes_of_interest <- c( #subcluster0
  "SPP1", "MMP13", "SFRP4",
  "CLU","IBSP", "SLC29A1"
)
genes_of_interest <- c( #subcluster2
  "COL2A1", "S100A1", "SNORC",
  "COL9A3", "HAPLN1", "GSTA1",
  "TNNT3", "BNIP3", "SCIN", "FGFBP2", "COL11A2",
  "ENO2", "VEGFA"
)
genes_of_interest <- c( #subcluster3
  "OASL", "C5orf46", "RSAD2",
  "IFIT2","IFIT1", "IFIT3",
  "ISG15", "TNFSF10", "OAS1", "MX1"
)
Gene_set <- list( #"GSVA expects list of gene sets"
  GSVA_sig = genes_of_interest
)
#gsva for scoring
library(GSVA)
expr <- as.matrix(genes)
storage.mode(expr) <- "numeric"
expr <- t(expr)

#deleting the dupes
expr_df <- as.data.frame(expr)
expr_df$gene <- rownames(expr_df)

library(dplyr)
expr_collapsed <- expr_df %>%
  group_by(gene) %>%
  summarise(across(everything(), max)) %>%   # or mean, depending on your preference
  as.data.frame()

rownames(expr_collapsed) <- expr_collapsed$gene
expr_collapsed$gene <- NULL
expr <- as.matrix(expr_collapsed)

#continue gsva
param <- ssgseaParam(
  exprData = expr,
  geneSets = Gene_set
)
gsva_result <- as.data.frame(t(gsva(param)))
gsva_result$cases.submitter_id <- substr(rownames(gsva_result), 1, nchar(rownames(gsva_result))-4)

#Loading clinical data
clinical<-fread("/Users/alfred/Downloads/clinical.cart.2026-08-20/clinical.tsv") #accessed from GDC portal
clinical <- clinical[!duplicated(clinical$cases.submitter_id), ] #remove duplicated patients
clinical <- clinical[, c(10,11,18,26,28,33,64, 112)] 
genes2 <- genes
genes2$cases.submitter_id<- substr(rownames(genes2), 1, nchar(rownames(genes2))-4)

# Remove invalid column names
keep <- !is.na(colnames(genes2)) & colnames(genes2) != ""
genes2 <- genes2[, keep]

# Ensure all column names are unique
colnames(genes2) <- make.unique(colnames(genes2))

# Verify
anyDuplicated(colnames(genes2))
sum(is.na(colnames(genes2)))
sum(colnames(genes2) == "")

#merge data
survival_data <- inner_join(clinical, genes2, by = "cases.submitter_id") #merge survival data with genes
survival_data <- inner_join(survival_data, gsva_result, by = "cases.submitter_id") #merge survival data with GSVA scores per patient

#Censoring Data
survival_data$deceased <- ifelse(survival_data$demographic.vital_status=="Alive", FALSE, TRUE)

#overall survival
survival_data$overall_survival <- ifelse(survival_data$demographic.vital_status=="Alive", 
                                                as.numeric(survival_data$diagnoses.days_to_last_follow_up), 
                                                as.numeric(survival_data$demographic.days_to_death))

#Patient stratification
survival_data$GSVA_sig <- as.numeric(survival_data$GSVA_sig) #scores merged as factor class, transform into numeric
tertiles<-quantile(survival_data$GSVA_sig, probs=c(0.33, 0.66))
for(i in 1:88){
  if(survival_data$GSVA_sig[i]<=tertiles[1]){
    survival_data$strata[i]<-"LOW"
  }else if(survival_data$GSVA_sig[i]>=tertiles[2]){
    survival_data$strata[i]<-"HIGH"
  }else{
    survival_data$strata[i]<-NA
  }
}

#creating the kaplan meier survival curve
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data=survival_data) #high vs medium vs low tert expression
ggsurvplot(fit, data=survival_data,  pval=T, risk.table=T, conf.int=F)
