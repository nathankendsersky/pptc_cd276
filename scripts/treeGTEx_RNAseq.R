#### workflow ####

## Treehouse Tumors vs. GTEx, CD276 TPM
## Treehouse: https://treehousegenomics.soe.ucsc.edu/public-data/#tumor_v11_polyA
### download Clinical Data and TPM Expression [log2(TPM+1)] or counts from April 2020, extract CD276 expresion, convert to TPM
## GTEx: https://gtexportal.org/home/datasets
### download Gene TPMs or counts, GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
### download Sample Annotations: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

# NOTE: if you don't have access to a computer with >64GB RAM, you can download the count files and merge)
# the following code contains portions that were submitted to a computing cluster (marked), this will be commented out

## Computing Cluster; treeGTEx_RNA_qsub.R (>64GB RAM)
# 1. Prepare DGEList, filter, TMM-normalize
## Local Computer
# 2. Subset for CD276, clean Histology; treeGTEx_RNAseq.R
# 3. Visualize in main_pptc_cd276.R

#### prepare DGEList, filter, TMM-normalize ####

# source("./scripts/treeGTEx_RNAseq_qsub.R")
## this requires >64GB of RAM of load / analyze these datasets

#### load packages and functions ####

cat(paste0(Sys.time(),": Loading R packages..."),file=logFile1,append=T,sep="\n")

library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tidyselect)
library(data.table)
library(dplyr)
library(edgeR)
library(limma)
library(ggplot2)
if (!requireNamespace("tximport",quietly = T))
  BiocManager::install("tximport")
library(tximport) # package for getting Kallisto results into R
if (!requireNamespace("ensembldb",quietly = T))
  BiocManager::install("ensembldb")
library(ensembldb) #helps deal with ensembl
if (!requireNamespace("EnsDb.Hsapiens.v86",quietly = T))
  BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package

'%ni%' <- Negate('%in%')
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

cat(paste0(Sys.time(),":   Successfully loaded R packages."),file=logFile1,append=T,sep="\n")

#### load data from qsub script ####

# load RNAseq data (log2TPM tmm normalized and filtered)
load(file=paste0(checkpoints.out,"myDGEList.filtered.norm.RData"))
cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=F)
cpm.filtered.norm.df <- as_tibble(cpm.filtered.norm, rownames = "geneID")
rm(cpm.filtered.norm)

# load annotation data
if (dataset == "treehouse") {
  tumor.annot <- read_tsv("/Volumes/target_nbl_ngs/Kendsersky/TumorCompendium_v11_PublicPolyA/clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv")
  ped.tumors <- c("Ewing sarcoma","rhabdoid tumor","rhabdomyosarcoma","alveolar rhabdomyosarcoma","embryonal rhabdomyosarcoma",
                  "hepatoblastoma","neuroblastoma","meningioma","osteosarcoma","wilms tumor")
  tumor.annot.sub <- tumor.annot[tumor.annot$disease %in% ped.tumors,]
  tumor.annot.sub.slim <- tumor.annot.sub[,c("th_sampleid","disease")]
  colnames(tumor.annot.sub.slim) <- c("SAMPID","Histology")
  
}

gtex.annot <- read_tsv("/Volumes/target_nbl_ngs/Kendsersky/GTEx_Analysis_2017-06-05_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
gtex.annot.slim <- gtex.annot[,c("SAMPID","SMTS")]
colnames(gtex.annot.slim) <- c("SAMPID","Histology")
tumor.gtex.annot <- rbind(tumor.annot.sub.slim,gtex.annot.slim)

# add gene symbols to these matrices
colnames(cpm.filtered.norm.df)[1] <- "ID"
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys=as.character(cpm.filtered.norm.df$ID), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneIDs <- geneIDs %>% dplyr::rename(ID = GENEID)

cpm.filtered.norm.symbol.df <- merge(cpm.filtered.norm.df,geneIDs,by="ID")
cpm.filtered.norm.symbol.df <- as_tibble(cpm.filtered.norm.symbol.df[,c(ncol(cpm.filtered.norm.symbol.df),2:(ncol(cpm.filtered.norm.symbol.df)-1))])

#### filter and clean dataframe ####

tree.gtex.cd276 <- cpm.filtered.norm.symbol.df %>%
  dplyr::filter(SYMBOL=="CD276")
tree.gtex.names <- as.vector(colnames(tree.gtex.cd276[2:length(tree.gtex.cd276)]))
tree.gtex.cd276 <- pivot_longer(tree.gtex.cd276,
                                    cols = as.name(tree.gtex.names[1]):as.name(tree.gtex.names[length(tree.gtex.names)]),
                                    names_to = "SAMPID",
                                    values_to = "CD276_tpm") 

tree.gtex.cd276 <- tree.gtex.cd276 %>%
  dplyr::select(SAMPID,CD276_tpm)

tree.gtex.cd276 <- merge(tree.gtex.cd276,tumor.gtex.annot,by="SAMPID")

# clean histology names 
tree.gtex.cd276$Histology <- gsub("rhabdoid tumor","Rhabdoid tumor",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("alveolar rhabdomyosarcoma","Alveolar RMS",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("embryonal rhabdomyosarcoma","Embryonal RMS",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("rhabdomyosarcoma","RMS",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("hepatoblastoma","Hepatoblastoma",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("neuroblastoma","Neuroblastoma",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("meningioma","Meningioma",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("osteosarcoma","Osteosarcoma",tree.gtex.cd276$Histology)
tree.gtex.cd276$Histology <- gsub("wilms tumor","Wilms tumor",tree.gtex.cd276$Histology)


