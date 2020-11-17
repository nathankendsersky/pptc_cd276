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
#### prepare DGEList, filter, TMM-normalize ####

# source("./scripts/treeGTEx_RNAseq_qsub.R")
## this requires >64GB of RAM of load / analyze these datasets

#### load data from qsub script ####

dataset.dir <- "treehouse_data/"
analysis.dir <- "solidtumor_analysis/"
home <- "/Volumes/target_nbl_ngs/Kendsersky/crc-project/tumor_datasets/"
checkpoints.out <- paste0(home,dataset.dir,analysis.dir,"checkpoints_out/")
resources.dir <- "/Volumes/target_nbl_ngs/Kendsersky/crc-project/resources/"
dataset <- gsub("_data/","",dataset.dir)
tumor <- "multi-tumor"

# load RNAseq data (log2TPM tmm normalized and filtered)
load(file=paste0(checkpoints.out,"myDGEList.filtered.norm.RData"))
tree.gtex.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=F)
tree.gtex.cpm.filtered.norm.df <- as_tibble(tree.gtex.cpm.filtered.norm, rownames = "GENEID")

# load annotation data, treehouse
tumor.annot <- read_tsv("/Volumes/target_nbl_ngs/Kendsersky/TumorCompendium_v11_PublicPolyA/clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv")
ped.tumors <- c("Ewing sarcoma","rhabdoid tumor","rhabdomyosarcoma","alveolar rhabdomyosarcoma","embryonal rhabdomyosarcoma",
                "hepatoblastoma","neuroblastoma","meningioma","osteosarcoma","wilms tumor")
tumor.annot.sub <- tumor.annot[tumor.annot$disease %in% ped.tumors,]
tumor.annot.sub.slim <- tumor.annot.sub[,c("th_sampleid","disease")]
colnames(tumor.annot.sub.slim) <- c("SAMPID","Histology")

# load annotation data, gtex
gtex.annot <- read_tsv("/Volumes/target_nbl_ngs/Kendsersky/GTEx_Analysis_2017-06-05_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
gtex.annot.slim <- gtex.annot[,c("SAMPID","SMTS")]
colnames(gtex.annot.slim) <- c("SAMPID","Histology")
tumor.gtex.annot <- rbind(tumor.annot.sub.slim,gtex.annot.slim)

# add gene symbols to these matrices
gene.convert <- ensembldb::select(EnsDb.Hsapiens.v86, keys=as.character(tree.gtex.cpm.filtered.norm.df$GENEID), 
                                  keytype = "GENEID", columns = c("SYMBOL","GENEID"))

tree.gtex.cpm.symbol.df <- merge(tree.gtex.cpm.filtered.norm.df,gene.convert,by="GENEID")
tree.gtex.cpm.symbol.df <- as_tibble(tree.gtex.cpm.symbol.df[,c(ncol(tree.gtex.cpm.symbol.df),2:(ncol(tree.gtex.cpm.symbol.df)-1))])

#### filter and clean dataframe ####

first.samp <- as.name(colnames(tree.gtex.cpm.symbol.df[2]))
last.samp <- as.name(colnames(tree.gtex.cpm.symbol.df[ncol(tree.gtex.cpm.symbol.df)]))

tree.gtex.cd276 <- tree.gtex.cpm.symbol.df %>%
  dplyr::filter(SYMBOL=="CD276") %>%
  pivot_longer(cols = first.samp:last.samp,
               names_to = "SAMPID",
               values_to = "TPM") %>%
  dplyr::select(SAMPID,TPM) %>%
  dplyr::rename(CD276_tpm=TPM)

# merge annotations
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

# save table
write.table(tree.gtex.cd276,"./data/RNA_data/Treehouse-GTEx-CD276-tpm.txt",sep = "\t", col.names = T, row.names = F, quote = F)