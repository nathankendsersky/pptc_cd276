#### workflow ####

# Data from Rokita et al. 2019
# RNA-seq data present on dbGaP (phs001437.v1.p1)
# Processed FPKM data at: https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147

# 1. Run kallisto on fastq.gz files; runKallisto.sh
# 2. Run tximport; pptc_RNAseq.R (chunk below)
# 3. Prepare DGEList, filter, TMM-normalize; pptc_RNAseq.R (chunk below)
# 4. Subset for CD276, clean Histology; pptc_RNAseq.R (chunk below, START HERE with RData file)
# 5. Visualize in main_pptc_cd276.R

#### load packages ####

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
setwd(".") # change to github working directory 

#### prepare annotation file ####

# manifest is all files in the kallisto output directory (aln)
manifest <- list.files(path="/Volumes/target_nbl_ngs/Kendsersky/b7h3/rnaseq_data/aln/",pattern=".aligned")

# declare tumor types of interest
solid.tumors <- c("ATRT","Ewing Sarcoma","Extracranial Rhabdoid","Fusion- RMS","Fusion+ RMS",
                  "Hepatoblastoma","Neuroblastoma","Osteosarcoma","Wilms")

# file from Rokita et al, 2019
pptc.annot.full <- read_delim(url("https://ndownloader.figshare.com/files/16603061"),delim = "\t",)
# filter for solid tumors with RNA data
pptc.annot <- pptc.annot.full %>%
  dplyr::filter(Have.fpkm.file == "yes", RNA.human.bam.filename %in% gsub(".aligned",".bam",manifest)) %>%
  dplyr::select(Model,Histology.Detailed,RNA.human.bam.filename) %>%
  dplyr::filter(Histology.Detailed %in% solid.tumors) %>%
  dplyr::rename(SAMPID=Model,Histology=Histology.Detailed)

# change case in names
pptc.annot$Histology <- gsub("Ewing Sarcoma","Ewing sarcoma",pptc.annot$Histology)
pptc.annot$Histology <- gsub("Wilms","Wilms tumor",pptc.annot$Histology)
pptc.annot$Histology[which(pptc.annot$SAMPID=="KT-12")] <- "Extracranial Rhabdoid" # change classification of KT-12

#### Run tximport ####

# from runKallisto.sh output, load the abundance.tsv file for each sample, fun tximport

path <- file.path(paste0("/Volumes/target_nbl_ngs/Kendsersky/b7h3/rnaseq_data/aln/",
                         gsub(".bam",".aligned",pptc.annot$RNA.human.bam.filename)),"abundance.tsv")
all(file.exists(path))
# set file paths to your mapped data
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
beep(2)

#### prepare dgelist, filter, TMM-normalized ####

pptc.names <- as.vector(pptc.annot$SAMPID)
pptc.myDGEList <- DGEList(Txi_gene$counts)

# use cpm to filter
pptc.cpm <- cpm(pptc.myDGEList)
pptc.min.histo <- min(table(pptc.annot$Histology)) # one Hepatoblastoma sample
pptc.keepers <- rowSums(pptc.cpm>1)>=pptc.min.histo
pptc.myDGEList.filtered <- pptc.myDGEList[pptc.keepers,]

# normalize with TMM
pptc.myDGEList.filtered.norm <- calcNormFactors(pptc.myDGEList.filtered, method = "TMM")
pptc.cpm.filtered.norm <- cpm(pptc.myDGEList.filtered.norm, log=FALSE)
pptc.cpm.filtered.norm.df <- tibble::rownames_to_column(data.frame(pptc.cpm.filtered.norm),"SYMBOL")
colnames(pptc.cpm.filtered.norm.df) <- c("SYMBOL", pptc.names)

#### subset for CD276, merge with pptc.annot ####

cpm.filtered.norm.df.pivot <- pivot_longer(pptc.cpm.filtered.norm.df,
                                           cols = as.name(pptc.names[1]):as.name(pptc.names[length(pptc.names)]),
                                           names_to = "SAMPID",
                                           values_to = "TPM")

pptc.cd276 <- cpm.filtered.norm.df.pivot %>%
  dplyr::filter(SYMBOL=="CD276") %>% 
  dplyr::select(SAMPID,TPM) %>%
  dplyr::rename(CD276_tpm=TPM)

pptc.cd276 <- merge(pptc.cd276,pptc.annot[,1:2],by="SAMPID")

# load samples from study
in.study <- scan(file="./data/Models-in-Study.txt",what=character())
pptc.cd276$InStudy <- ifelse(pptc.cd276$SAMPID %in% in.study,"In Study","Not In Study")

# save table
write.table(pptc.cd276,"./data/RNA_data/PPTC-CD276-tpm.txt",sep = "\t", col.names = T, row.names = F, quote = F)
