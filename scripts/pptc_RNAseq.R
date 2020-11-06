#### workflow ####

# Data from Rokita et al. 2019
# RNA-seq data present on dbGaP (phs001437.v1.p1)
# Processed FPKM data at: https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147

# 1. Run kallisto on fastq.gz files; runKallisto.sh
# 2. Run tximport; pptc_RNAseq.R (chunk below)
# 3. Prepare DGEList, filter, TMM-normalize; pptc_RNAseq.R (chunk below)
# 4. Subset for CD276, clean Histology; pptc_RNAseq.R (chunk below)
# 5. Visualize in main_pptc_cd276.R

#### load packages ####

library(ggplot2)
library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
setwd("./") # change to github working directory 

#### Run tximport ####

# from runKallisto, load the abundance.tsv file for each sample, fun tximport
path <- file.path(paste0("/Volumes/target_nbl_ngs/Kendsersky/b7h3/rnaseq_data/aln/",manifest), "abundance.tsv")
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


tumor.counts <- tibble::rownames_to_column(as.data.frame(Txi_gene$counts),"Gene")
sample.names <- as.vector(targets$sample)
colnames(tumor.counts) <- c("Gene",sample.names)



pptc.rna <- readRDS("")
