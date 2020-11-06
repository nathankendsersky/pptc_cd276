#### load packages ####

library(ggplot2)
library(data.table)
library(tidyverse)
library(cowplot)
setwd("./") # change to github working directory 

#### figure 1 ####

# figure 1a
## Treehouse Tumors vs. GTEx, CD276 TPM
## Treehouse: https://treehousegenomics.soe.ucsc.edu/public-data/#tumor_v11_polyA
### download Clinical Data and TPM Expression [log2(TPM+1)] from April 2020, extract CD276 expresion, convert to TPM
## GTEx: https://gtexportal.org/home/datasets
### download Gene TPMs, GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
### download Sample Annotations: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

# load Treehouse/GTEx merged tpm df
tree.gtex.tpm <- read.delim("data/Treehouse-GTEx-CD276-tpm.txt",header = T,sep = "\t")
# set order of Histologies/Tissues in plot
tree.gtex.tpm$Histology <- factor(tree.gtex.tpm$Histology, levels = c("Ewing sarcoma","Rhabdoid tumor","RMS","Embryonal RMS","Alveolar RMS","Hepatoblastoma",
                                                                      "Meningioma","Neuroblastoma","Osteosarcoma","Wilms tumor",
                                                                      "Adipose Tissue","Adrenal Gland",
                                                                      "Blood","Blood Vessel","Bone Marrow","Bladder","Brain","Breast",
                                                                      "Colon","Esophagus","Fallopian Tube","Heart","Kidney","Liver","Lung",
                                                                      "Muscle","Nerve","Ovary","Pancreas","Pituitary",
                                                                      "Prostate","Small Intestine","Salivary Gland","Skin","Stomach","Spleen",
                                                                      "Testis","Thyroid","Uterus","Cervix Uteri","Vagina"),order=T)

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologys/Tissues in plot
histo.colors.plot1a <- histo.colors %>%
  dplyr::filter(Histology %in% tree.gtex.tpm$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig1a <- ggplot(tree.gtex.tpm.color, aes(x=Histology,y=CD276_tpm,fill=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.size=0.5) + 
  theme_bw() + 
  scale_fill_manual(values = histo.colors.plot1a) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_blank(), axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "none") + 
  ylab("CD276 TPM")
fig1a


# figure 1b
## PPTC,CD276 TPM
### data from Rokita et al. 2019
### https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147
### download 2019-02-15-pptc_rnaseq_hg38_matrix_244.RData; TPM

# load PPTC tpm df
pptc.tpm <- read.delim("data/PPTC-CD276-tpm.txt",header = T,sep = "\t")

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")

# extract histology colors for Histologys/Tissues in plot
histo.colors.plot1b <- histo.colors %>%
  dplyr::filter(Histology %in% pptc.tpm$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig1b <- ggplot(pptc.tpm, aes(x=Histology,y=CD276_tpm,fill=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.shape = NA) + 
  theme_bw() + 
  scale_fill_manual(values = histo.colors.plot1b) +
  geom_jitter(shape=16, position=position_jitter(0.2),size=1,aes(color=InStudy),) + 
  scale_color_manual(values = c("black",adjustcolor("lightgray",alpha.f=0.5))) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_blank(), axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "none") + 
  ylab("CD276 TPM")
fig1b


#### figure 2 ####

# figure 2a

# load H-Score values from NBL TMA
nbl.ihc <- read.delim("data/PPTC-NBL-CD276-IHC.txt",header = T,sep = "\t")

# set NBL order from low to high
nbl.order <- nbl.ihc %>%
  group_by(SAMPFACTOR) %>%
  summarise(mean=mean(H.Score)) %>%
  arrange(mean) %>%
  pull(SAMPFACTOR)
nbl.ihc$SAMPFACTOR <- factor(nbl.ihc$SAMPFACTOR,levels=nbl.order)

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologys/Tissues in plot
histo.colors.plot2a <- histo.colors %>%
  dplyr::filter(Histology %in% nbl.ihc$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig2a <- ggplot(nbl.ihc,aes(x=SAMPFACTOR,y=H.Score,color=InStudy,fill=Histology)) + 
  geom_boxplot(color="black",outlier.size = 0.5) + 
  geom_point(size=1,aes(color=InStudy)) +
  scale_color_manual(values = c("black",adjustcolor("darkgray",alpha.f=0.5))) +
  scale_fill_manual(values = histo.colors.plot2a) +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(axis.text.x = element_text(size = 7, angle = 60, hjust=1, color="black"),
        axis.title.y = element_text(size=7),axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 7,colour = "black"),legend.position = "none")
fig2a

# figure 2b


