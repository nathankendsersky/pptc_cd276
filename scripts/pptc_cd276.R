#### load packages ####

library(ggplot2)
library(data.table)
library(tidyverse)
library(cowplot)
setwd("./")

#### figure 1 ####

tree.gtex.tpm <- read.delim("data/Treehouse-GTEx-CD276-tpm.txt",header = T,sep = "\t")
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")

tree.gtex.tpm$Histology <- factor(tree.gtex.tpm$Histology, levels = c("Ewing sarcoma","Rhabdoid tumor","RMS","Embryonal RMS","Alveolar RMS","Hepatoblastoma",
                                                                      "Meningioma","Neuroblastoma","Osteosarcoma","Wilms tumor",
                                                                      "Adipose Tissue","Adrenal Gland",
                                                                      "Blood","Blood Vessel","Bone Marrow","Bladder","Brain","Breast",
                                                                      "Colon","Esophagus","Fallopian Tube","Heart","Kidney","Liver","Lung",
                                                                      "Muscle","Nerve","Ovary","Pancreas","Pituitary",
                                                                      "Prostate","Small Intestine","Salivary Gland","Skin","Stomach","Spleen",
                                                                      "Testis","Thyroid","Uterus","Cervix Uteri","Vagina"),order=T)
histo.colors.plot1a <- histo.colors %>%
  dplyr::filter(Histology %in% tree.gtex.tpm$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig1a <- ggplot(tree.gtex.tpm.color, aes(x=Histology,y=CD276_tpm,fill=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.size=0.5) + 
  theme_bw() + 
  scale_fill_manual(values = histo.colors.plot1a) +
  #geom_jitter(shape=16,position = position_jitter(0.2),alpha=0.8) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_blank(), axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "none") + 
  ylab("CD276 TPM")
fig1a

