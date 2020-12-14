#### load packages ####

library(ggplot2)
library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
setwd("./") # change to github working directory 

#### figure 1a ####

## Treehouse Tumors vs. GTEx, CD276 TPM
## Treehouse: https://treehousegenomics.soe.ucsc.edu/public-data/#tumor_v11_polyA
### download Clinical Data and TPM Expression [log2(TPM+1)] from April 2020, extract CD276 expresion, convert to TPM
## GTEx: https://gtexportal.org/home/datasets
### download Gene TPMs, GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
### download Sample Annotations: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

# load Treehouse/GTEx merged tpm df
tree.gtex.tpm <- read.delim("data/RNA_data/Treehouse-GTEx-CD276-tpm.txt",header = T,sep = "\t")
# set order of Histologies/Tissues in plot
tree.gtex.tpm$Histology <- factor(tree.gtex.tpm$Histology, levels = c("Ewing sarcoma","Rhabdoid tumor","RMS","Embryonal RMS","Alveolar RMS","Hepatoblastoma",
                                                                      "Meningioma","Neuroblastoma","Osteosarcoma","Wilms tumor",
                                                                      "Adipose Tissue","Adrenal Gland",
                                                                      "Blood","Blood Vessel","Bone Marrow","Bladder","Brain","Breast",
                                                                      "Colon","Esophagus","Fallopian Tube","Heart","Kidney","Liver","Lung",
                                                                      "Muscle","Nerve","Ovary","Pancreas","Pituitary",
                                                                      "Prostate","Small Intestine","Salivary Gland","Skin","Stomach","Spleen",
                                                                      "Testis","Thyroid","Uterus","Cervix Uteri","Vagina"),order=T)
tree.gtex.tpm$Dataset <- factor(tree.gtex.tpm$Dataset, levels = c("Human Tumors (TCCI)", "Human Normal (GTEx)"))


# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")

# extract histology colors for Histologys/Tissues in plot
histo.colors.plot1a <- histo.colors %>%
  dplyr::filter(Histology %in% tree.gtex.tpm$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig1a <- ggplot(tree.gtex.tpm, aes(x=Histology,y=CD276_tpm,fill=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.size=0.5) + 
  facet_grid(cols = vars(Dataset), scales = "free", space = "free") +
  theme_bw() + 
  scale_fill_manual(values = histo.colors.plot1a) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_blank(), axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "none") + 
  ylab("CD276 Transcripts per Million (tpm)")
fig1a

#### figure 1b ####

## PPTC,CD276 TPM
### data from Rokita et al. 2019
### FPKM data available at https://figshare.com/projects/Genomic_landscape_of_childhood_cancer_patient-derived_xenograft_models/38147
### Analysis: 

# load PPTC tpm df
pptc.tpm <- read.delim("data/RNA_data/PPTC-CD276-tpm.txt",header = T,sep = "\t")

# set order for Histologies in plot
pptc.tpm$Histology <- factor(pptc.tpm$Histology, levels = c("ATRT", "Ewing sarcoma", "Extracranial Rhabdoid","Embryonal RMS",
                                                            "Alveolar RMS", "Hepatoblastoma","Neuroblastoma","Osteosarcoma","Wilms tumor"))

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")

# extract histology colors for Histologys/Tissues in plot
histo.colors.plot1b <- histo.colors %>%
  dplyr::filter(Histology %in% pptc.tpm$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig1b <- ggplot(pptc.tpm, aes(x=Histology,y=CD276_tpm,fill=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.shape = NA) + 
  facet_grid(cols = vars(pptc.tpm$Dataset)) +
  theme_bw() + 
  scale_fill_manual(values = histo.colors.plot1b) +
  geom_jitter(shape=16, position=position_jitter(0.2),size=1,aes(color=InStudy),) + 
  scale_color_manual(values = c("black",adjustcolor("lightgray",alpha.f=0.5))) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_blank(), axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "none") + 
  ylab("CD276 Transcripts per Million (tpm)")
fig1b


#### figure 2a ####

# load H-Score values from NBL TMA
nbl.ihc <- read.delim("data/IHC_data/IHC-NBL-CD276.txt",header = T,sep = "\t")

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
  geom_boxplot(color="black",outlier.shape = NA) + 
  geom_point(size=0.5,aes(color=InStudy)) +
  scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = histo.colors.plot2a) +
  ggtitle("Neuroblastoma PDX TMA") +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 10,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 60, hjust=1, color="black"),
        axis.title.y = element_text(size=6),axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 6,colour = "black"),legend.position = "none")
fig2a

#### figure 2b ####

# load H-Score values from NBL TMA
os.ihc <- read.delim("data/IHC_data/IHC-OS-CD276.txt",header = T,sep = "\t")

# set NBL order from low to high
os.order <- os.ihc %>%
  group_by(SAMPFACTOR) %>%
  summarise(mean=mean(H.Score)) %>%
  arrange(mean) %>%
  pull(SAMPFACTOR)
os.ihc$SAMPFACTOR <- factor(os.ihc$SAMPFACTOR,levels=os.order)

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologys/Tissues in plot
histo.colors.plot2b <- histo.colors %>%
  dplyr::filter(Histology %in% os.ihc$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig2b <- ggplot(os.ihc,aes(x=SAMPFACTOR,y=H.Score,color=InStudy,fill=Histology)) + 
  geom_boxplot(color="black",outlier.shape = NA) + 
  geom_point(size=0.5,aes(color=InStudy)) +
  scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = histo.colors.plot2b) +
  ggtitle("Osteosarcoma PDX TMA") +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 10,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 60, hjust=1, color="black"),
        axis.title.y = element_text(size=6),axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 6,colour = "black"),legend.position = "none")
fig2b

#### figure 2c ####

# load H-Score values from NBL TMA
mix.ihc <- read.delim("data/IHC_data/IHC-MIXED-CD276.txt",header = T,sep = "\t")

# set NBL order from low to high
mix.order <- mix.ihc %>%
  group_by(SAMPFACTOR) %>%
  summarise(mean=mean(H.Score)) %>%
  arrange(mean) %>%
  pull(SAMPFACTOR)
mix.ihc$SAMPFACTOR <- factor(mix.ihc$SAMPFACTOR,levels=mix.order)

mix.ihc$Histology <- factor(mix.ihc$Histology, levels = c("ATRT", "Ewing sarcoma", "Extracranial Rhabdoid","Embryonal RMS",
                                                          "Alveolar RMS", "Hepatoblastoma","Neuroblastoma","Osteosarcoma","PXA","Wilms tumor"))
mix.ihc$Histology.Short <- factor(mix.ihc$Histology.Short, levels = c("ATRT", "Ewing sarcoma", "Rhabdoid","ERMS",
                                                                      "ARMS", "Hepatoblastoma","Neuroblastoma","Osteosarcoma","PXA","Wilms tumor"))

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologys/Tissues in plot
histo.colors.plot2c <- histo.colors %>%
  dplyr::filter(Histology %in% mix.ihc$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig2c <- ggplot(mix.ihc,aes(x=SAMPFACTOR,y=H.Score,color=InStudy,fill=Histology)) + 
  geom_boxplot(color="black",outlier.shape = NA) + 
  geom_point(size=0.5,aes(color=InStudy)) +
  scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = histo.colors.plot2c) +
  facet_grid(cols = vars(mix.ihc$Histology.Short), scales = "free", space = "free") +
  ggtitle("Mixed Histology PDX TMA") +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 10,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 60, hjust=1, color="black"),
        axis.title.y = element_text(size=6),axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 6,colour = "black"),legend.position = "none")
fig2c

#### figure 3a ####

# load N>1 relative tumor volume
rtv <- read.delim("data/PPTC-Ngreater1-RTV.txt",header = T,sep = "\t")
# calculate % change in minRTV
rtv$percent_minRTV <- round((rtv$minRTV*100)-100,digits = 2)

# set order high to low
rtv.order <- rtv %>%
  group_by(SAMPID, Histology) %>%
  summarise(med_percent_minRTV=median(percent_minRTV)) %>%
  arrange(desc(med_percent_minRTV),Histology) %>%
  pull(SAMPID) %>%
  unique()
rtv$SAMPID <- factor(rtv$SAMPID,levels = rtv.order)


# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologys/Tissues in plot
histo.colors.plot3a <- histo.colors %>%
  dplyr::filter(Histology %in% rtv$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()


fig3a <- ggplot(rtv, aes(x=SAMPID, y=percent_minRTV, fill=Histology)) + 
  theme_bw() +   
  stat_summary(fun = median, geom="bar",color="black") + 
  scale_fill_manual(values = histo.colors.plot3a) + 
  geom_point() +
  facet_grid(cols = vars(rtv$Histology), scales = "free", space = "free") +
  ylab("% Change in Minimum Relative Tumor Volume (RTV)") + xlab("Tumor Model") +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=6),axis.title.x=element_text(size=8),
        axis.text.y = element_text(size = 6,colour = "black"),
        legend.title = element_text(size=6), legend.text = element_text(size=6),
        legend.position = c(0.6, 0.75),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
fig3a

#### save figures as PDFs ####

pdf("figures/pptc_cd276_Figure1.pdf",width = 8.5,height = 4)
plot_grid(fig1a,fig1b, labels=c("A","B"), align = 'h', nrow = 1, rel_widths = c(1/1,1/3), label_size=12)
dev.off()

pdf("figures/pptc_cd276_Figure2.pdf",width = 11,height = 8.5)
ggdraw() +
  draw_plot(fig2c, 0,.67,1,.33) +
  draw_plot(fig2a,0,.25,1,.42) +
  draw_plot(fig2b, 0,0,.5,.25) +
  draw_plot_label(c("A","B","C"), c(0, 0, 0), c(0, 0.67, 0.25), size = 12)
dev.off()

pdf("figures/pptc_cd276_Figure3.pdf",width = 5,height = 5)
plot_grid(fig3a, labels=c("A"), align = 'h', nrow = 1, rel_widths = c(1), label_size=12)
dev.off()



