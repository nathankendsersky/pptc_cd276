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
tree.gtex.tpm$Dataset <- gsub(" \\(TCCI\\)","",tree.gtex.tpm$Dataset)
tree.gtex.tpm$Dataset <- factor(tree.gtex.tpm$Dataset, levels = c("Human Tumors", "Human Normal (GTEx)"))
# table(tree.gtex.tpm$Dataset)

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
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 8,colour = "black"),
        strip.text = element_text(size = 8),
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
  scale_color_manual(values = c("black",adjustcolor("lightgray",alpha.f=.5))) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 8,colour = "black"),
        strip.text = element_text(size = 8),
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
# extract histology colors for Histologies/Tissues in plot
histo.colors.plot2a <- histo.colors %>%
  dplyr::filter(Histology %in% nbl.ihc$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig2a <- ggplot(nbl.ihc,aes(x=SAMPFACTOR,y=H.Score,color=InStudy,fill=Histology)) + 
  geom_boxplot(color="black",outlier.shape = NA) + 
  geom_point(size=1,aes(color=InStudy)) +
  scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = histo.colors.plot2a) +
  facet_grid(cols = vars(nbl.ihc$Histology), scales = "free", space = "free") +
  #ggtitle("Neuroblastoma PDX TMA") +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 45, hjust=1, color="black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 8,colour = "black"),
        strip.text = element_text(size=8),
        legend.position = "none")
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
# extract histology colors for Histologiess/Tissues in plot
histo.colors.plot2b <- histo.colors %>%
  dplyr::filter(Histology %in% os.ihc$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

fig2b <- ggplot(os.ihc,aes(x=SAMPFACTOR,y=H.Score,color=InStudy,fill=Histology)) + 
  geom_boxplot(color="black",outlier.shape = NA) + 
  geom_point(size=1,aes(color=InStudy)) +
  scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = histo.colors.plot2b) +
  facet_grid(cols = vars(os.ihc$Histology), scales = "free", space = "free") +
  #ggtitle("Osteosarcoma PDX TMA") +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 45, hjust=1, color="black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 8,colour = "black"),
        strip.text = element_text(size=8),
        legend.position = "none")
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
  geom_point(size=1,aes(color=InStudy)) +
  scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = histo.colors.plot2c) +
  facet_grid(cols = vars(mix.ihc$Histology.Short), scales = "free", space = "free") +
  #ggtitle("Mixed Histology PDX TMA") +
  xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 45, hjust=1, color="black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 8,colour = "black"),
        strip.text = element_text(size=8),
        legend.position = "none")
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

# prepare factors
rtv$Histology <- factor(rtv$Histology, levels = c("ATRT", "Ewing sarcoma", "Extracranial Rhabdoid","Embryonal RMS",
                                                  "Alveolar RMS", "Hepatoblastoma","Neuroblastoma","Osteosarcoma","PXA","Wilms tumor"))
rtv$Histology.Short <- factor(rtv$Histology.Short, levels = c("ATRT", "EWS", "Rhabdoid","ERMS",
                                                                      "ARMS", "Hepatoblastoma","Neuroblastoma","Osteosarcoma","PXA","WT"))

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
  facet_grid(cols = vars(rtv$Histology.Short), scales = "free", space = "free") +
  ylab("% Change in minimum RTV") +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust=1, color="black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 6,colour = "black"),
        legend.title = element_text(size=6), 
        legend.text = element_text(size=6),
        strip.text = element_text(size=6),
        legend.position = "none",panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
fig3a

#### figure 3b ####

# load N=1 relative tumor volume
smt.rtv <- read.delim("data/PPTC-SMT-RTV.txt",header = T,sep = "\t")
# calculate % change in minRTV
smt.rtv$percent_minRTV <- round((smt.rtv$minRTV*100)-100,digits = 2)

# set order high to low
smt.rtv.order <- smt.rtv %>%
  group_by(SAMPID, Histology) %>%
  arrange(desc(percent_minRTV),Histology) %>%
  pull(SAMPID) %>%
  unique()
smt.rtv$SAMPID <- factor(smt.rtv$SAMPID,levels = smt.rtv.order)

smt.rtv$Histology <- factor(smt.rtv$Histology, levels = c("ATRT","Ewing sarcoma","Extracranial Rhabdoid","Embryonal RMS","Alveolar RMS","Hepatoblastoma",
                                                          "Meningioma","Neuroblastoma","Osteosarcoma","PXA","Wilms tumor"))
smt.rtv$Histology.Short <- factor(smt.rtv$Histology.Short, levels = c("ATRT","Ewing sarcoma","Rhabdoid","Embryonal RMS","Alveolar RMS","HB",
                                                                      "MT","NBL","OS","PXA","Wilms tumor"))

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologies/Tissues in plot
histo.colors.plot3b <- histo.colors %>%
  dplyr::filter(Histology %in% smt.rtv$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()


fig3b <- ggplot(smt.rtv, aes(x=SAMPID, y=percent_minRTV, fill=Histology)) + 
  theme_bw() +   
  stat_summary(fun = median, geom="bar",color="black") + 
  scale_fill_manual(values = histo.colors.plot3b) + 
  geom_point() +
  facet_grid(cols = vars(smt.rtv$Histology.Short), scales = "free", space = "free") +
  ylab("% Change in minimum RTV") +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust=1, color="black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 6,colour = "black"),
        legend.title = element_text(size=6), 
        legend.text = element_text(size=6),
        strip.text = element_text(size=6),
        legend.position = "none",panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
fig3b

#### supplemental figure 1 ####

# summarize IHC data
nbl.ihc.slim <- nbl.ihc %>%
  group_by(SAMPFACTOR,Histology) %>%
  summarise(medianHScore=median(H.Score)) %>%
  dplyr::rename(SAMPID=SAMPFACTOR)

mix.ihc.slim <- mix.ihc %>%
  group_by(SAMPFACTOR,Histology) %>%
  summarise(medianHScore=median(H.Score)) %>%
  dplyr::rename(SAMPID=SAMPFACTOR)

os.ihc.slim <- os.ihc %>%
  group_by(SAMPFACTOR,Histology) %>%
  summarise(medianHScore=median(H.Score)) %>%
  dplyr::rename(SAMPID=SAMPFACTOR)

# row bind IHC data
ihc.allsamp <- rbind(nbl.ihc.slim,rbind(mix.ihc.slim,os.ihc.slim))

# select cols in tpm df, then merge
pptc.tpm.slim <- pptc.tpm %>%
  dplyr::select(SAMPID,CD276_tpm,Histology)
ihc.tma.final <- merge(pptc.tpm.slim, ihc.allsamp, by=c("SAMPID","Histology"))

# perform correlation between RNA (tpm) and protein (H Score)
rna.pro.cor <- cor.test(ihc.tma.final$medianHScore, ihc.tma.final$CD276_tpm, method="pearson")

histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")
# extract histology colors for Histologys/Tissues in plot
histo.colors.sfig1 <- histo.colors %>%
  dplyr::filter(Histology %in% ihc.tma.final$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

sfig1 <- ggplot(ihc.tma.final, aes(x = medianHScore, y = CD276_tpm)) +
  geom_point(size=3,aes(color=factor(Histology))) +
  geom_smooth(method=lm,se=TRUE,color="black") +
  scale_color_manual(values=histo.colors.sfig1) + 
  xlab("CD276 Protein (Median H-Score)") + ylab("CD276 mRNA (Transcripts per Million, tpm)") +
  labs(title = "Correlation between CD276 mRNA and Protein",
       subtitle = paste0("Pearson Correlation: ",round(rna.pro.cor$estimate,3),
                         "\nP value: ", rna.pro.cor$p.value),
       color = "Histology") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8,color="black"),
        axis.title.y = element_text(size=8,color="black"),
        axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "bottom",legend.text = element_text(size = 8,colour = "black"))
sfig1

#### supplemental figure for overall response - RNA ####

obj.resp <- read_tsv("data/Objective-Responses.txt")
obj.resp$Objective_Response <- factor(obj.resp$Objective_Response,levels=c("PD1","SD","PR","CR","MCR"))

# merge with pptc.tpm
pptc.tpm.resp <- merge(pptc.tpm,obj.resp,by="SAMPID")

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")

# extract histology colors for Histologys/Tissues in plot
histo.colors.sfig.resp.a <- histo.colors %>%
  dplyr::filter(Histology %in% pptc.tpm.resp$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

pptc.tpm.resp$Histology <- factor(pptc.tpm.resp$Histology, levels = c("ATRT","Ewing sarcoma","Extracranial Rhabdoid","Embryonal RMS","Alveolar RMS","Hepatoblastoma",
                                                          "Meningioma","Neuroblastoma","Osteosarcoma","PXA","Wilms tumor"))

sfig.resp.a <- ggplot(pptc.tpm.resp, aes(x=Objective_Response,y=CD276_tpm,color=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.shape = NA) + 
  theme_bw() + 
  scale_color_manual(values = histo.colors.sfig.resp.a) +
  geom_jitter(shape=16, position=position_jitter(0.2),size=2) + 
  #scale_color_manual(values = c("red",adjustcolor("black",alpha.f=0.5))) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_text(size = 8,colour = "black"), 
        axis.text.y = element_text(size = 8,colour = "black"),
        legend.position = "none") + 
  ylab("CD276 Transcripts per Million (TPM)") + xlab("Objective Response Category")
sfig.resp.a

#### supplemental figure for overall response - protein ####

# merge with ihc.allsamp
ihc.allsamp.resp <- merge(ihc.allsamp,obj.resp,by="SAMPID")

# load Histology Colors
histo.colors <- read.delim("data/Histology-Colors.txt",header=T,sep = "\t")

# extract histology colors for Histologys/Tissues in plot
histo.colors.sfig.resp.b <- histo.colors %>%
  dplyr::filter(Histology %in% ihc.allsamp.resp$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

ihc.allsamp.resp$Histology <- factor(ihc.allsamp.resp$Histology, levels = c("ATRT","Ewing sarcoma","Extracranial Rhabdoid","Embryonal RMS","Alveolar RMS","Hepatoblastoma",
                                                                      "Meningioma","Neuroblastoma","Osteosarcoma","PXA","Wilms tumor"))

sfig.resp.b <- ggplot(ihc.allsamp.resp, aes(x=Objective_Response,y=medianHScore,color=Histology)) + 
  geom_boxplot(color="black",alpha=.9,outlier.shape = NA) + 
  theme_bw() + 
  scale_color_manual(values = histo.colors.sfig.resp.b) +
  geom_jitter(shape=16, position=position_jitter(0.2),size=2) + 
  #scale_color_manual(values = c("red",adjustcolor("black",alpha.f=0.5))) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1,colour = "black"),
        axis.title.y = element_text(size=8),axis.title.x=element_text(size = 8,colour = "black"), 
        axis.text.y = element_text(size = 8,colour = "black")) + 
  ylab("CD276 Protein Expression (Median H Score)") + xlab("Objective Response Category")
sfig.resp.b

#### stop ####
#### save figures as PDF ####

pdf("figures/pptc_cd276_Figure1.pdf",width = 6.5,height = 3)
plot_grid(fig1a,fig1b, labels=c("A","B"), align = 'h', nrow = 1, rel_widths = c(1/1,1/3), label_size=12)
dev.off()

pdf("figures/pptc_cd276_Figure2.pdf",width = 6.5,height = 4.875)
ggdraw() +
  draw_plot(fig2c, 0,.67,1,.33) +
  draw_plot(fig2a,0,.27,1,.39) +
  draw_plot(fig2b, 0,0,.5,.28) +
  draw_plot_label(c("A","B","C"), c(0, 0, 0), c(1, 0.67, 0.3), size = 12)
dev.off()

pdf("figures/pptc_cd276_Figure3.pdf",width = 6.5,height = 5)
ggdraw() +
  draw_plot(fig3a,0,.5,.6,.5) +
  draw_plot(fig3b, 0,0,1,.5) +
  draw_plot_label(c("A","C"), c(0, 0), c(1, .5), size = 12)
dev.off()

pdf("figures/pptc_cd276_SFig1.pdf",width = 6,height = 6)
sfig1
dev.off()

pdf("figures/pptc_cd276_SFig9.pdf",width = 8,height = 5)
plot_grid(sfig.resp.a,sfig.resp.b, labels=c("A","B"), align = 'h', nrow = 1, rel_widths = c(.65,1), label_size=12)
dev.off()

#### prepare legend ####

mix.ihc.slim
nbl.ihc.slim
os.ihc.slim

total.ihc.slim <- rbind(mix.ihc.slim,rbind(nbl.ihc.slim,os.ihc.slim))

total.ihc.plot <- histo.colors %>%
  dplyr::filter(Histology %in% total.ihc.slim$Histology) %>%
  dplyr::pull(Color) %>%
  as.character()

total.ihc.ggplot <- ggplot(total.ihc.slim,aes(x=SAMPID,y=medianHScore,fill=Histology)) + 
  geom_boxplot(color="black",outlier.shape = NA) + 
  geom_point(size=1) +
  #scale_color_manual(values = c("black","darkgray")) +
  scale_fill_manual(values = total.ihc.plot) +
  #facet_grid(cols = vars(total.ihc.slim$Histology.Short), scales = "free", space = "free") +
  #ggtitle("Mixed Histology PDX TMA") +
  #xlab("") + ylab("CD276 H-Score") +
  theme_bw() + ylim(0,300) +
  theme(title = element_text(size = 12,colour = "black"),
        axis.text.x = element_text(size = 6, angle = 45, hjust=1, color="black"),
        axis.title.y = element_text(size=8),
        axis.title.x=element_blank(), 
        axis.text.y = element_text(size = 8,colour = "black"),
        strip.text = element_text(size=8))
total.ihc.ggplot

# pdf("~/Desktop/figure-for-legend.pdf",width=5,height=5)
# total.ihc.ggplot
# dev.off()
