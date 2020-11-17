#### delcare datasets / tumor type / output directors ####

home <- "/mnt/isilon/maris_lab/target_nbl_ngs/Kendsersky/crc-project/tumor_datasets/"
dataset.dir <- "treehouse_data/"
analysis.dir <- "solidtumor_analysis/"
dataset <- gsub("_data/","",dataset.dir)
tumor <- "multi-tumor"

log.out <- paste0(home,dataset.dir,analysis.dir,"logs_out/")
checkpoints.out <- paste0(home,dataset.dir,analysis.dir,"checkpoints_out/")

logFile1 <- paste0(log.out,"Step1_log.txt")
if (!file.exists(logFile1)) {
  file.create(logFile1)
}

cat(paste0(Sys.time(),": RUNNING STEP1: DATA IMPORT ======================================="),file=logFile1,append=T,sep="\n")
cat(paste0(Sys.time(),": Analyzing ",tumor," samples in ", dataset," dataset."),file=logFile1,append=T,sep="\n")

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

#### load tumor and GTEx datasets ####

cat(paste0(Sys.time(),": Loading RNAseq datasets..."),file=logFile1,append=T,sep="\n")

if (dataset == "treehouse") {
  
  # dataset includes treehouse tumor compendium counts tsv
  tumor.counts <- read_tsv("/mnt/isilon/maris_lab/target_nbl_ngs/Kendsersky/TumorCompendium_v11_PublicPolyA/TumorCompendium_v11_PolyA_ensembl_expected_count_58581genes_2020-04-09.tsv")
  tumor.annot <- read_tsv("/mnt/isilon/maris_lab/target_nbl_ngs/Kendsersky/TumorCompendium_v11_PublicPolyA/clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv")
  
  cat(paste0(Sys.time(),":   Successfully loaded tumor dataset."),file=logFile1,append=T,sep="\n")
  
}

# gtex database, counts
gtex.counts <- read_tsv("/mnt/isilon/maris_lab/target_nbl_ngs/Kendsersky/GTEx_Analysis_2017-06-05_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.txt")
gtex.annot <- read_tsv("/mnt/isilon/maris_lab/target_nbl_ngs/Kendsersky/GTEx_Analysis_2017-06-05_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
## remove name, and change "Description" to "Gene"
gtex.counts <- gtex.counts[,-2]
colnames(gtex.counts)[1] <- "Gene"

cat(paste0(Sys.time(),":   Successfully loaded gtex dataset."),file=logFile1,append=T,sep="\n")


#### combining tumor and gtex datasets ####

cat(paste0(Sys.time(),": Subseting and merging tumor/gtex count matrices..."),file=logFile1,append=T,sep="\n")


# subset tumor histologies in annotation file

if (dataset == "treehouse" & tumor == "multi-tumor") {
  ped.tumors <- c("Ewing sarcoma","rhabdoid tumor","rhabdomyosarcoma","alveolar rhabdomyosarcoma","embryonal rhabdomyosarcoma",
                  "hepatoblastoma","neuroblastoma","meningioma","osteosarcoma","wilms tumor")
  tumor.annot.sub <- tumor.annot[tumor.annot$disease %in% ped.tumors,]
  
  # subset tumor histology in counts file
  tumor.counts.sub <- tumor.counts[,c("Gene",tumor.annot.sub$th_sampleid)]
  
  tumor.annot.sub.slim <- tumor.annot.sub[,c("th_sampleid","disease")]
  colnames(tumor.annot.sub.slim) <- c("id","tissue")
  
}


# prepare annotation df rms and gtex; tumor.gtex.annot
gtex.annot.slim <- gtex.annot[,c("SAMPID","SMTS")]
colnames(gtex.annot.slim) <- c("id","tissue")
tumor.gtex.annot <- rbind(tumor.annot.sub.slim,gtex.annot.slim)
dim(tumor.gtex.annot)

# prepare counts matrix with tumor data
tumor.counts.sub$Gene <- gsub("\\..*","",tumor.counts.sub$Gene)

# prepare counts matrix for gtex... remove DUP genes WITH rowSums == 0
gtex.counts$Gene <- gsub("\\..*","",gtex.counts$Gene)
dup.genes <- gtex.counts$Gene[!duplicated(gtex.counts$Gene)]
gtex.counts <- gtex.counts[!(gtex.counts$Gene %in% dup.genes & rowSums(gtex.counts[,-1])==0),]

cat(paste0(Sys.time(),":   Identified ",length(unique(intersect(gtex.counts$Gene,tumor.counts.sub$Gene)))," common genes."),file=logFile1,append=T,sep="\n")

# merge into a matrix
tumor.gtex.mat <- matrix.please(merge(tumor.counts.sub,gtex.counts,by="Gene"))
sample.in.common <- unique(intersect(colnames(tumor.gtex.mat),tumor.gtex.annot$id))
tumor.gtex.mat <- tumor.gtex.mat[,sample.in.common]
tumor.gtex.annot <- tumor.gtex.annot %>%
  dplyr::filter(id %in% sample.in.common)

cat(paste0(Sys.time(),":   Dimensions of merged tumor/gtex matrix is ",nrow(tumor.gtex.mat)," x ",ncol(tumor.gtex.mat),"."),file=logFile1,append=T,sep="\n")
cat(paste0(Sys.time(),":   Number of unique genes in matrix is ",length(unique(rownames(tumor.gtex.mat))),"."),file=logFile1,append=T,sep="\n")

#### perform filtering and normalization ####

cat(paste0(Sys.time(),": Filtering and normalizing myDGEList..."),file=logFile1,append=T,sep="\n")

## counts, unfiltered, non-normalied
myDGEList <- DGEList(tumor.gtex.mat)
# use the 'cpm' function from EdgeR to get counts per million
cpm <- cpm(myDGEList)

cat(paste0(Sys.time(),":   Number of genes is ",nrow(myDGEList$counts),"."),file=logFile1,append=T,sep="\n")

## filter and normalize
min.group.num <- min(table(tumor.gtex.annot$tissue)) # smallest group is 9, Fallopian Tube
min.group.name <- names(which(table(tumor.gtex.annot$tissue) == min.group.num))

cat(paste0(Sys.time(),":   Smallest group is ",min.group.name," with ",min.group.num," samples."),file=logFile1,append=T,sep="\n")
keepers <- rowSums(cpm>=1)>=min.group.num

# to  memory, remove cpm after finding "keepers" and bring back myDGEList
myDGEList.filtered <- myDGEList[keepers,]
cat(paste0(Sys.time(),":   After filtering, number of genes is ",nrow(myDGEList.filtered$counts),"."),file=logFile1,append=T,sep="\n")

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

cat(paste0(Sys.time(),":   Successfully performed calcNormFactors() with TMM method."),file=logFile1,append=T,sep="\n")

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
dup.genes <- log2.cpm.filtered.norm.df$geneID[duplicated(log2.cpm.filtered.norm.df$geneID)]
log2.cpm.filtered.norm.df.dups <- log2.cpm.filtered.norm.df %>%
  dplyr::filter(geneID %in% dup.genes)

cat(paste0(Sys.time(),":   There are ",length(dup.genes)," duplicated genes in the DGEList."),file=logFile1,append=T,sep="\n")

#### save some checkpoints ####
save(myDGEList.filtered.norm, file = paste0(checkpoints.out,"myDGEList.filtered.norm.RData"))
write.table(log2.cpm.filtered.norm.df, file = paste0(checkpoints.out,"log2.cpm.filtered.norm.df"),quote=F,row.names = F,col.names = T,sep = "\t")

save.image(file=paste0(checkpoints.out,Sys.Date(),"_",dataset,"-",tumor,"-workspaceFromScript1.RData"))
cat(paste0(Sys.time(),": Saved myDGEList.filtered.norm and log2.cpm.filtered.norm.df in ",checkpoints.out,"."),file=logFile1,append=T,sep="\n")

mem <- gc(full=T)
mem.max.gb <- (mem[1,6]+mem[2,6])/1024
cat(paste0(Sys.time(),": Total memory usage (gc) was ",mem.max.gb," GB."),file=logFile1,append=T,sep="\n")