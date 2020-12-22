#########################################################################
### CLIP Benchmarking                                                 ###
### Identification of cancer driver genes                             ###
#########################################################################
source("clip-meta/R/clip/clip_functions.R")
source("clip-meta/R/reproducibility_analysis/reproducibility_analysis_functions.R")

# known cancer driver genes
# https://www.sciencedirect.com/science/article/pii/S009286741830237X#app2
gold.std.genes <- readxl::read_xlsx("clip-meta/data/1-s2.0-S009286741830237X-mmc1.xlsx", sheet=2, skip=2)

#subset to cancer driver genes detected Breast cancer and PanCancer dataset
gold.std.genes <- gold.std.genes[gold.std.genes$Cancer %in% c("BRCA", "PANCAN"), ]
gold.std.genes <- unique(gold.std.genes$Gene)

#Load rCCS data from CLIP analysis for breast cancer cell lines
load(file="clip-meta/processed_files/clip/CLIP_rCCS_Breast.RData")

# Breast cancer subtype info
brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)

# Driver genes identified by CLIP for each breast cancer subtype
brca.subtype <- brca.data.view
brca.subtype$ER_pos <- brca.subtype$subtype_three_receptor
brca.subtype$ER_pos[which(brca.subtype$ER_pos != "ER")]  = 0
brca.subtype$ER_pos[which(brca.subtype$ER_pos == "ER")]  = 1

brca.subtype$HER2_pos <- brca.subtype$subtype_three_receptor
brca.subtype$HER2_pos[which(brca.subtype$HER2_pos != "ERBB2")]  = 0
brca.subtype$HER2_pos[which(brca.subtype$HER2_pos == "ERBB2")]  = 1

brca.subtype$TNBC <- brca.subtype$subtype_three_receptor
brca.subtype$TNBC[which(brca.subtype$TNBC != "TNBC")]  = 0
brca.subtype$TNBC[which(brca.subtype$TNBC == "TNBC")]  = 1

brca.subtype <- brca.subtype[,c(2,8:10)]
brca.subtype[,-1] <- apply(brca.subtype[,-1],2,as.numeric)

# Statistically enriched driver genes identified by CLIP
clip.subtype.drivers <- list()
for (c in 2:ncol(brca.subtype)){
  subtype.compare.name <- colnames(brca.subtype)[c]
  subtype.pos.cells <- brca.subtype$CellName[which(brca.subtype[,c]==1)]
  subtype.neg.cells <- brca.subtype$CellName[brca.subtype[,c]==0]
  data <- brca.rCCSmat
  rownames(data) <- data$GeneID
  data <- data[,-1]
  t.csp.genes.FE <- PerformFisherExact(data, subtype.pos.cells, subtype.neg.cells, "greater")
  t.csp.genes.FE <- t.csp.genes.FE[order(names(t.csp.genes.FE))]
  t.csp.genes.FE <- cbind(GeneID=names(t.csp.genes.FE), Pval =t.csp.genes.FE)
  t.csp.genes.FE <- as.data.frame(t.csp.genes.FE, stringsAsFactors = F)
  clip.subtype.drivers[[subtype.compare.name]] <- t.csp.genes.FE
}
clip.subtype.drivers <- Reduce(function(...) merge(...,by = "GeneID", all=T), clip.subtype.drivers)
colnames(clip.subtype.drivers)[-1] <- colnames(brca.subtype)[-1]
clip.subtype.drivers <- as.data.frame(clip.subtype.drivers, stringsAsFactors =F)
clip.subtype.drivers[,-1] <- apply(clip.subtype.drivers[,-1],2,as.numeric)
View(clip.subtype.drivers)

# rCCS genes for each subtype
er_pos_rCCS <- na.omit(clip.subtype.drivers$GeneID[clip.subtype.drivers$ER_pos < 0.05])
her2_pos_rCCS <- na.omit(clip.subtype.drivers$GeneID[clip.subtype.drivers$HER2_pos < 0.05])
tnbc_pos_rCCS <- na.omit(clip.subtype.drivers$GeneID[clip.subtype.drivers$TNBC < 0.05])

# OVerlap with true positive drivers
clip.er.drivers <- intersect(er_pos_rCCS, gold.std.genes) 
clip.her2.drivers <- intersect(her2_pos_rCCS, gold.std.genes) 
clip.tnbc.drivers <- intersect(tnbc_pos_rCCS, gold.std.genes) 
clip.unique <- unique(c(clip.er.drivers, clip.her2.drivers, clip.tnbc.drivers))


### Benchmarking with Differential analysis #######

#Perform Differential Expression analysis on Transcriptome data
# load Gene expression data
path = "clip-meta/data/"
modality = "GEXP"
GEXP.list <- loadData(path, modality)
for (i in 1:length(GEXP.list)){
  dname <- names(GEXP.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = GEXP.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(GEXP.list)

#Log2 normalization of RNAseq studies
gexp.broad.log <-  apply(gexp.broad,2,function(x) log2(x+1))
gexp.gcsi.log <-  apply(gexp.gcsi,2,function(x) log2(x+1))
gexp.ohsu.log <-  apply(gexp.ohsu,2,function(x) log2(x+1))
gexp.list <- list(BROAD=gexp.broad.log, UHN= gexp.uhn, gCSI=gexp.gcsi.log, SANGER=gexp.sanger,OHSU=gexp.ohsu.log )
gexp.list <- lapply(gexp.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
gexp.list <- lapply(gexp.list,QuantNormScale)
driver.GEXP.DE <- PerformDE(brca.subtype, gexp.list, "GEXP",  ref.drivers = gold.std.genes )

#Perform DE analysis on Functional data
# load Functional screening data
path = "clip-meta/data/"
modality = "FUNC"
FUNC.list <- loadData(path, modality)
for (i in 1:length(FUNC.list)){
  dname <- names(FUNC.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = FUNC.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(FUNC.list)

func.list <- list(BROAD_ACHILLES=func.broad_achilles, BROAD_AVANA=func.broad_avana,
                  SANGER=func.sanger, UHN=func.uhn, DRIVE=func.drive )
func.list <- lapply(func.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
func.list <- lapply(func.list,QuantNormScale)
driver.FUNC.DE <- PerformDE(brca.subtype, func.list, "FUNC",  ref.drivers = gold.std.genes )

#Perform DE analysis in Methylation data
path = "clip-meta/data/"
modality = "METH"
METH.list <- loadData(path, modality)
for (i in 1:length(METH.list)){
  dname <- names(METH.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = METH.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(METH.list)

meth.list <- list(BROAD=meth.broad,SANGER=meth.sanger,OHSU=meth.ohsu )
meth.list <- lapply(meth.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
meth.list <- lapply(meth.list,QuantNormScale)
driver.METH.DE <- PerformDE(brca.subtype, meth.list, "METH",  ref.drivers = gold.std.genes )

#Perform DE analysis in PEXP data
# load Protein expression data
path = "clip-meta/data/"
modality = "PEXP"
PEXP.list <- loadData(path, modality)
for (i in 1:length(PEXP.list)){
  dname <- names(PEXP.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = PEXP.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(PEXP.list)

pexp.mghcc_breast.log <-  apply(pexp.mghcc_breast,2,log2)
pexp.list <- list(MGHCC=pexp.mghcc_breast.log, CCLE=pexp.ccle)
pexp.list <- lapply(pexp.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
pexp.list <- lapply(pexp.list,QuantNormScale)
driver.PEXP.DE <- PerformDE(brca.subtype, pexp.list, "PEXP",  ref.drivers = gold.std.genes )

#Perform DE analysis in CNV data
# load Copy number variation data
path = "clip-meta/data/"
modality = "CNV"
CNV.list <- loadData(path, modality)
for (i in 1:length(CNV.list)){
  dname <- names(CNV.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = CNV.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(CNV.list)

cnv.list <- list(BROAD=cnv.broad,SANGER=cnv.sanger,OHSU=cnv.ohsu, gCSI=cnv.gcsi, UHN=cnv.uhn )
cnv.list <- lapply(cnv.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
driver.CNV.DE <- PerformDE(brca.subtype, cnv.list, "CNV",  ref.drivers = gold.std.genes )

# Recall of unique true drivers
drivers.DE <- c(driver.GEXP.DE, driver.FUNC.DE, driver.METH.DE, driver.PEXP.DE, driver.CNV.DE)
x <- seq(1,15,3)
drivers.unique <- c()
for (i in x){
  k <- i + 2
  uniq.drivers <- unique(unlist(drivers.DE[i:k]))
  drivers.unique <- c(drivers.unique, length(uniq.drivers))
}
names(drivers.unique) <- c("GEXP", "FUNC", "METH", "PEXP", "CNV")
drivers.unique <- c(drivers.unique, CLIP=length(clip.unique))
drivers.unique <- c(drivers.unique, CLIP=27)
drivers.unique <- drivers.unique[c(6,3,5,1,4,2)]
drivers.unique.prop <- drivers.unique/201

barplot(drivers.unique.prop, beside = F, las=2,ylab="TP fraction", ylim=c(0,0.14),
        col=c("gray", "steelblue1", "darkgoldenrod3", 'royalblue4', "peachpuff", "salmon2"),
        border=c("gray", "steelblue1", "darkgoldenrod3", 'royalblue4', "peachpuff", "salmon2"))




### Benchmarking with MOFA+ on BROAD dataset #####

library(reticulate)
use_python("/usr/local/bin/python3.8", required = TRUE)

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

path = "clip-meta/data/"

# load Methylation profiles
modality = "METH"
METH.list <- loadData(path, modality)
for (i in 1:length(METH.list)){
  dname <- names(METH.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = METH.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(METH.list)

# load Mutation profiles
modality = "MUT"
MUT.list <- loadData(path, modality)
for (i in 1:length(MUT.list)){
  dname <- names(MUT.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = MUT.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(MUT.list)

# load Copy number variation data
modality = "CNV"
CNV.list <- loadData(path, modality)
for (i in 1:length(CNV.list)){
  dname <- names(CNV.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = CNV.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(CNV.list)

# load Gene expression data
modality = "GEXP"
GEXP.list <- loadData(path, modality)
for (i in 1:length(GEXP.list)){
  dname <- names(GEXP.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = GEXP.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(GEXP.list)

# load Protein expression data
modality = "PEXP"
PEXP.list <- loadData(path, modality)
for (i in 1:length(PEXP.list)){
  dname <- names(PEXP.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = PEXP.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(PEXP.list)

# laod Protein phohphorylation data
modality = "PHOS"
PHOS.list <- loadData(path, modality)
for (i in 1:length(PHOS.list)){
  dname <- names(PHOS.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = PHOS.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(PHOS.list)

# load Functional screening data
modality = "FUNC"
FUNC.list <- loadData(path, modality)
for (i in 1:length(FUNC.list)){
  dname <- names(FUNC.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = FUNC.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(FUNC.list)

# load Drug sensitivity data
modality = "DSS"
DSS.list <- loadData(path, modality)
for (i in 1:length(DSS.list)){
  dname <- names(DSS.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = DSS.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(DSS.list)

# load Target Addiction data
modality = "TAS"
TAS.list <- loadData(path, modality)
for (i in 1:length(TAS.list)){
  dname <- names(TAS.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = TAS.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(TAS.list)

#log normalize RNAseq profiles before MOFA+ run
gexp.broad.log <- apply(gexp.broad,2,function(x) log2(x+1))

#
broad.data <- list(CNV=as.matrix(cnv.broad),
                  GEXP=as.matrix(gexp.broad.log),
                  METH=as.matrix(meth.broad),
                  MUT=as.matrix(mut.broad),
                  PEXP=as.matrix(pexp.broad),
                  FUNC_CRISPR=as.matrix(func.broad_avana),
                  FUNC_RNAI=as.matrix(func.broad_achilles))
#PHOS and TAS modality not included because of low dimensions

#Subset each dataset to coding genes
homo.genes <- read.table("clip-meta/data/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

#Subset to Breast cell lines and coding genes
brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)
breast.cells <- brca.data.view$CellName[!is.na(brca.data.view$subtype_three_receptor)]
broad.data <- lapply(broad.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(broad.data, dim)

#Subset MUT and CNV data to variable genes
broad.data$MUT <- broad.data$MUT[rowSums(broad.data$MUT)!=0, ]
broad.data$CNV <- broad.data$CNV[complete.cases(broad.data$CNV), ]
broad.data$CNV <- broad.data$CNV[rowSums(broad.data$CNV)!=0, ]

#Quantile normalize continuous datasets
broad.data$GEXP <- QuantNormScale(broad.data$GEXP)
broad.data$METH <- QuantNormScale(broad.data$METH)
broad.data$FUNC_CRISPR <- QuantNormScale(broad.data$FUNC_CRISPR)
broad.data$FUNC_RNAI <- QuantNormScale(broad.data$FUNC_RNAI)
broad.data$PEXP <- QuantNormScale(broad.data$PEXP)

cell.freq <- table(unlist(sapply(broad.data, colnames)))
# n <- length(broad.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
broad.brca <- lapply(broad.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(broad.brca, dim)

#Create MOFA object
MOFAobject.broad <- create_mofa(broad.brca)
data_opts <- get_default_data_options(MOFAobject.broad)
model_opts <- get_default_model_options(MOFAobject.broad)

#Add subtype information
sample_metadata <- MOFAobject.broad@samples_metadata
sample_metadata <- merge(sample_metadata, brca.data.view[,c(2,6)], by.x="sample", by.y="CellName")
sample_metadata$subtype_three_receptor[is.na(sample_metadata$subtype_three_receptor)] <- "Unclassified"
sample_metadata <- sample_metadata[,-2]
colnames(sample_metadata)[2] <- "group"
sample_metadata$group <- as.factor(sample_metadata$group)

#Run MOFA with default parameters
MOFAobject.broad <- create_mofa(broad.brca, groups=sample_metadata$group)
data_opts <- get_default_data_options(MOFAobject.broad)
model_opts <- get_default_model_options(MOFAobject.broad)
train_opts <- get_default_training_options(MOFAobject.broad)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.broad <- prepare_mofa(MOFAobject.broad,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.broad <- run_mofa(MOFAobject.broad)

#Total Variance explained per modality
head(MOFAobject.broad@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.broad, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.broad, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.broad, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.broad, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

#factor correlation
plot_factor_cor(MOFAobject.broad)

#Top Factors for each subgroup
#ER+: Factor1,2,6
#ERBB2+: Factor1,3,4,6
#TNBC: Factor3,4,5,7,8

#Plot factor scores per group
plot_factor(MOFAobject.broad, 
            factor = c(1:15),
            color_by = "group")

#Plot top feature weight scores per modality
plot_weights(MOFAobject.broad,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.broad.factors <- get_factors(MOFAobject.broad, as.data.frame = T)
mofa.broad.weights <- get_weights(MOFAobject.broad, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.broad.top.features <- ExtractTopFeatures(mofa.broad.weights)

#Overlap with known driver genes
mofa.broad.top.features.gold <- mofa.broad.top.features[mofa.broad.top.features$FeatureID %in% gold.std.genes,  ]

#subset to top factors and views
factors.vec <- as.character(unique(mofa.broad.top.features.gold$factor))
subset_combos <- rbind(c('Factor1', "CNV" ),
                       c('Factor1', "GEXP" ),
                       c('Factor1', "METH" ),
                       c('Factor1', "PEXP" ),
                       c('Factor2', "CNV" ),
                       c('Factor2', "GEXP" ),
                       c('Factor2', "METH" ),
                       c('Factor2', "FUNC_RNAI" ),
                       c('Factor3', "FUNC_CRISPR" ),
                       c('Factor3', "FUNC_RNAI" ),
                       c('Factor4', "FUNC_CRISPR" ),
                       c('Factor5', "FUNC_CRISPR" ),
                       c('Factor5', "FUNC_RNAI" ),
                       c('Factor6', "FUNC_CRISPR" ),
                       c('Factor7', "CNV" ),
                       c('Factor7', "GEXP" ),
                       c('Factor7', "METH" ),
                       c('Factor8', "GEXP" ))

mofa.broad.top.features.gold.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.broad.top.features.gold.sub <- mofa.broad.top.features.gold[mofa.broad.top.features.gold$factor == c[1] & mofa.broad.top.features.gold$view == c[2], ]
  mofa.broad.top.features.gold.ev <- rbind(mofa.broad.top.features.gold.ev, mofa.broad.top.features.gold.sub)
}

### Benchmarking with MOFA+ on SANGER dataset #####

#
sanger.data <- list(CNV=as.matrix(cnv.sanger),
                    GEXP=as.matrix(gexp.sanger),
                    METH=as.matrix(meth.sanger),
                    MUT=as.matrix(mut.sanger),
                    FUNC=as.matrix(func.sanger))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
sanger.data <- lapply(sanger.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(sanger.data, dim)

#Subset MUT and CNV data to variable genes
sanger.data$MUT <- sanger.data$MUT[rowSums(sanger.data$MUT)!=0, ]
sanger.data$CNV <- sanger.data$CNV[complete.cases(sanger.data$CNV), ]
sanger.data$CNV <- sanger.data$CNV[rowSums(sanger.data$CNV)!=0, ]

#Quantile normalize continuous datasets
sanger.data$GEXP <- QuantNormScale(sanger.data$GEXP)
sanger.data$METH <- QuantNormScale(sanger.data$METH)
sanger.data$FUNC <- QuantNormScale(sanger.data$FUNC)

cell.freq <- table(unlist(sapply(sanger.data, colnames)))
# n <- length(sanger.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
sanger.brca <- lapply(sanger.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(sanger.brca, dim)

#Create MOFA object
MOFAobject.sanger <- create_mofa(sanger.brca)
data_opts <- get_default_data_options(MOFAobject.sanger)
model_opts <- get_default_model_options(MOFAobject.sanger)

#Add subtype information
sample_metadata <- MOFAobject.sanger@samples_metadata
sample_metadata <- merge(sample_metadata, brca.data.view[,c(2,6)], by.x="sample", by.y="CellName")
sample_metadata$subtype_three_receptor[is.na(sample_metadata$subtype_three_receptor)] <- "Unclassified"
sample_metadata <- sample_metadata[,-2]
colnames(sample_metadata)[2] <- "group"
sample_metadata$group <- as.factor(sample_metadata$group)

#Run MOFA with default parameters
MOFAobject.sanger <- create_mofa(sanger.brca, groups=sample_metadata$group)
data_opts <- get_default_data_options(MOFAobject.sanger)
model_opts <- get_default_model_options(MOFAobject.sanger)
train_opts <- get_default_training_options(MOFAobject.sanger)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.sanger <- prepare_mofa(MOFAobject.sanger,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.sanger <- run_mofa(MOFAobject.sanger)

#Total Variance explained per modality
head(MOFAobject.sanger@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.sanger, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.sanger, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.sanger, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.sanger, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


#factor correlation
plot_factor_cor(MOFAobject.sanger)

#Plot factor scores per group
plot_factor(MOFAobject.sanger, 
            factor = c(1:15),
            color_by = "group")


#Plot top feature weight scores per modality
plot_weights(MOFAobject.sanger,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.sanger.factors <- get_factors(MOFAobject.sanger, as.data.frame = T)
mofa.sanger.weights <- get_weights(MOFAobject.sanger, as.data.frame = T)


#Top feature weights for each factor and overlap with driver genes
mofa.sanger.top.features <- ExtractTopFeatures(mofa.sanger.weights)

#Overlap with known driver genes
mofa.sanger.top.features.gold <- mofa.sanger.top.features[mofa.sanger.top.features$FeatureID %in% gold.std.genes,  ]

#subset to top factors and views
factors.vec <- as.character(unique(mofa.sanger.top.features.gold$factor))
subset_combos <- rbind(c('Factor1', "METH" ),
                       c('Factor2', "CNV" ),
                       c('Factor2', "FUNC" ),
                       c('Factor3', "CNV" ),
                       c('Factor3', "FUNC" ),
                       c('Factor4', "GEXP" ),
                       c('Factor4', "METH" ),
                       c('Factor5', "CNV" ),
                       c('Factor5', "GEXP" ),
                       c('Factor5', "METH" ),
                       c('Factor6', "CNV" ),
                       c('Factor6', "METH" ),
                       c('Factor6', "FUNC" ),
                       c('Factor7', "CNV" ),
                       c('Factor8', "CNV" ),
                       c('Factor9', "CNV" ),
                       c('Factor10', "CNV" ),
                       c('Factor10', "GEXP" ))

mofa.sanger.top.features.gold.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.sanger.top.features.gold.sub <- mofa.sanger.top.features.gold.ev[mofa.sanger.top.features.gold.ev$factor == c[1] & mofa.sanger.top.features.gold$view == c[2], ]
  mofa.sanger.top.features.gold.ev <- rbind(mofa.sanger.top.features.gold.sub, mofa.sanger.top.features.gold.ev)
}


### Benchmarking with MOFA+ on OHSU dataset #####

#
ohsu.data <- list(CNV=as.matrix(cnv.ohsu),
                  GEXP=as.matrix(gexp.ohsu),
                  METH=as.matrix(meth.ohsu),
                  MUT=as.matrix(mut.ohsu))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
ohsu.data <- lapply(ohsu.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(ohsu.data, dim)

#Subset MUT and CNV data to variable genes
ohsu.data$MUT <- ohsu.data$MUT[rowSums(ohsu.data$MUT)!=0, ]
ohsu.data$CNV <- ohsu.data$CNV[complete.cases(ohsu.data$CNV), ]
ohsu.data$CNV <- ohsu.data$CNV[rowSums(ohsu.data$CNV)!=0, ]

#Quantile normalize continuous datasets
ohsu.data$GEXP <- QuantNormScale(ohsu.data$GEXP)
ohsu.data$METH <- QuantNormScale(ohsu.data$METH)

cell.freq <- table(unlist(sapply(ohsu.data, colnames)))
# n <- length(ohsu.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
ohsu.brca <- lapply(ohsu.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(ohsu.brca, dim)

#Create MOFA object
MOFAobject.ohsu <- create_mofa(ohsu.brca)
data_opts <- get_default_data_options(MOFAobject.ohsu)
model_opts <- get_default_model_options(MOFAobject.ohsu)

#Add subtype information
sample_metadata <- MOFAobject.ohsu@samples_metadata
sample_metadata <- merge(sample_metadata, brca.data.view[,c(2,6)], by.x="sample", by.y="CellName")
sample_metadata$subtype_three_receptor[is.na(sample_metadata$subtype_three_receptor)] <- "Unclassified"
sample_metadata <- sample_metadata[,-2]
colnames(sample_metadata)[2] <- "group"
sample_metadata$group <- as.factor(sample_metadata$group)

#Run MOFA with default parameters
MOFAobject.ohsu <- create_mofa(ohsu.brca, groups=sample_metadata$group)
data_opts <- get_default_data_options(MOFAobject.ohsu)
model_opts <- get_default_model_options(MOFAobject.ohsu)
train_opts <- get_default_training_options(MOFAobject.ohsu)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ohsu <- prepare_mofa(MOFAobject.ohsu,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.ohsu <- run_mofa(MOFAobject.ohsu)

#Total Variance explained per modality
head(MOFAobject.ohsu@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.ohsu, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.ohsu, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.ohsu, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.ohsu, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


#Plot factor scores per group
plot_factor(MOFAobject.ohsu, 
            factor = c(1:15),
            color_by = "group")


#Plot top feature weight scores per modality
plot_weights(MOFAobject.ohsu,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.ohsu.factors <- get_factors(MOFAobject.ohsu, as.data.frame = T)
mofa.ohsu.weights <- get_weights(MOFAobject.ohsu, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.ohsu.top.features <- ExtractTopFeatures(mofa.ohsu.weights)

#Overlap with known driver genes
mofa.ohsu.top.features.gold <- mofa.ohsu.top.features[mofa.ohsu.top.features$FeatureID %in% gold.std.genes,  ]


#subset to top factors and views
factors.vec <- as.character(unique(mofa.ohsu.top.features.gold$factor))
subset_combos <- rbind(c('Factor1', "CNV" ),
                       c('Factor1', "GEXP" ),
                       c('Factor2', "CNV" ),
                       c('Factor2', "GEXP" ),
                       c('Factor2', "METH" ),
                       c('Factor3', "CNV" ),
                       c('Factor3', "GEXP" ),
                       c('Factor3', "METH" ),
                       c('Factor4', "CNV" ),
                       c('Factor4', "GEXP" ),
                       c('Factor4', "METH" ),
                       c('Factor5', "CNV" ),
                       c('Factor5', "GEXP" ),
                       c('Factor5', "METH" ),
                       c('Factor6', "CNV" ),
                       c('Factor6', "GEXP" ),
                       c('Factor6', "METH" ),
                       c('Factor7', "CNV" ),
                       c('Factor7', "GEXP" ),
                       c('Factor8', "CNV" ),
                       c('Factor8', "GEXP" ),
                       c('Factor8', "METH" ),
                       c('Factor9', "CNV" ),
                       c('Factor10', "CNV" ))

mofa.ohsu.top.features.gold.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.ohsu.top.features.gold.sub <- mofa.ohsu.top.features.gold.ev[mofa.ohsu.top.features.gold.ev$factor == c[1] & mofa.gdsc.top.features.gold$view == c[2], ]
  mofa.ohsu.top.features.gold.ev <- rbind(mofa.ohsu.top.features.gold.sub, mofa.ohsu.top.features.gold.ev)
}

### Benchmarking with MOFA+ on UHN dataset #####

#
uhn.data <- list(CNV=as.matrix(cnv.uhn),
                 GEXP=as.matrix(gexp.uhn),
                 FUNC=as.matrix(func.uhn))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
uhn.data <- lapply(uhn.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(uhn.data, dim)

#Subset MUT and CNV data to variable genes
uhn.data$CNV <- uhn.data$CNV[complete.cases(uhn.data$CNV), ]
uhn.data$CNV <- uhn.data$CNV[rowSums(uhn.data$CNV)!=0, ]

#Quantile normalize continuous datasets
uhn.data$GEXP <- QuantNormScale(uhn.data$GEXP)
uhn.data$FUNC <- QuantNormScale(uhn.data$FUNC)

cell.freq <- table(unlist(sapply(uhn.data, colnames)))
# n <- length(uhn.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
uhn.brca <- lapply(uhn.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(uhn.brca, dim)

#Create MOFA object
MOFAobject.uhn <- create_mofa(uhn.brca)
data_opts <- get_default_data_options(MOFAobject.uhn)
model_opts <- get_default_model_options(MOFAobject.uhn)

#Add subtype information
sample_metadata <- MOFAobject.uhn@samples_metadata
sample_metadata <- merge(sample_metadata, brca.data.view[,c(2,6)], by.x="sample", by.y="CellName")
sample_metadata$subtype_three_receptor[is.na(sample_metadata$subtype_three_receptor)] <- "Unclassified"
sample_metadata <- sample_metadata[,-2]
colnames(sample_metadata)[2] <- "group"
sample_metadata$group <- as.factor(sample_metadata$group)

#Run MOFA with default parameters
MOFAobject.uhn <- create_mofa(uhn.brca, groups=sample_metadata$group)
data_opts <- get_default_data_options(MOFAobject.uhn)
model_opts <- get_default_model_options(MOFAobject.uhn)
train_opts <- get_default_training_options(MOFAobject.uhn)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.uhn <- prepare_mofa(MOFAobject.uhn,
                               data_options = data_opts,
                               model_options = model_opts,
                               training_options = train_opts)
MOFAobject.uhn <- run_mofa(MOFAobject.uhn)


#Total Variance explained per modality
head(MOFAobject.uhn@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.uhn, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.uhn, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.uhn, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.uhn, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


#Plot factor scores per group
plot_factor(MOFAobject.uhn, 
            factor = c(1:15),
            color_by = "group")

#Plot top feature weight scores per modality
plot_weights(MOFAobject.uhn,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.uhn.factors <- get_factors(MOFAobject.uhn, as.data.frame = T)
mofa.uhn.weights <- get_weights(MOFAobject.uhn, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.uhn.top.features <- ExtractTopFeatures(mofa.uhn.weights)

#Overlap with known driver genes
mofa.uhn.top.features.gold <- mofa.uhn.top.features[mofa.uhn.top.features$FeatureID %in% gold.std.genes,  ]

#subset to top factors and views
factors.vec <- as.character(unique(mofa.uhn.top.features.gold$factor))
subset_combos <- rbind(c('Factor1', "METH" ),
                       c('Factor2', "CNV" ),
                       c('Factor2', "FUNC" ),
                       c('Factor3', "CNV" ),
                       c('Factor3', "FUNC" ),
                       c('Factor4', "GEXP" ),
                       c('Factor4', "METH" ),
                       c('Factor5', "CNV" ),
                       c('Factor5', "GEXP" ),
                       c('Factor5', "METH" ),
                       c('Factor6', "CNV" ),
                       c('Factor6', "METH" ),
                       c('Factor6', "FUNC" ),
                       c('Factor7', "CNV" ),
                       c('Factor8', "CNV" ),
                       c('Factor9', "CNV" ),
                       c('Factor10', "CNV" ),
                       c('Factor10', "GEXP" ))

mofa.uhn.top.features.gold.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.uhn.top.features.gold.sub <- mofa.uhn.top.features.gold.ev[mofa.uhn.top.features.gold.ev$factor == c[1] & mofa.gdsc.top.features.gold$view == c[2], ]
  mofa.uhn.top.features.gold.ev <- rbind(mofa.uhn.top.features.gold.sub, mofa.uhn.top.features.gold.ev)
}
mofa.uhn.top.features.gold <- mofa.uhn.top.features.gold[mofa.uhn.top.features.gold$factor %in% c(), ]



save(MOFAobject.uhn, MOFAobject.ccle,
     MOFAobject.ohsu, MOFAobject.gdsc,
     file="clip-meta/data//MOFA_BRCA_Subtype.RData", compress = "bzip2")




#########################################################################
### ECHDC1 identificaton                                              ###
#########################################################################
##### Using a simpler approach #####
# load Methylation profiles
source("clip-meta/R/reproducibility_analysis/reproducibility_analysis_functions.R")
path = "clip-meta/data/"
modality = "METH"
METH.list <- loadData(path, modality)
for (i in 1:length(METH.list)){
  dname <- names(METH.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = METH.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(METH.list)


brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)
breast.cells <- brca.data.view$CellName

#Subset to coding genes and breast cell lines
meth.broad <- meth.broad[intersect(homo.genes.coding, rownames(meth.broad)),  
                         sort(intersect(breast.cells,colnames(meth.broad)))]
meth.sanger <- meth.sanger[intersect(homo.genes.coding, rownames(meth.sanger)),
                       sort(intersect(breast.cells,colnames(meth.sanger)))]
meth.ohsu <- meth.ohsu[intersect(homo.genes.coding, rownames(meth.ohsu)),
                       sort(intersect(breast.cells,colnames(meth.ohsu)))]

# Quantile normalize each dataset
meth.broad.qnorm <- QuantNormScale(meth.broad)
meth.sanger.qnorm <- QuantNormScale(meth.sanger)
meth.ohsu.qnorm <- QuantNormScale(meth.ohsu)

#Z-scaling
# Broad
meth.brca.broad.zscaled <- apply(t(meth.broad.qnorm),2,function(e) scale(e, center = T))
meth.brca.broad.zscaled <- t(meth.brca.broad.zscaled)
rownames(meth.brca.broad.zscaled) <- rownames(meth.broad.qnorm)
colnames(meth.brca.broad.zscaled) <- colnames(meth.broad.qnorm)

# SANGER
meth.brca.sanger.zscaled <- apply(t(meth.sanger.qnorm),2,function(e) scale(e, center = T))
meth.brca.sanger.zscaled <- t(meth.brca.sanger.zscaled)
rownames(meth.brca.sanger.zscaled) <- rownames(meth.sanger.qnorm)
colnames(meth.brca.sanger.zscaled) <- colnames(meth.sanger.qnorm)

#OHSU
meth.brca.ohsu.zscaled <- apply(t(meth.ohsu.qnorm),2,function(e) scale(e, center = T))
meth.brca.ohsu.zscaled <- t(meth.brca.ohsu.zscaled)
rownames(meth.brca.ohsu.zscaled) <- rownames(meth.ohsu.qnorm)
colnames(meth.brca.ohsu.zscaled) <- colnames(meth.ohsu.qnorm)

# Count outlier frequency
meth.brca.broad.zscaled.long <- melt(as.matrix(meth.brca.broad.zscaled))
meth.brca.broad.zscaled.long$Keep = "F"
meth.brca.broad.zscaled.long$Keep[abs(meth.brca.broad.zscaled.long$value) >=sd(meth.brca.broad.zscaled.long$value, na.rm = T)*1.66] = "T"
meth.brca.broad.zscaled.long <- meth.brca.broad.zscaled.long[meth.brca.broad.zscaled.long$Var2 %in% breast.cells, ]
ccs.freq <- as.data.frame(table(meth.brca.broad.zscaled.long$Var1,meth.brca.broad.zscaled.long$Keep ))
ccs.freq <- ccs.freq[ccs.freq$Var2 != "F",]
ccs.freq[ccs.freq$Var1 == "ECHDC1", ]

# SANGER
meth.brca.sanger.zscaled.long <- melt(as.matrix(meth.brca.sanger.zscaled))
meth.brca.sanger.zscaled.long$Keep = "F"
meth.brca.sanger.zscaled.long$Keep[abs(meth.brca.sanger.zscaled.long$value) >=sd(meth.brca.sanger.zscaled.long$value, na.rm = T)*1.66] = "T"
meth.brca.sanger.zscaled.long <- meth.brca.sanger.zscaled.long[meth.brca.sanger.zscaled.long$Var2 %in% breast.cells, ]
ccs.freq <- as.data.frame(table(meth.brca.sanger.zscaled.long$Var1,meth.brca.sanger.zscaled.long$Keep ))
ccs.freq <- ccs.freq[ccs.freq$Var2 != "F",]
ccs.freq[ccs.freq$Var1 == "ECHDC1", ]

# OHSU
meth.brca.ohsu.zscaled.long <- melt(as.matrix(meth.brca.ohsu.zscaled))
meth.brca.ohsu.zscaled.long$Keep = "F"
meth.brca.ohsu.zscaled.long$Keep[abs(meth.brca.ohsu.zscaled.long$value) >=sd(meth.brca.ohsu.zscaled.long$value, na.rm = T)*1.66] = "T"
meth.brca.ohsu.zscaled.long <- meth.brca.ohsu.zscaled.long[meth.brca.ohsu.zscaled.long$Var2 %in% breast.cells, ]
ccs.freq <- as.data.frame(table(meth.brca.ohsu.zscaled.long$Var1,meth.brca.ohsu.zscaled.long$Keep ))
ccs.freq <- ccs.freq[ccs.freq$Var2 != "F",]
ccs.freq[ccs.freq$Var1 == "ECHDC1", ]

freq.ccs <- c(24,6,3)
names(freq.ccs) <- c("CLIP", "CCLE", "SANGER")
barplot(freq.ccs, las=2, ylab= "CCS frequency", border = NA, ylim = c(0,25))

##### MOFA+ BROAD #####

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)


source("clip-meta/R/reproducibility_analysis/reproducibility_analysis_functions.R")

path = "clip-meta/data/"

# load Methylation profiles
modality = "METH"
METH.list <- loadData(path, modality)
for (i in 1:length(METH.list)){
  dname <- names(METH.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = METH.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(METH.list)

# load Mutation profiles
modality = "MUT"
MUT.list <- loadData(path, modality)
for (i in 1:length(MUT.list)){
  dname <- names(MUT.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = MUT.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(MUT.list)

# load Copy number variation data
modality = "CNV"
CNV.list <- loadData(path, modality)
for (i in 1:length(CNV.list)){
  dname <- names(CNV.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = CNV.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(CNV.list)

# load Gene expression data
modality = "GEXP"
GEXP.list <- loadData(path, modality)
for (i in 1:length(GEXP.list)){
  dname <- names(GEXP.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = GEXP.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(GEXP.list)

# load Protein expression data
modality = "PEXP"
PEXP.list <- loadData(path, modality)
for (i in 1:length(PEXP.list)){
  dname <- names(PEXP.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = PEXP.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(PEXP.list)

# laod Protein phohphorylation data
modality = "PHOS"
PHOS.list <- loadData(path, modality)
for (i in 1:length(PHOS.list)){
  dname <- names(PHOS.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = PHOS.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(PHOS.list)

# load Functional screening data
modality = "FUNC"
FUNC.list <- loadData(path, modality)
for (i in 1:length(FUNC.list)){
  dname <- names(FUNC.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = FUNC.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(FUNC.list)

# load Drug sensitivity data
modality = "DSS"
DSS.list <- loadData(path, modality)
for (i in 1:length(DSS.list)){
  dname <- names(DSS.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = DSS.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(DSS.list)

# load Target Addiction data
modality = "TAS"
TAS.list <- loadData(path, modality)
for (i in 1:length(TAS.list)){
  dname <- names(TAS.list)[i]
  dname <- tolower(paste(modality, dname, sep="."))
  mat = TAS.list[[i]]
  assign(dname, as.data.frame(mat, stringsAsFactors = F))
  rm(mat)
}
rm(TAS.list)

#
homo.genes <- read.table("clip-meta/data/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

# log normalize RNASeq RPKM data before MOFA+ run
gexp.broad.log <- apply(gexp.broad,2,function(x) log2(x+1))

#
broad.data <- list(CNV=as.matrix(cnv.broad),
                  GEXP=as.matrix(gexp.broad.log),
                  METH=as.matrix(meth.broad),
                  MUT=as.matrix(mut.broad),
                  PEXP=as.matrix(pexp.broad),
                  FUNC_CRISPR=as.matrix(func.broad_avana),
                  FUNC_RNAI=as.matrix(func.broad_achilles))

#Subset to Breast cell lines and coding genes
brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)
breast.cells <- brca.data.view$CellName[!is.na(brca.data.view$subtype_three_receptor)]
broad.data <- lapply(broad.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(broad.data, dim)

#Check density distribution
par(mfrow=c(2,4))
pDensity <- lapply(broad.data, function(x){
  plot(density(as.matrix(x), na.rm = T))
})

#Subset MUT and CNV data to variable genes
broad.data$MUT <- broad.data$MUT[rowSums(broad.data$MUT)!=0, ]
broad.data$CNV <- broad.data$CNV[complete.cases(broad.data$CNV), ]
broad.data$CNV <- broad.data$CNV[rowSums(broad.data$CNV)!=0, ]

#Quantile normalize continuous datasets
broad.data$GEXP <- QuantNormScale(broad.data$GEXP)
broad.data$METH <- QuantNormScale(broad.data$METH)
broad.data$FUNC_CRISPR <- QuantNormScale(broad.data$FUNC_CRISPR)
broad.data$FUNC_RNAI <- QuantNormScale(broad.data$FUNC_RNAI)
broad.data$PEXP <- QuantNormScale(broad.data$PEXP)

cell.freq <- table(unlist(sapply(broad.data, colnames)))
# n <- length(broad.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
broad.brca <- lapply(broad.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(broad.brca, dim)

#Create MOFA object
MOFAobject.broad <- create_mofa(broad.brca)
data_opts <- get_default_data_options(MOFAobject.broad)
model_opts <- get_default_model_options(MOFAobject.broad)

#Run MOFA with default parameters
train_opts <- get_default_training_options(MOFAobject.broad)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.broad <- prepare_mofa(MOFAobject.broad,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.broad <- run_mofa(MOFAobject.broad)


#Total Variance explained per modality
head(MOFAobject.broad@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.broad, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.broad, x="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.broad, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.broad, x="view")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

#factor correlation
plot_factor_cor(MOFAobject.broad)

#Top Factors for each subgroup
#ER+: Factor1,2,6
#ERBB2+: Factor1,3,4,6
#TNBC: Factor3,4,5,7,8

#Plot factor scores per group
plot_factor(MOFAobject.broad, 
            factor = c(1:15))

#Plot top feature weight scores per modality
plot_weights(MOFAobject.broad,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.broad.factors <- get_factors(MOFAobject.broad, as.data.frame = T)
mofa.broad.weights <- get_weights(MOFAobject.broad, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.broad.top.features <- ExtractTopFeatures(mofa.broad.weights)

#subset to top factors and views
factors.vec <- as.character(unique(mofa.broad.top.features$factor))
subset_combos <- rbind(c('Factor1', "GEXP" ),
                       c('Factor1', "PEXP" ),
                       c('Factor1', "METH" ),
                       c('Factor2', "GEXP" ),
                       c('Factor2', "PEXP" ),
                       c('Factor2', "METH" ),
                       c('Factor3', "GEXP" ),
                       c('Factor3', "PEXP" ),
                       c('Factor3', "METH" ),
                       c('Factor5', "GEXP" ),
                       c('Factor5', "PEXP" ),
                       c('Factor5', "METH" ),
                       c('Factor6', "GEXP" ),
                       c('Factor6', "PEXP" ),
                       c('Factor6', "METH" ),
                       c('Factor7', "GEXP" ),
                       c('Factor7', "PEXP" ),
                       c('Factor7', "METH" ),
                       c('Factor8', "GEXP" ),
                       c('Factor8', "PEXP" ),
                       c('Factor8', "METH" ))

mofa.broad.top.features.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.broad.top.features.sub <- mofa.broad.top.features[mofa.broad.top.features$factor == c[1] & mofa.broad.top.features$view == c[2], ]
  mofa.broad.top.features.ev <- rbind(mofa.broad.top.features.ev, mofa.broad.top.features.sub)
}

##### MOFA+ SANGER  #####

#
sanger.data <- list(CNV=as.matrix(cnv.sanger),
                  GEXP=as.matrix(gexp.sanger),
                  METH=as.matrix(meth.sanger),
                  MUT=as.matrix(mut.sanger),
                  FUNC=as.matrix(func.sanger))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
sanger.data <- lapply(sanger.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(sanger.data, dim)

#Subset MUT and CNV data to variable genes
sanger.data$MUT <- sanger.data$MUT[rowSums(sanger.data$MUT)!=0, ]
sanger.data$CNV <- sanger.data$CNV[complete.cases(sanger.data$CNV), ]
sanger.data$CNV <- sanger.data$CNV[rowSums(sanger.data$CNV)!=0, ]

#Quantile normalize continuous datasets
sanger.data$GEXP <- QuantNormScale(sanger.data$GEXP)
sanger.data$METH <- QuantNormScale(sanger.data$METH)
sanger.data$FUNC <- QuantNormScale(sanger.data$FUNC)

cell.freq <- table(unlist(sapply(sanger.data, colnames)))
# n <- length(sanger.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
sanger.brca <- lapply(sanger.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(sanger.brca, dim)

#Create MOFA object
MOFAobject.sanger <- create_mofa(sanger.brca)
data_opts <- get_default_data_options(MOFAobject.sanger)
model_opts <- get_default_model_options(MOFAobject.sanger)
train_opts <- get_default_training_options(MOFAobject.sanger)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.sanger <- prepare_mofa(MOFAobject.sanger,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.sanger <- run_mofa(MOFAobject.sanger)

#Total Variance explained per modality
head(MOFAobject.sanger@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.sanger, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.sanger, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.sanger, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.sanger, x="view")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
()

#factor correlation
plot_factor_cor(MOFAobject.sanger)

#Plot factor scores per group
plot_factor(MOFAobject.sanger, 
            factor = c(1:15))

#Plot top feature weight scores per modality
plot_weights(MOFAobject.sanger,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.sanger.factors <- get_factors(MOFAobject.sanger, as.data.frame = T)
mofa.sanger.weights <- get_weights(MOFAobject.sanger, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.sanger.top.features <- ExtractTopFeatures(mofa.sanger.weights)

#subset to top factors and views
factors.vec <- as.character(unique(mofa.sanger.top.features$factor))
subset_combos <- rbind(c('Factor1', "GEXP" ),
                       c('Factor1', "METH" ),
                       c('Factor2', "GEXP" ),
                       c('Factor2', "METH" ),
                       c('Factor3', "GEXP" ),
                       c('Factor3', "METH" ),
                       c('Factor4', "GEXP" ),
                       c('Factor4', "METH" ),
                       c('Factor6', "GEXP" ),
                       c('Factor6', "METH" ))
mofa.sanger.top.features.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.sanger.top.features.sub <- mofa.sanger.top.features[mofa.sanger.top.features$factor == c[1] & mofa.sanger.top.features$view == c[2], ]
  mofa.sanger.top.features.ev <- rbind(mofa.sanger.top.features.sub, mofa.sanger.top.features.ev)
}


##### MOFA+ OHSU  #####

#
ohsu.data <- list(CNV=as.matrix(cnv.ohsu),
                  GEXP=as.matrix(gexp.ohsu),
                  METH=as.matrix(meth.ohsu),
                  MUT=as.matrix(mut.ohsu))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
ohsu.data <- lapply(ohsu.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(ohsu.data, dim)

#Subset MUT and CNV data to variable genes
ohsu.data$MUT <- ohsu.data$MUT[rowSums(ohsu.data$MUT)!=0, ]
ohsu.data$CNV <- ohsu.data$CNV[complete.cases(ohsu.data$CNV), ]
ohsu.data$CNV <- ohsu.data$CNV[rowSums(ohsu.data$CNV)!=0, ]

#Quantile normalize continuous datasets
ohsu.data$GEXP <- QuantNormScale(ohsu.data$GEXP)
ohsu.data$METH <- QuantNormScale(ohsu.data$METH)

cell.freq <- table(unlist(sapply(ohsu.data, colnames)))
# n <- length(ohsu.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
ohsu.brca <- lapply(ohsu.data, function(x){
  overlap.cells <- intersect(names(cell.freq), colnames(x))
  y <- x[, overlap.cells]
  nonoverlap.cells <- setdiff(names(cell.freq), overlap.cells)
  y1 <- matrix(NA, nrow(y), length(nonoverlap.cells)) 
  rownames(y1) <- rownames(y)
  colnames(y1) <- nonoverlap.cells
  yFinal <- cbind(y, y1)
  yFinal <- yFinal[, order(colnames(yFinal))]
  return(yFinal)
})
lapply(ohsu.brca, dim)

#Create MOFA object
MOFAobject.ohsu <- create_mofa(ohsu.brca)
data_opts <- get_default_data_options(MOFAobject.ohsu)
model_opts <- get_default_model_options(MOFAobject.ohsu)
train_opts <- get_default_training_options(MOFAobject.ohsu)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ohsu <- prepare_mofa(MOFAobject.ohsu,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.ohsu <- run_mofa(MOFAobject.ohsu)




#Total Variance explained per modality
head(MOFAobject.ohsu@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.ohsu, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.ohsu, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.ohsu, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.ohsu, x="view")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

#Plot factor scores per group
plot_factor(MOFAobject.ohsu, 
            factor = c(1:15),
            color_by = "group")

#Plot top feature weight scores per modality
plot_weights(MOFAobject.ohsu,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.ohsu.factors <- get_factors(MOFAobject.ohsu, as.data.frame = T)
mofa.ohsu.weights <- get_weights(MOFAobject.ohsu, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.ohsu.top.features <- ExtractTopFeatures(mofa.ohsu.weights)

#subset to top factors and views
factors.vec <- as.character(unique(mofa.ohsu.top.features$factor))
subset_combos <- rbind(c('Factor1', "GEXP" ),
                       c('Factor1', "METH" ),
                       c('Factor2', "GEXP" ),
                       c('Factor2', "METH" ),
                       c('Factor3', "GEXP" ),
                       c('Factor3', "METH" ),
                       c('Factor4', "GEXP" ),
                       c('Factor4', "METH" ),
                       c('Factor5', "GEXP" ),
                       c('Factor5', "METH" ),
                       c('Factor6', "GEXP" ),
                       c('Factor6', "METH" ))

mofa.ohsu.top.features.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.ohsu.top.features.sub <- mofa.ohsu.top.features[mofa.ohsu.top.features$factor == c[1] & mofa.ohsu.top.features$view == c[2], ]
  mofa.ohsu.top.features.ev <- rbind(mofa.ohsu.top.features.sub, mofa.ohsu.top.features.ev)
}
