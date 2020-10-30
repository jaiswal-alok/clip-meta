#########################################################################
### CLIP Benchmarking                                                 ###
### Identification of cancer driver genes                             ###
#########################################################################
source("clip-meta/R/clip/clip_functions.R")

# known cancer driver genes
# https://www.sciencedirect.com/science/article/pii/S009286741830237X#app2
gold.std.genes <- readxl::read_xlsx("1-s2.0-S009286741830237X-mmc1.xlsx", sheet=2, skip=2)

#subset to cancer driver genes detected Breast cancer and PanCancer dataset
gold.std.genes <- gold.std.genes[gold.std.genes$Cancer %in% c("BRCA", "PANCAN"), ]
gold.std.genes <- unique(gold.std.genes$Gene)

#Load rCCS data from CLIP analysis for breast cancer cell lins
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
load("clip-meta/data/GEXP.RData")
#Log2 normalization of RNAseq studies
gexp.broad.log <-  apply(gexp.broad,2,function(x) log2(x+1))
gexp.gcsi.log <-  apply(gexp.gcsi,2,function(x) log2(x+1))
gexp.ohsu.log <-  apply(gexp.ohsu,2,function(x) log2(x+1))
gexp.list <- list(BROAD=gexp.broad.log, UHN= gexp.uhn, gCSI=gexp.gcsi.log, GDSC=gexp.gdsc,OHSU=gexp.ohsu.log )
gexp.list <- lapply(gexp.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
gexp.list <- lapply(gexp.list,QuantNormScale)
driver.GEXP.DE <- PerformDE(brca.subtype, gexp.list, "GEXP",  ref.drivers = gold.std.genes )

#Perform DE analysis on Functional data
load("clip-meta/data/FUNC.RData")
func.list <- list(BROAD_ACHILLES=func.broad_achilles, BROAD_AVANA=func.broad_avana,
                  GDSC=func.gdsc, UHN=func.uhn, DRIVE=func.drive )
func.list <- lapply(func.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
func.list <- lapply(func.list,QuantNormScale)
driver.FUNC.DE <- PerformDE(brca.subtype, func.list, "FUNC",  ref.drivers = gold.std.genes )

#Perform DE analysis in Methylation data
load("clip-meta/data/METH.RData")
meth.list <- list(BROAD=meth.broad,GDSC=meth.gdsc,OHSU=meth.ohsu )
meth.list <- lapply(meth.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
meth.list <- lapply(meth.list,QuantNormScale)
driver.METH.DE <- PerformDE(brca.subtype, meth.list, "METH",  ref.drivers = gold.std.genes )

#Perform DE analysis in PEXP data
load("clip-meta/data/PEXP.RData")
pexp.mghcc_breast.log <-  apply(pexp.mghcc_breast,2,log2)
pexp.list <- list(MGHCC=pexp.mghcc_breast.log, CCLE=pexp.ccle)
pexp.list <- lapply(pexp.list, function(data){
  data <- data[intersect(homo.genes.coding, rownames(data)),  sort(intersect(breast.cells,colnames(data)))]
})
pexp.list <- lapply(pexp.list,QuantNormScale)
driver.PEXP.DE <- PerformDE(brca.subtype, pexp.list, "PEXP",  ref.drivers = gold.std.genes )

#Perform DE analysis in CNV data
load("clip-meta/data/CNV.RData")
cnv.list <- list(BROAD=cnv.broad,GDSC=cnv.gdsc,OHSU=cnv.ohsu, gCSI=cnv.gcsi, UHN=cnv.uhn )
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




### Benchmarking with MOFA+ on CCLE dataset #####

library(reticulate)
use_python("/usr/local/bin/python3.8", required = TRUE)

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

load("clip-meta/data/CNV.RData")
load("clip-meta/data/GEXP.RData")
load("clip-meta/data/FUNC.RData")
load("clip-meta/data/PEXP.RData")
load("clip-meta/data/METH.RData")
load("clip-meta/data/MUT.RData")
load("clip-meta/data/PHOS.RData")
load("clip-meta/data/TAS.RData")

gexp.broad.log <- apply(gexp.broad,2,function(x) log2(x+1))

#
ccle.data <- list(CNV=as.matrix(cnv.broad),
                  GEXP=as.matrix(gexp.broad.log),
                  METH=as.matrix(meth.broad),
                  MUT=as.matrix(mut.broad),
                  PEXP=as.matrix(pexp.ccle),
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
ccle.data <- lapply(ccle.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(ccle.data, dim)

#Subset MUT and CNV data to variable genes
ccle.data$MUT <- ccle.data$MUT[rowSums(ccle.data$MUT)!=0, ]
ccle.data$CNV <- ccle.data$CNV[complete.cases(ccle.data$CNV), ]
ccle.data$CNV <- ccle.data$CNV[rowSums(ccle.data$CNV)!=0, ]

#Quantile normalize continuous datasets
ccle.data$GEXP <- QuantNormScale(ccle.data$GEXP)
ccle.data$METH <- QuantNormScale(ccle.data$METH)
ccle.data$FUNC_CRISPR <- QuantNormScale(ccle.data$FUNC_CRISPR)
ccle.data$FUNC_RNAI <- QuantNormScale(ccle.data$FUNC_RNAI)
ccle.data$PEXP <- QuantNormScale(ccle.data$PEXP)

cell.freq <- table(unlist(sapply(ccle.data, colnames)))
# n <- length(ccle.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
ccle.brca <- lapply(ccle.data, function(x){
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
lapply(ccle.brca, dim)

#Create MOFA object
MOFAobject.ccle <- create_mofa(ccle.brca)
data_opts <- get_default_data_options(MOFAobject.ccle)
model_opts <- get_default_model_options(MOFAobject.ccle)

#Add subtype information
sample_metadata <- MOFAobject.ccle@samples_metadata
sample_metadata <- merge(sample_metadata, brca.data.view[,c(2,6)], by.x="sample", by.y="CellName")
sample_metadata$subtype_three_receptor[is.na(sample_metadata$subtype_three_receptor)] <- "Unclassified"
sample_metadata <- sample_metadata[,-2]
colnames(sample_metadata)[2] <- "group"
sample_metadata$group <- as.factor(sample_metadata$group)

#Run MOFA with default parameters
MOFAobject.ccle <- create_mofa(ccle.brca, groups=sample_metadata$group)
data_opts <- get_default_data_options(MOFAobject.ccle)
model_opts <- get_default_model_options(MOFAobject.ccle)
train_opts <- get_default_training_options(MOFAobject.ccle)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ccle <- prepare_mofa(MOFAobject.ccle,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.ccle <- run_mofa(MOFAobject.ccle)

#Total Variance explained per modality
head(MOFAobject.ccle@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.ccle, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.ccle, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.ccle, title="")

#Variance explained per group, per factor, per modality
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_CCLE_Factor_VarExp.pdf", height = 5, width = 25)
p <- plot_variance_explained(MOFAobject.ccle, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
dev.off()

#factor correlation
plot_factor_cor(MOFAobject.ccle)

#Top Factors for each subgroup
#ER+: Factor1,2,6
#ERBB2+: Factor1,3,4,6
#TNBC: Factor3,4,5,7,8

#Plot factor scores per group
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_CCLE_FactorScores.pdf", height = 5, width = 25)
plot_factor(MOFAobject.ccle, 
            factor = c(1:15),
            color_by = "group")
dev.off()

#Plot top feature weight scores per modality
plot_weights(MOFAobject.ccle,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.ccle.factors <- get_factors(MOFAobject.ccle, as.data.frame = T)
mofa.ccle.weights <- get_weights(MOFAobject.ccle, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.ccle.top.features <- ExtractTopFeatures(mofa.ccle.weights)

#Overlap with known driver genes
mofa.ccle.top.features.gold <- mofa.ccle.top.features[mofa.ccle.top.features$FeatureID %in% gold.std.genes,  ]

#subset to top factors and views
factors.vec <- as.character(unique(mofa.ccle.top.features.gold$factor))
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

mofa.ccle.top.features.gold.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.ccle.top.features.gold.sub <- mofa.ccle.top.features.gold[mofa.ccle.top.features.gold$factor == c[1] & mofa.ccle.top.features.gold$view == c[2], ]
  mofa.ccle.top.features.gold.ev <- rbind(mofa.ccle.top.features.gold.ev, mofa.ccle.top.features.gold.sub)
}

### Benchmarking with MOFA+ on GDSC dataset #####

#
gdsc.data <- list(CNV=as.matrix(cnv.gdsc),
                  GEXP=as.matrix(gexp.gdsc),
                  METH=as.matrix(meth.gdsc),
                  MUT=as.matrix(mut.gdsc),
                  FUNC=as.matrix(func.gdsc))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
gdsc.data <- lapply(gdsc.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(gdsc.data, dim)

#Subset MUT and CNV data to variable genes
gdsc.data$MUT <- gdsc.data$MUT[rowSums(gdsc.data$MUT)!=0, ]
gdsc.data$CNV <- gdsc.data$CNV[complete.cases(gdsc.data$CNV), ]
gdsc.data$CNV <- gdsc.data$CNV[rowSums(gdsc.data$CNV)!=0, ]

#Quantile normalize continuous datasets
gdsc.data$GEXP <- QuantNormScale(gdsc.data$GEXP)
gdsc.data$METH <- QuantNormScale(gdsc.data$METH)
gdsc.data$FUNC <- QuantNormScale(gdsc.data$FUNC)

cell.freq <- table(unlist(sapply(gdsc.data, colnames)))
# n <- length(gdsc.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
gdsc.brca <- lapply(gdsc.data, function(x){
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
lapply(gdsc.brca, dim)

#Create MOFA object
MOFAobject.gdsc <- create_mofa(gdsc.brca)
data_opts <- get_default_data_options(MOFAobject.gdsc)
model_opts <- get_default_model_options(MOFAobject.gdsc)

#Add subtype information
sample_metadata <- MOFAobject.gdsc@samples_metadata
sample_metadata <- merge(sample_metadata, brca.data.view[,c(2,6)], by.x="sample", by.y="CellName")
sample_metadata$subtype_three_receptor[is.na(sample_metadata$subtype_three_receptor)] <- "Unclassified"
sample_metadata <- sample_metadata[,-2]
colnames(sample_metadata)[2] <- "group"
sample_metadata$group <- as.factor(sample_metadata$group)

#Run MOFA with default parameters
MOFAobject.gdsc <- create_mofa(gdsc.brca, groups=sample_metadata$group)
data_opts <- get_default_data_options(MOFAobject.gdsc)
model_opts <- get_default_model_options(MOFAobject.gdsc)
train_opts <- get_default_training_options(MOFAobject.gdsc)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.gdsc <- prepare_mofa(MOFAobject.gdsc,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.gdsc <- run_mofa(MOFAobject.gdsc)

#Total Variance explained per modality
head(MOFAobject.gdsc@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.gdsc, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.gdsc, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.gdsc, title="")

#Variance explained per group, per factor, per modality
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_gdsc_Factor_VarExp.pdf", height = 5, width = 25)
p <- plot_variance_explained(MOFAobject.gdsc, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
dev.off()

#factor correlation
plot_factor_cor(MOFAobject.gdsc)

#Plot factor scores per group
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_gdsc_FactorScores.pdf", height = 5, width = 25)
plot_factor(MOFAobject.gdsc, 
            factor = c(1:15),
            color_by = "group")
dev.off()

#Plot top feature weight scores per modality
plot_weights(MOFAobject.gdsc,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.gdsc.factors <- get_factors(MOFAobject.gdsc, as.data.frame = T)
mofa.gdsc.weights <- get_weights(MOFAobject.gdsc, as.data.frame = T)


#Top feature weights for each factor and overlap with driver genes
mofa.gdsc.top.features <- ExtractTopFeatures(mofa.gdsc.weights)

#Overlap with known driver genes
mofa.gdsc.top.features.gold <- mofa.gdsc.top.features[mofa.gdsc.top.features$FeatureID %in% gold.std.genes,  ]

#subset to top factors and views
factors.vec <- as.character(unique(mofa.gdsc.top.features.gold$factor))
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

mofa.gdsc.top.features.gold.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.gdsc.top.features.gold.sub <- mofa.gdsc.top.features.gold.ev[mofa.gdsc.top.features.gold.ev$factor == c[1] & mofa.gdsc.top.features.gold$view == c[2], ]
  mofa.gdsc.top.features.gold.ev <- rbind(mofa.gdsc.top.features.gold.sub, mofa.gdsc.top.features.gold.ev)
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
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_ohsu_Factor_VarExp.pdf", height = 5, width = 25)
p <- plot_variance_explained(MOFAobject.ohsu, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
dev.off()

#Plot factor scores per group
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_ohsu_FactorScores.pdf", height = 5, width = 25)
plot_factor(MOFAobject.ohsu, 
            factor = c(1:15),
            color_by = "group")
dev.off()

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
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_uhn_Factor_VarExp.pdf", height = 5, width = 25)
p <- plot_variance_explained(MOFAobject.uhn, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
dev.off()

#Plot factor scores per group
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_BRCA_uhn_FactorScores.pdf", height = 5, width = 25)
plot_factor(MOFAobject.uhn, 
            factor = c(1:15),
            color_by = "group")
dev.off()

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
     file="~/Desktop/FIMM_Work/CLIP_Datasets/MOFA_BRCA_Subtype.RData", compress = "bzip2")


load("~/Desktop/FIMM_Work/CLIP_Datasets/MOFA_BRCA_Subtype.RData")


#########################################################################
### ECHDC1 identificaton                                              ###
#########################################################################
##### Using a simpler approach #####
load("/clip-meta/data/METH.RData")

brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)
breast.cells <- brca.data.view$CellName

#Subset to coding genes and breast cell lines
meth.broad <- meth.broad[intersect(homo.genes.coding, rownames(meth.broad)),  
                         sort(intersect(breast.cells,colnames(meth.broad)))]
meth.gdsc <- meth.gdsc[intersect(homo.genes.coding, rownames(meth.gdsc)),
                       sort(intersect(breast.cells,colnames(meth.gdsc)))]
meth.ohsu <- meth.ohsu[intersect(homo.genes.coding, rownames(meth.ohsu)),
                       sort(intersect(breast.cells,colnames(meth.ohsu)))]

# Quantile normalize each dataset
meth.broad.qnorm <- QuantNormScale(meth.broad)
meth.gdsc.qnorm <- QuantNormScale(meth.gdsc)
meth.ohsu.qnorm <- QuantNormScale(meth.ohsu)

#Z-scaling
# Broad
meth.brca.broad.zscaled <- apply(t(meth.broad.qnorm),2,function(e) scale(e, center = T))
meth.brca.broad.zscaled <- t(meth.brca.broad.zscaled)
rownames(meth.brca.broad.zscaled) <- rownames(meth.broad.qnorm)
colnames(meth.brca.broad.zscaled) <- colnames(meth.broad.qnorm)

# GDSC
meth.brca.gdsc.zscaled <- apply(t(meth.gdsc.qnorm),2,function(e) scale(e, center = T))
meth.brca.gdsc.zscaled <- t(meth.brca.gdsc.zscaled)
rownames(meth.brca.gdsc.zscaled) <- rownames(meth.gdsc.qnorm)
colnames(meth.brca.gdsc.zscaled) <- colnames(meth.gdsc.qnorm)

#OHSU
meth.brca.ohsu.zscaled <- apply(t(meth.ohsu.qnorm),2,function(e) scale(e, center = T))
meth.brca.ohsu.zscaled <- t(meth.brca.ohsu.zscaled)
rownames(meth.brca.ohsu.zscaled) <- rownames(meth.ohsu.qnorm)
colnames(meth.brca.ohsu.zscaled) <- colnames(meth.ohsu.qnorm)

# Coont outlier frequency
meth.brca.broad.zscaled.long <- melt(as.matrix(meth.brca.broad.zscaled))
meth.brca.broad.zscaled.long$Keep = "F"
meth.brca.broad.zscaled.long$Keep[abs(meth.brca.broad.zscaled.long$value) >=sd(meth.brca.broad.zscaled.long$value, na.rm = T)*1.66] = "T"
meth.brca.broad.zscaled.long <- meth.brca.broad.zscaled.long[meth.brca.broad.zscaled.long$Var2 %in% breast.cells, ]
ccs.freq <- as.data.frame(table(meth.brca.broad.zscaled.long$Var1,meth.brca.broad.zscaled.long$Keep ))
ccs.freq <- ccs.freq[ccs.freq$Var2 != "F",]
ccs.freq[ccs.freq$Var1 == "ECHDC1", ]

# GDSC
meth.brca.gdsc.zscaled.long <- melt(as.matrix(meth.brca.gdsc.zscaled))
meth.brca.gdsc.zscaled.long$Keep = "F"
meth.brca.gdsc.zscaled.long$Keep[abs(meth.brca.gdsc.zscaled.long$value) >=sd(meth.brca.gdsc.zscaled.long$value, na.rm = T)*1.66] = "T"
meth.brca.gdsc.zscaled.long <- meth.brca.gdsc.zscaled.long[meth.brca.gdsc.zscaled.long$Var2 %in% breast.cells, ]
ccs.freq <- as.data.frame(table(meth.brca.gdsc.zscaled.long$Var1,meth.brca.gdsc.zscaled.long$Keep ))
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

##### MOFA+ CCLE #####

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)


load("clip-meta/data/CNV.RData")
load("clip-meta/data/GEXP.RData")
load("clip-meta/data/FUNC.RData")
load("clip-meta/data/PEXP.RData")
load("clip-meta/data/METH.RData")
load("clip-meta/data/MUT.RData")
load("clip-meta/data/PHOS.RData")
load("clip-meta/data/TAS.RData")

homo.genes <- read.table("clip-meta/data/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3


gexp.broad.log <- apply(gexp.broad,2,function(x) log2(x+1))
#
ccle.data <- list(CNV=as.matrix(cnv.broad),
                  GEXP=as.matrix(gexp.broad.log),
                  METH=as.matrix(meth.broad),
                  MUT=as.matrix(mut.broad),
                  PEXP=as.matrix(pexp.ccle),
                  FUNC_CRISPR=as.matrix(func.broad_avana),
                  FUNC_RNAI=as.matrix(func.broad_achilles))

#Subset to Breast cell lines and coding genes
brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)
breast.cells <- brca.data.view$CellName[!is.na(brca.data.view$subtype_three_receptor)]
ccle.data <- lapply(ccle.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(ccle.data, dim)

#Check density distribution
par(mfrow=c(2,4))
pDensity <- lapply(ccle.data, function(x){
  plot(density(as.matrix(x), na.rm = T))
})

#Subset MUT and CNV data to variable genes
ccle.data$MUT <- ccle.data$MUT[rowSums(ccle.data$MUT)!=0, ]
ccle.data$CNV <- ccle.data$CNV[complete.cases(ccle.data$CNV), ]
ccle.data$CNV <- ccle.data$CNV[rowSums(ccle.data$CNV)!=0, ]

#Quantile normalize continuous datasets
ccle.data$GEXP <- QuantNormScale(ccle.data$GEXP)
ccle.data$METH <- QuantNormScale(ccle.data$METH)
ccle.data$FUNC_CRISPR <- QuantNormScale(ccle.data$FUNC_CRISPR)
ccle.data$FUNC_RNAI <- QuantNormScale(ccle.data$FUNC_RNAI)
ccle.data$PEXP <- QuantNormScale(ccle.data$PEXP)

cell.freq <- table(unlist(sapply(ccle.data, colnames)))
# n <- length(ccle.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
ccle.brca <- lapply(ccle.data, function(x){
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
lapply(ccle.brca, dim)

#Create MOFA object
MOFAobject.ccle <- create_mofa(ccle.brca)
data_opts <- get_default_data_options(MOFAobject.ccle)
model_opts <- get_default_model_options(MOFAobject.ccle)

#Run MOFA with default parameters
train_opts <- get_default_training_options(MOFAobject.ccle)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ccle <- prepare_mofa(MOFAobject.ccle,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.ccle <- run_mofa(MOFAobject.ccle)


#Total Variance explained per modality
head(MOFAobject.ccle@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.ccle, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.ccle, x="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.ccle, title="")

#Variance explained per group, per factor, per modality
p <- plot_variance_explained(MOFAobject.ccle, x="view")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

#factor correlation
plot_factor_cor(MOFAobject.ccle)

#Top Factors for each subgroup
#ER+: Factor1,2,6
#ERBB2+: Factor1,3,4,6
#TNBC: Factor3,4,5,7,8

#Plot factor scores per group
plot_factor(MOFAobject.ccle, 
            factor = c(1:15))

#Plot top feature weight scores per modality
plot_weights(MOFAobject.ccle,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.ccle.factors <- get_factors(MOFAobject.ccle, as.data.frame = T)
mofa.ccle.weights <- get_weights(MOFAobject.ccle, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.ccle.top.features <- ExtractTopFeatures(mofa.ccle.weights)

#subset to top factors and views
factors.vec <- as.character(unique(mofa.ccle.top.features$factor))
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

mofa.ccle.top.features.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.ccle.top.features.sub <- mofa.ccle.top.features[mofa.ccle.top.features$factor == c[1] & mofa.ccle.top.features$view == c[2], ]
  mofa.ccle.top.features.ev <- rbind(mofa.ccle.top.features.ev, mofa.ccle.top.features.sub)
}

##### MOFA+ GDSC  #####

#
gdsc.data <- list(CNV=as.matrix(cnv.gdsc),
                  GEXP=as.matrix(gexp.gdsc),
                  METH=as.matrix(meth.gdsc),
                  MUT=as.matrix(mut.gdsc),
                  FUNC=as.matrix(func.gdsc))                  
#PHOS and TAS modality not included because of low dimensions

#Subset to Breast cell lines and coding genes
gdsc.data <- lapply(gdsc.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
  return(x)
} )
lapply(gdsc.data, dim)

#Subset MUT and CNV data to variable genes
gdsc.data$MUT <- gdsc.data$MUT[rowSums(gdsc.data$MUT)!=0, ]
gdsc.data$CNV <- gdsc.data$CNV[complete.cases(gdsc.data$CNV), ]
gdsc.data$CNV <- gdsc.data$CNV[rowSums(gdsc.data$CNV)!=0, ]

#Quantile normalize continuous datasets
gdsc.data$GEXP <- QuantNormScale(gdsc.data$GEXP)
gdsc.data$METH <- QuantNormScale(gdsc.data$METH)
gdsc.data$FUNC <- QuantNormScale(gdsc.data$FUNC)

cell.freq <- table(unlist(sapply(gdsc.data, colnames)))
# n <- length(gdsc.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subset to cell line profiled in all modalities
#For MOFA inout
gdsc.brca <- lapply(gdsc.data, function(x){
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
lapply(gdsc.brca, dim)

#Create MOFA object
MOFAobject.gdsc <- create_mofa(gdsc.brca)
data_opts <- get_default_data_options(MOFAobject.gdsc)
model_opts <- get_default_model_options(MOFAobject.gdsc)
train_opts <- get_default_training_options(MOFAobject.gdsc)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.gdsc <- prepare_mofa(MOFAobject.gdsc,
                                data_options = data_opts,
                                model_options = model_opts,
                                training_options = train_opts)
MOFAobject.gdsc <- run_mofa(MOFAobject.gdsc)

#Total Variance explained per modality
head(MOFAobject.gdsc@cache$variance_explained$r2_total[[1]])

#Variance explained per modality, per group
plot_variance_explained(MOFAobject.gdsc, plot_total = T)[[2]]
plot_variance_explained(MOFAobject.gdsc, x="group", y="factor", plot_total = T)

#Variance explained per factor per modality, per group
plot_variance_explained(MOFAobject.gdsc, title="")

#Variance explained per group, per factor, per modality
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_ECHDC1_gdsc_Factor_VarExp.pdf", height = 5, width = 6)
p <- plot_variance_explained(MOFAobject.gdsc, x="view")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
dev.off()

#factor correlation
plot_factor_cor(MOFAobject.gdsc)

#Plot factor scores per group
pdf("~/Desktop/FIMM_Work/MOFA/MOFA_ECHDC1_gdsc_FactorScores.pdf", height = 5, width = 25)
plot_factor(MOFAobject.gdsc, 
            factor = c(1:15))
dev.off()

#Plot top feature weight scores per modality
plot_weights(MOFAobject.gdsc,
             view = "GEXP",
             factor = c(1,2,6),
             nfeatures = 25,     
             scale = T)           


#Get factors and weights
mofa.gdsc.factors <- get_factors(MOFAobject.gdsc, as.data.frame = T)
mofa.gdsc.weights <- get_weights(MOFAobject.gdsc, as.data.frame = T)

#Top feature weights for each factor and overlap with driver genes
mofa.gdsc.top.features <- ExtractTopFeatures(mofa.gdsc.weights)

#subset to top factors and views
factors.vec <- as.character(unique(mofa.gdsc.top.features$factor))
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
mofa.gdsc.top.features.ev <- c()                       
for(i in 1:nrow(subset_combos)){
  c = subset_combos[i,]
  mofa.gdsc.top.features.sub <- mofa.gdsc.top.features[mofa.gdsc.top.features$factor == c[1] & mofa.gdsc.top.features$view == c[2], ]
  mofa.gdsc.top.features.ev <- rbind(mofa.gdsc.top.features.sub, mofa.gdsc.top.features.ev)
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
