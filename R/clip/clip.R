###########################################################################
### Step 1a: Quantile Normalization and Scaling of continuous datasets  ###
### Outlier Evidence scores (OES)                                      ###
###########################################################################
source("clip-meta/R/clip/clip_functions.R")

# Read table of human genes and subset to coding genes
homo.genes <- read.table("clip-meta/data/Human_coding_genes.tsv", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

#Download all h5 format processed datasets from: https://doi.org/10.6084/m9.figshare.13473168
#Save datasets in clip-meta/data folder

### ------------- METHYLATION profiles ------------------ ###
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

#Quantile normalization
METH.QUANT.LIST <- list(meth.nci60, meth.ohsu, meth.sanger, meth.broad)
names(METH.QUANT.LIST) <- c("METH.NCI60", "METH.OHSU", "METH.SANGER", "METH.BROAD")
METH.QUANT.LIST <- lapply(METH.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} ) #subset to coding genes
METH.QUANT.LIST <- lapply(METH.QUANT.LIST,QuantNormScale) #
lapply(METH.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- TRANSCRIPTOME profiles ------------------ ###
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
GEXP.QUANT.LIST <- list(gexp.uhn, gexp.ohsu, gexp.sanger, gexp.broad, as.matrix(gexp.nci60), gexp.gcsi)
names(GEXP.QUANT.LIST) <- c("GEXP.UHN", "GEXP.OHSU", "GEXP.SANGER", "GEXP.BROAD", "GEXP.NCI60","GEXP.gCSI")
GEXP.QUANT.LIST <- lapply(GEXP.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} ) #subset to coding genes
GEXP.QUANT.LIST <- lapply(GEXP.QUANT.LIST,QuantNormScale)  #Quantile normalization
lapply(GEXP.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- PROTEOME profiles ------------------ ###

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

PEXP.QUANT.LIST <- list(pexp.uw_tnbc, pexp.mghcc_breast, pexp.nci60, pexp.mipb_hgsoc, pexp.broad)
names(PEXP.QUANT.LIST) <- c("PEXP.UW_TNBC", "PEXP.MGHCC_BREAST", "PEXP.NCI60", "PEXP.MPIB_HGSOC", "PEXP.BROAD")
PEXP.QUANT.LIST <- lapply(PEXP.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} ) #subset to coding genes
PEXP.QUANT.LIST <- lapply(PEXP.QUANT.LIST,QuantNormScale)  #Quantile normalization
lapply(PEXP.QUANT.LIST, function(x) x[1:5, 1:5]) 

### ------------- PHOSPHOPROTEOME profiles ------------------ ###
# laod Protein phohphorylation data
path = "clip-meta/data/"
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

#average intensities of multiple phospho-sites in each gene
phos.broad.avg <- ResidueAvg(phos.broad)
phos.mclp.avg <- ResidueAvg(phos.mclp)
phos.ohsu.avg <- ResidueAvg(phos.ohsu)
phos.uhn.avg <- ResidueAvg(phos.uhn)

PHOS.QUANT.LIST <- list(phos.uhn.avg, phos.ohsu.avg, phos.mclp.avg, phos.nci60, phos.broad.avg)
names(PHOS.QUANT.LIST) <- c("PHOS.UHN", "PHOS.OHSU", "PHOS.MCLP", "PHOS.NCI60", "PHOS.BROAD")
PHOS.QUANT.LIST <- lapply(PHOS.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} ) #subset to coding genes
PHOS.QUANT.LIST <- lapply(PHOS.QUANT.LIST,QuantNormScale)  #Quantile normalization
lapply(PHOS.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- FUNCTIONAL profiles ------------------ ###
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

FUNC.QUANT.LIST <- list(func.drive, func.broad_achilles, func.uhn, func.broad_avana, func.broad_gecko, func.broad_aml, func.sanger)
names(FUNC.QUANT.LIST) <- c("FUNC.DRIVE", "FUNC.ACHILLES", "FUNC.UHN", "FUNC.AVANA", "FUNC.GECKO_BROAD", "FUNC.GECKO_AML", "FUNC.SANGER")
FUNC.QUANT.LIST <- lapply(FUNC.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
}) #subset to coding genes
FUNC.QUANT.LIST <- lapply(FUNC.QUANT.LIST,QuantNormScale)  #Quantile normalization
lapply(FUNC.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- TAS profiles ------------------ ###
# load Target Addiction data
path = "clip-meta/data/"
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

TAS.QUANT.LIST <- list(tas.fimm, tas.ohsu, tas.gcsi, tas.sanger, tas.broad_ccle, tas.broad_ctrp)
names(TAS.QUANT.LIST) <- c("TAS.FIMM", "TAS.OHSU", "TAS.gCSI", "TAS.SANGER", "TAS.BROAD_CCLE", "TAS.BROAD_CTRP")
TAS.QUANT.LIST <- lapply(TAS.QUANT.LIST, function(x){
  
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
}) #subset to coding genes
TAS.QUANT.LIST <- lapply(TAS.QUANT.LIST,QuantNormScale)  #Quantile normalization
lapply(TAS.QUANT.LIST, function(x) x[1:5, 1:5])

######################################################################
### Step 1b: Proportion scores in binary datasets                  ###
######################################################################

### ------------- CNV profiles ------------------ ###
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


# Binarize copy number profiles 
broad.cnv.bin <- BinarizeCNV(cnv.broad, threshDEL = 1, threshAMP = 1)
sanger.cnv.bin <- BinarizeCNV(cnv.sanger, threshDEL = 0.5, threshAMP = 0.5)
gcsi.cnv.bin <- BinarizeCNV(cnv.gcsi, threshDEL = 0.5, threshAMP = 0.5)
nci60.cnv.bin <- BinarizeCNV(cnv.nci60, threshDEL = 0.4, threshAMP = 0.4)
uhn.cnv.bin <- BinarizeCNV(cnv.uhn, threshDEL = 0.4, threshAMP = 0.4)
ohsu.cnv.bin <- BinarizeCNV(cnv.ohsu, threshDEL = 0.5, threshAMP = 0.5)

# Estimate proportion scores
broad.cnv.prop <- lapply(broad.cnv.bin, CalculatePropScores)
sanger.cnv.prop <- lapply(sanger.cnv.bin, CalculatePropScores)
gcsi.cnv.prop <- lapply(gcsi.cnv.bin, CalculatePropScores)
nci60.cnv.prop <- lapply(nci60.cnv.bin, CalculatePropScores)
uhn.cnv.prop <- lapply(uhn.cnv.bin, CalculatePropScores)
ohsu.cnv.prop <- lapply(ohsu.cnv.bin, CalculatePropScores)

#Subset to set of coding genes
CNV.DEL.LIST <- list(uhn.cnv.prop$DEL, broad.cnv.prop$DEL, nci60.cnv.prop$DEL, ohsu.cnv.prop$DEL, gcsi.cnv.prop$DEL, sanger.cnv.prop$DEL)
CNV.DEL.LIST <- lapply(CNV.DEL.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
names(CNV.DEL.LIST) <- c("CNV_DEL.UHN", "CNV_DEL.BROAD", "CNV_DEL.NCI60", "CNV_DEL.OHSU", "CNV_DEL.gCSI", "CNV_DEL.SANGER")
lapply(CNV.DEL.LIST, function(x) x[1:5, 1:5])

CNV.AMP.LIST <- list(uhn.cnv.prop$AMP, broad.cnv.prop$AMP, nci60.cnv.prop$AMP, ohsu.cnv.prop$AMP, gcsi.cnv.prop$AMP, sanger.cnv.prop$AMP)
CNV.AMP.LIST <- lapply(CNV.AMP.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
names(CNV.AMP.LIST) <- c("CNV_AMP.UHN", "CNV_AMP.BROAD", "CNV_AMP.NCI60", "CNV_AMP.OHSU", "CNV_AMP.gCSI", "CNV_AMP.SANGER")
lapply(CNV.AMP.LIST, function(x) x[1:5, 1:5])

### ------------- MUT profiles ------------------ ###
# load Mutation profiles
path = "clip-meta/data/"
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


# Estimate proportion scores
broad.mut.prop <- CalculatePropScores(mut.broad)
sanger.mut.prop <- CalculatePropScores(mut.sanger)
gcsi.mut.prop <- CalculatePropScores(mut.gcsi)
nci60.mut.prop <- CalculatePropScores(mut.nci60)
uhn.mut.prop <- CalculatePropScores(mut.uhn)
ohsu.mut.prop <- CalculatePropScores(mut.ohsu)

#Subset to set of coding genes
MUT.LIST <- list(uhn.mut.prop, sanger.mut.prop, broad.mut.prop, ohsu.mut.prop, nci60.mut.prop,gcsi.mut.prop)
MUT.LIST <- lapply(MUT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
names(MUT.LIST) <- c("MUT.UHN", "MUT.SANGER", "MUT.BROAD", "MUT.OHSU", "MUT.NCI60","MUT.gCSI")
lapply(MUT.LIST, function(x) x[1:5, 1:5])

#Save quantile  normalized datasets
save(file="clip-meta/processed_files/clip/Normalized_Modalities.RData",compress="bzip2",
     FUNC.QUANT.LIST, METH.QUANT.LIST, GEXP.QUANT.LIST, PEXP.QUANT.LIST, 
     PHOS.QUANT.LIST, TAS.QUANT.LIST, CNV.DEL.LIST, CNV.AMP.LIST, MUT.LIST)

######################################################################
### Step 2: CLIP - meta analysis of breast cancer cell lines       ###
######################################################################
load("clip-meta/processed_files/clip/Normalized_Modalities.RData")

# Read Meta data on Breast cancer cell lines
brca.data.view <- readxl::read_xlsx("clip-meta/data/BRCA_Cellline_metadata.xlsx", sheet =1)
breast.cells <- brca.data.view$CellName

# Meta-analysis and selection of CCS genes in each modality
brca.mut.meta     <- IntegrateCCS(MUT.LIST, 
                                  Modality= "MUT",     
                                  cellNames= breast.cells, 
                                  SingleSelectThreshold= 10, 
                                  SelectIndx=4, 
                                  SelectThreshold=0.1)

brca.cnv_del.meta <- IntegrateCCS(CNV.DEL.LIST, 
                                  Modality= "CNV_DEL", 
                                  cellNames= breast.cells,
                                  SingleSelectThreshold= 10,
                                  SelectIndx=4, 
                                  SelectThreshold=0.1)

brca.cnv_amp.meta <- IntegrateCCS(CNV.AMP.LIST, 
                                  Modality= "CNV_AMP", 
                                  cellNames= breast.cells, 
                                  SingleSelectThreshold= 10, 
                                  SelectIndx=4, 
                                  SelectThreshold=0.1)


# Meta-analysis and selection of CCS genes in  each modality
# Rank-product integration of CCS scores 
brca.meth.meta <- IntegrateCCS(METH.QUANT.LIST, 
                               Modality= "METH", 
                               cellNames= breast.cells, 
                               SingleSelectThreshold= 100, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

brca.gexp.meta <- IntegrateCCS(GEXP.QUANT.LIST,
                               Modality= "GEXP", 
                               cellNames= breast.cells, 
                               SingleSelectThreshold= 100, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

brca.pexp.meta <- IntegrateCCS(PEXP.QUANT.LIST,
                               Modality= "PEXP", 
                               cellNames= breast.cells, 
                               SingleSelectThreshold= 100, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

brca.phos.meta <- IntegrateCCS(PHOS.QUANT.LIST,
                               Modality= "PHOS", 
                               cellNames= breast.cells, 
                               SingleSelectThreshold= 10, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

brca.func.meta <- IntegrateCCS(FUNC.QUANT.LIST,  
                               Modality = "FUNC", 
                               cellNames= breast.cells, 
                               SingleSelectThreshold= 200, 
                               SelectIndx=4, 
                               SelectThreshold=0.25)

brca.tas.meta <- IntegrateCCS(TAS.QUANT.LIST,
                               Modality= "TAS", 
                               cellNames= breast.cells, 
                               SingleSelectThreshold= 10, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

save(file = "clip-meta/processed_files/clip//CLIP_Integrated_Breast.RData",
     brca.mut.meta, brca.cnv_del.meta, brca.cnv_amp.meta, brca.meth.meta, brca.gexp.meta,
     brca.pexp.meta, brca.phos.meta, brca.func.meta, brca.tas.meta, compress="bzip2")

#########################################################################
### Step 3: CLIP - Rank Product Integration                          ####
### Identification of rCCS genes                                     ####
#########################################################################
load("clip-meta/processed_files/clip//CLIP_Integrated_Breast.RData")

# Retriever CCS genes for each cell line
brca.rCCS.list <- list(FUNC    = brca.func.meta,
                       TAS     = brca.tas.meta,
                       METH    = brca.meth.meta, 
                       GEXP    = brca.gexp.meta,
                       PEXP    = brca.pexp.meta,
                       PHOS    = brca.phos.meta, 
                       CNV_DEL = brca.cnv_del.meta,
                       CNV_AMP = brca.cnv_amp.meta, 
                       MUT     = brca.mut.meta)

brca.rCCS.genes.list <- list()
for ( c in breast.cells){
  cellName = c
  brca.rCCS.genes.list[[cellName]] <- list()

  for ( i in 1:length(brca.rCCS.list)){
    modality <- names(brca.rCCS.list)[i]
    brca.rCCS.genes.list[[cellName]][[modality]] <- list()
    
    #Table 1 = CCS_DOWN genes (for FUNC, CCS_UP)
    modality1 = paste(modality, "_T1", sep="")
    geneNames_1 = brca.rCCS.list[[i]]$cellIntResults[[c]]$Table1.genes
    if(length(geneNames_1)==0) geneNames_T1 <- cbind(modality=NA, GeneID=NA)
    if(length(geneNames_1)>0) geneNames_T1 <- cbind(modality=modality1, GeneID=geneNames_1)
    
    #Table 2 = CCS_UP genes (for FUNC, CCS_DOWN)
    modality2 = paste(modality, "_T2", sep="")
    geneNames_2 = brca.rCCS.list[[i]]$cellIntResults[[c]]$Table2.genes
    if(length(geneNames_2)==0) geneNames_T2 <- cbind(modality=NA, GeneID=NA)
    if(length(geneNames_2)>0) geneNames_T2 <- cbind(modality=modality2, GeneID=geneNames_2)
    
    geneNames_T1_T2 <- rbind(geneNames_T1, geneNames_T2)
    geneNames_T1_T2 <- as.data.frame(geneNames_T1_T2, stringsAsFactors = F)
    geneNames_T1_T2 <- geneNames_T1_T2[complete.cases(geneNames_T1_T2), ]
    
    brca.rCCS.genes.list[[cellName]][[modality]] <- geneNames_T1_T2
  }
}

# Add CCS_UP or CCS_DOWN direction
all.modes <- c("FUNC_UP", "TAS_UP", "TAS_DOWN", "METH_UP", "METH_DOWN", 
               "GEXP_UP", "GEXP_DOWN", "PEXP_UP", "PEXP_DOWN",
               "PHOS_UP", "PHOS_DOWN", "CNV_DEL", "CNV_AMP", "MUT")

brca.rCCSmat.list <- lapply(brca.rCCS.genes.list,RearrangeCCS)

brca.rCCSmat <- ProcessCCSlist(brca.rCCSmat.list, keepMode = NA)
View(brca.rCCSmat[,c("GeneID","Sum")])

save(file="clip-meta/processed_files/clip/CLIP_rCCS_Breast.RData", brca.rCCSmat.list, brca.rCCSmat, compress="bzip2")


#########################################################################
### Run CLIP on all cell lines                                        ###
#########################################################################
load("clip-meta/processed_files/clip/Normalized_Modalities.RData")

# Subest to cell lines with data from >=6 modalities
# Read table with data availability info on STUDY and MODALITY for each cell line
cells.data.view <- read.table("clip-meta/data/CELLLINExSTUDY.MODALITY.txt", sep="\t", header = T, stringsAsFactors = F)
rownames(cells.data.view) <- cells.data.view$CellName
cells.data.view <- unique(cells.data.view)
cells.data.view[,c(57:63,65)] <- apply(cells.data.view[,c(57:63,65)],2, as.numeric)
cells.data.view[,c(57:63,65)] <- apply(cells.data.view[,c(57:63,65)],2, function(x){
  x[x>0] <- 1
  return(x)
})
cells.data.view$TotalViewCount1 <- apply(cells.data.view[,c(57:63,65)],1,sum)
hist(cells.data.view$TotalViewCount1)
cells.data.view <- cells.data.view[cells.data.view$TotalViewCount1 >= 6, ]
compreh.cells <- cells.data.view$CellName #1047


# Meta-analysis and selection of CCS genes in each modality
all.mut.meta     <- IntegrateCCS(MUT.LIST, 
                                  Modality= "MUT",     
                                  cellNames= compreh.cells, 
                                  SingleSelectThreshold= 10, 
                                  SelectIndx=4, 
                                  SelectThreshold=0.1)

all.cnv_del.meta <- IntegrateCCS(CNV.DEL.LIST, 
                                  Modality= "CNV_DEL", 
                                  cellNames= compreh.cells,
                                  SingleSelectThreshold= 10,
                                  SelectIndx=4, 
                                  SelectThreshold=0.1)

all.cnv_amp.meta <- IntegrateCCS(CNV.AMP.LIST, 
                                  Modality= "CNV_AMP", 
                                  cellNames= compreh.cells, 
                                  SingleSelectThreshold= 10, 
                                  SelectIndx=4, 
                                  SelectThreshold=0.1)


# Meta-analysis and selection of CCS genes in  each modality
# Rank-product integration of CCS scores 
all.meth.meta <- IntegrateCCS(METH.QUANT.LIST, 
                               Modality= "METH", 
                               cellNames= compreh.cells, 
                               SingleSelectThreshold= 100, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

all.gexp.meta <- IntegrateCCS(GEXP.QUANT.LIST,
                               Modality= "GEXP", 
                               cellNames= compreh.cells, 
                               SingleSelectThreshold= 100, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

all.pexp.meta <- IntegrateCCS(PEXP.QUANT.LIST,
                               Modality= "PEXP", 
                               cellNames= compreh.cells, 
                               SingleSelectThreshold= 100, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

all.phos.meta <- IntegrateCCS(PHOS.QUANT.LIST,
                               Modality= "PHOS", 
                               cellNames= compreh.cells, 
                               SingleSelectThreshold= 10, 
                               SelectIndx=4, 
                               SelectThreshold=0.1)

all.func.meta <- IntegrateCCS(FUNC.QUANT.LIST,  
                               Modality = "FUNC", 
                               cellNames= compreh.cells, 
                               SingleSelectThreshold= 200, 
                               SelectIndx=4, 
                               SelectThreshold=0.25)

all.tas.meta <- IntegrateCCS(TAS.QUANT.LIST,
                              Modality= "TAS", 
                              cellNames= compreh.cells, 
                              SingleSelectThreshold= 10, 
                              SelectIndx=4, 
                              SelectThreshold=0.1)

save(file = "clip-meta/processed_files/clip//CLIP_Integrated_All.RData",
     all.mut.meta, all.cnv_del.meta, all.cnv_amp.meta, all.meth.meta, all.gexp.meta,
     all.pexp.meta, all.phos.meta, all.func.meta, all.tas.meta, compress="bzip2")


all.rCCS.list <- list(FUNC    = all.func.meta,
                       TAS     = all.tas.meta,
                       METH    = all.meth.meta, 
                       GEXP    = all.gexp.meta,
                       PEXP    = all.pexp.meta,
                       PHOS    = all.phos.meta, 
                       CNV_DEL = all.cnv_del.meta,
                       CNV_AMP = all.cnv_amp.meta, 
                       MUT     = all.mut.meta)

# rearrange CCS genes
all.rCCS.genes.list <- list()
for ( c in compreh.cells){
  cellName = c
  all.rCCS.genes.list[[cellName]] <- list()
  
  for ( i in 1:length(all.rCCS.list)){
    modality <- names(all.rCCS.list)[i]
    all.rCCS.genes.list[[cellName]][[modality]] <- list()
    
    #Table 1 = CCS_DOWN genes (for FUNC, CCS_UP)
    modality1 = paste(modality, "_T1", sep="")
    geneNames_1 = all.rCCS.list[[i]]$cellIntResults[[c]]$Table1.genes
    if(length(geneNames_1)==0) geneNames_T1 <- cbind(modality=NA, GeneID=NA)
    if(length(geneNames_1)>0) geneNames_T1 <- cbind(modality=modality1, GeneID=geneNames_1)
    
    #Table 2 = CCS_UP genes (for FUNC, CCS_DOWN)
    modality2 = paste(modality, "_T2", sep="")
    geneNames_2 = all.rCCS.list[[i]]$cellIntResults[[c]]$Table2.genes
    if(length(geneNames_2)==0) geneNames_T2 <- cbind(modality=NA, GeneID=NA)
    if(length(geneNames_2)>0) geneNames_T2 <- cbind(modality=modality2, GeneID=geneNames_2)
    
    geneNames_T1_T2 <- rbind(geneNames_T1, geneNames_T2)
    geneNames_T1_T2 <- as.data.frame(geneNames_T1_T2, stringsAsFactors = F)
    geneNames_T1_T2 <- geneNames_T1_T2[complete.cases(geneNames_T1_T2), ]
    
    all.rCCS.genes.list[[cellName]][[modality]] <- geneNames_T1_T2
  }
}

# Add CCS direction
all.modes <- c("FUNC_UP", "TAS_UP", "TAS_DOWN", "METH_UP", "METH_DOWN", 
               "GEXP_UP", "GEXP_DOWN", "PEXP_UP", "PEXP_DOWN",
               "PHOS_UP", "PHOS_DOWN", "CNV_DEL", "CNV_AMP", "MUT")

all.rCCSmat.list <- lapply(all.rCCS.genes.list, function(x){
  
  xMat <- do.call(rbind,x)
  xMat$Var <- 1
  xMat <- unique(xMat)
  
  #for FUNCTIONAL - opposite directions
  xMat$modality[xMat$modality %in% "FUNC_T2"] <-  "FUNC_DOWN"
  xMat$modality[xMat$modality %in% "FUNC_T1"] <-  "FUNC_UP"
  
  #For CNV 
  xMat$modality[xMat$modality %in% "CNV_AMP_T1"] <-  "CNV_AMP"
  xMat$modality[xMat$modality %in% "CNV_DEL_T1"] <-  "CNV_DEL"
  xMat$modality[xMat$modality %in% "MUT_T1"] <-  "MUT"
  
  xMat$modality <- gsub("_T1", "_DOWN", xMat$modality) 
  xMat$modality <- gsub("_T2", "_UP", xMat$modality) 
  
  #Remove FUNC_DOWN
  xMat <- xMat[xMat$modality %ni% "FUNC_DOWN",]
  
  #Remove TAS_DOWN
  xMat <- xMat[xMat$modality %ni% "TAS_DOWN",]
  
  #Add dummy variables for consistency
  xMat.modes <- unique(xMat$modality)
  xMat.modes.NotPresent <- setdiff(all.modes, xMat.modes)
  if(length(xMat.modes.NotPresent) > 0){
    for (p in xMat.modes.NotPresent){
      addRow <- c(p, xMat$GeneID[1], 0) 
      xMat <- rbind(xMat, addRow)
    }
  }
  
  #rCCS matrix for each cellline
  library(reshape2)
  xMat.Wide <- dcast(xMat, GeneID~modality, fill=0)
  xMat.Wide <- xMat.Wide[,c("GeneID", "FUNC_UP", "TAS_UP", "MUT","CNV_AMP","CNV_DEL","METH_UP","METH_DOWN", "GEXP_UP","GEXP_DOWN", "PEXP_UP", "PEXP_DOWN", "PHOS_UP", "PHOS_DOWN" )]
  xMat.Wide[,-1] <- apply(xMat.Wide[,-1],2,as.numeric)
  xMat.Wide$Sum <- rowSums(xMat.Wide[,-1])
  return(xMat.Wide)
})

all.rCCSmat <- ProcessCCSlist(all.rCCSmat.list, keepMode = NA)
View(all.rCCSmat[,c("GeneID","Sum")])

save(file="clip-meta/processed_files/clip/CLIP_rCCS_All.RData", all.rCCSmat.list, all.rCCSmat, compress="bzip2")
