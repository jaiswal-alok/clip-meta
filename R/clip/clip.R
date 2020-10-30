###########################################################################
### Step 1a: Quantile Normalization and Scaling of continuous datasets  ###
### Outlier Evidence scores (OES)                                      ###
###########################################################################
source("clip-meta/R/clip/clip_functions.R")

#coding genes
homo.genes <- read.table("clip-meta/data/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

### ------------- METHYLATION profiles ------------------ ###
load("clip-meta/data/METH.RData")
METH.QUANT.LIST <- list(meth.nci60, meth.ohsu, meth.gdsc, meth.broad)
names(METH.QUANT.LIST) <- c("METH.NCI60", "METH.OHSU", "METH.GDSC", "METH.BROAD")
METH.QUANT.LIST <- lapply(METH.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} ) #subset to coding genes
METH.QUANT.LIST <- lapply(METH.QUANT.LIST,QuantNormScale)
lapply(METH.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- TRANSCRIPTOME profiles ------------------ ###
load("clip-meta/data/GEXP.RData")
GEXP.QUANT.LIST <- list(gexp.uhn, gexp.ohsu, gexp.gdsc, gexp.broad, as.matrix(gexp.nci60), gexp.gcsi)
names(GEXP.QUANT.LIST) <- c("GEXP.UHN", "GEXP.OHSU", "GEXP.GDSC", "GEXP.BROAD", "GEXP.NCI60","GEXP.gCSI")
GEXP.QUANT.LIST <- lapply(GEXP.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
GEXP.QUANT.LIST <- lapply(GEXP.QUANT.LIST,QuantNormScale)
lapply(GEXP.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- PROTEOME profiles ------------------ ###
load("clip-meta/data/PEXP.RData")
PEXP.QUANT.LIST <- list(pexp.uw_tnbc, pexp.mghcc_breast, pexp.nci60, pexp.mipb_hgsoc, pexp.ccle)
names(PEXP.QUANT.LIST) <- c("PEXP.UW_TNBC", "PEXP.MGHCC_BREAST", "PEXP.NCI60", "PEXP.MPIB_HGSOC", "PEXP.CCLE")
PEXP.QUANT.LIST <- lapply(PEXP.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
PEXP.QUANT.LIST <- lapply(PEXP.QUANT.LIST,QuantNormScale)
lapply(PEXP.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- PHOSPHOPROTEOME profiles ------------------ ###
load("clip-meta/data/PHOS.RData")
#average intensities of multiple phospho-sites in each gene
phos.ccle.avg <- ResidueAvg(phos.ccle)
phos.mclp.avg <- ResidueAvg(phos.mclp)
phos.ohsu.avg <- ResidueAvg(phos.ohsu)
phos.uhn.avg <- ResidueAvg(phos.uhn)

PHOS.QUANT.LIST <- list(phos.uhn.avg, phos.ohsu.avg, phos.mclp.avg, phos.nci60, phos.ccle.avg)
names(PHOS.QUANT.LIST) <- c("PHOS.UHN", "PHOS.OHSU", "PHOS.MCLP", "PHOS.NCI60", "PHOS.BROAD")
PHOS.QUANT.LIST <- lapply(PHOS.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
PHOS.QUANT.LIST <- lapply(PHOS.QUANT.LIST,QuantNormScale)
lapply(PHOS.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- FUNCTIONAL profiles ------------------ ###
load("clip-meta/data/FUNC.RData")
FUNC.QUANT.LIST <- list(func.drive, func.broad_achilles, func.uhn, func.broad_avana, func.broad_gecko, func.broad_aml, func.gdsc)
names(FUNC.QUANT.LIST) <- c("FUNC.DRIVE", "FUNC.ACHILLES", "FUNC.UHN", "FUNC.AVANA", "FUNC.GECKO_BROAD", "FUNC.GECKO_AML", "FUNC.GDSC")
FUNC.QUANT.LIST <- lapply(FUNC.QUANT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
})
FUNC.QUANT.LIST <- lapply(FUNC.QUANT.LIST,QuantNormScale)
lapply(FUNC.QUANT.LIST, function(x) x[1:5, 1:5])

### ------------- TAS profiles ------------------ ###
load("clip-meta/data/TAS.RData")
TAS.QUANT.LIST <- list(tas.fimm, tas.ohsu, tas.gcsi, tas.gdsc, tas.broad_ccle, tas.broad_ctrp)
names(TAS.QUANT.LIST) <- c("TAS.FIMM", "TAS.OHSU", "TAS.gCSI", "TAS.GDSC", "TAS.BROAD_CCLE", "TAS.BROAD_CTRP")
TAS.QUANT.LIST <- lapply(TAS.QUANT.LIST, function(x){
  
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
})
TAS.QUANT.LIST <- lapply(TAS.QUANT.LIST,QuantNormScale)
lapply(TAS.QUANT.LIST, function(x) x[1:5, 1:5])

######################################################################
### Step 1b: Proportion scores in binary datasets                  ###
######################################################################

### ------------- CNV profiles ------------------ ###
load("clip-meta/data/CNV.RData")

# Binarize copy number profiles 
ccle.cnv.bin <- BinarizeCNV(cnv.broad, threshDEL = 1, threshAMP = 1)
gdsc.cnv.bin <- BinarizeCNV(cnv.gdsc, threshDEL = 0.5, threshAMP = 0.5)
gcsi.cnv.bin <- BinarizeCNV(cnv.gcsi, threshDEL = 0.5, threshAMP = 0.5)
nci60.cnv.bin <- BinarizeCNV(cnv.nci60, threshDEL = 0.4, threshAMP = 0.4)
uhn.cnv.bin <- BinarizeCNV(cnv.uhn, threshDEL = 0.4, threshAMP = 0.4)
ohsu.cnv.bin <- BinarizeCNV(cnv.ohsu, threshDEL = 0.5, threshAMP = 0.5)

# Estimate proportion scores
ccle.cnv.prop <- lapply(ccle.cnv.bin, CalculatePropScores)
gdsc.cnv.prop <- lapply(gdsc.cnv.bin, CalculatePropScores)
gcsi.cnv.prop <- lapply(gcsi.cnv.bin, CalculatePropScores)
nci60.cnv.prop <- lapply(nci60.cnv.bin, CalculatePropScores)
uhn.cnv.prop <- lapply(uhn.cnv.bin, CalculatePropScores)
ohsu.cnv.prop <- lapply(ohsu.cnv.bin, CalculatePropScores)


CNV.DEL.LIST <- list(uhn.cnv.prop$DEL, ccle.cnv.prop$DEL, nci60.cnv.prop$DEL, ohsu.cnv.prop$DEL, gcsi.cnv.prop$DEL, gdsc.cnv.prop$DEL)
CNV.DEL.LIST <- lapply(CNV.DEL.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
names(CNV.DEL.LIST) <- c("CNV_DEL.UHN", "CNV_DEL.BROAD", "CNV_DEL.NCI60", "CNV_DEL.OHSU", "CNV_DEL.gCSI", "CNV_DEL.GDSC")
lapply(CNV.DEL.LIST, function(x) x[1:5, 1:5])

CNV.AMP.LIST <- list(uhn.cnv.prop$AMP, ccle.cnv.prop$AMP, nci60.cnv.prop$AMP, ohsu.cnv.prop$AMP, gcsi.cnv.prop$AMP, gdsc.cnv.prop$AMP)
CNV.AMP.LIST <- lapply(CNV.AMP.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
names(CNV.AMP.LIST) <- c("CNV_AMP.UHN", "CNV_AMP.BROAD", "CNV_AMP.NCI60", "CNV_AMP.OHSU", "CNV_AMP.gCSI", "CNV_AMP.GDSC")
lapply(CNV.AMP.LIST, function(x) x[1:5, 1:5])

### ------------- MUt profiles ------------------ ###
load("clip-meta/data/MUT.RData")

# Estimate proportion scores
broad.mut.prop <- CalculatePropScores(mut.broad)
gdsc.mut.prop <- CalculatePropScores(mut.gdsc)
gcsi.mut.prop <- CalculatePropScores(mut.gcsi)
nci60.mut.prop <- CalculatePropScores(mut.nci60)
uhn.mut.prop <- CalculatePropScores(mut.uhn)
ohsu.mut.prop <- CalculatePropScores(mut.ohsu)

MUT.LIST <- list(uhn.mut.prop, gdsc.mut.prop, broad.mut.prop, ohsu.mut.prop, nci60.mut.prop,gcsi.mut.prop)
MUT.LIST <- lapply(MUT.LIST, function(x){
  y <- x[rownames(x) %in% homo.genes.coding,]
  return(y)
} )
names(MUT.LIST) <- c("MUT.UHN", "MUT.GDSC", "MUT.BROAD", "MUT.OHSU", "MUT.NCI60","MUT.gCSI")
lapply(MUT.LIST, function(x) x[1:5, 1:5])

save(file="clip-meta/processed_files/clip/Normalized_Modalities.RData",compress="bzip2",
     FUNC.QUANT.LIST, METH.QUANT.LIST, GEXP.QUANT.LIST, PEXP.QUANT.LIST, 
     PHOS.QUANT.LIST, TAS.QUANT.LIST, CNV.DEL.LIST, CNV.AMP.LIST, MUT.LIST)

######################################################################
### Step 2: CLIP - meta analysis of breast cancer cell lines       ###
######################################################################
load("clip-meta/processed_files/clip/Normalized_Modalities.RData")

# Breast cancer cell lines
brca.data.view <- readxl::read_xlsx("clip-meta/data/breastcancer_celllines.xlsx", sheet =1)
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

brca.rCCS.list <- list(FUNC    = brca.func.meta,
                       TAS     = brca.tas.meta,
                       METH    = brca.meth.meta, 
                       GEXP    = brca.gexp.meta,
                       PEXP    = brca.pexp.meta,
                       PHOS    = brca.phos.meta, 
                       CNV_DEL = brca.cnv_del.meta,
                       CNV_AMP = brca.cnv_amp.meta, 
                       MUT     = brca.mut.meta)

# rearrange CCS genes
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

# Add CCS direction
all.modes <- c("FUNC_UP", "TAS_UP", "TAS_DOWN", "METH_UP", "METH_DOWN", 
               "GEXP_UP", "GEXP_DOWN", "PEXP_UP", "PEXP_DOWN",
               "PHOS_UP", "PHOS_DOWN", "CNV_DEL", "CNV_AMP", "MUT")

brca.rCCSmat.list <- lapply(brca.rCCS.genes.list, function(x){

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

brca.rCCSmat <- ProcessCCSlist(brca.rCCSmat.list, keepMode = NA)
View(brca.rCCSmat[,c("GeneID","Sum")])

save(file="clip-meta/processed_files/clip/CLIP_rCCS_Breast.RData", brca.rCCSmat.list, brca.rCCSmat, compress="bzip2")


#########################################################################
### Run CLIP on all cell lines                                        ###
#########################################################################
load("clip-meta/processed_files/clip/Normalized_Modalities.RData")

# Cell lines with data from >=6 modalities
cells.data.view <- read.table("~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis/Analysis/For_Paper/All_Datasets_Cells_ViewCounts_v2.txt", sep="\t", header = T, stringsAsFactors = F)
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
