###########################################################################
### Step 1a: Quantile Normalization and Scaling of continuous datasets  ###
### Cell line specificity (CCS) scores                                  ###
###########################################################################
source("clip-meta/R/clip/clip_functions.R")

#Subset each dataset to coding genes
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
} )
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

### ------------- CNV profiles ------------------ ###
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
### Step 2: CLIP integration of breast cancer cell lines           ###
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
### Step 3: CLIP integration - rCCS genes - Breast cancer cell lines ####
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
### CLIP Benchmarking                                                 ###
### Identification of cancer driver genes                             ###
#########################################################################
# known cancer driver genes
# https://www.sciencedirect.com/science/article/pii/S009286741830237X#app2
gold.std.genes <- readxl::read_xlsx("~/Desktop/FIMM_Work/CLIP_Review/1-s2.0-S009286741830237X-mmc1.xlsx", sheet=2, skip=2)
gold.std.genes <- readxl::read_xlsx("1-s2.0-S009286741830237X-mmc1.xlsx", sheet=2, skip=2)

#subset to cancer driver genes detected Breast cancer and PanCancer dataset
gold.std.genes <- gold.std.genes[gold.std.genes$Cancer %in% c("BRCA", "PANCAN"), ]
gold.std.genes <- unique(gold.std.genes$Gene)

#Load CCS data from CLIP analysis
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




### Benchmarking with MOFA on CCLE dataset #####

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

### Benchmarking with MOFA on GDSC dataset #####

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


### Benchmarking with MOFA on OHSU dataset #####

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

### Benchmarking with MOFA on UHN dataset #####

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
##### MOFA CCLE #####

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

load("~/Desktop/FIMM_Work/CLIP_Datasets/CNV.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/GEXP.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/FUNC.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/METH.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/MUT.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/PEXP.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/PHOS.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/TAS.RData")

homo.genes <- read.table("~/Desktop/Data/Housekeeping/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

#
ccle.data <- list(CNV=as.matrix(cnv.broad),
                  GEXP=as.matrix(gexp.broad),
                  METH=as.matrix(meth.broad),
                  MUT=as.matrix(mut.broad),
                  FUNC_CRISPR=as.matrix(func.broad_avana),
                  FUNC_RNAI=as.matrix(func.broad_achilles))
#PHOS=as.matrix(phos.ccle),
#TAS=as.matrix(tas.broad_ctrp))

#Subset to Breast cell lines and coding genes
load("~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis/BRCA_CSP_Gene_Mat_Prot.RData")
breast.cells <- colnames(brca.CSPgenes.Mat.Wide)[- 104]
ccle.data <- lapply(ccle.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
} )
lapply(ccle.data, dim)

#Check density distribution
par(mfrow=c(2,4))
pDensity <- lapply(ccle.data, function(x){
  plot(density(as.matrix(x), na.rm = T))
})

#Subset MUT data to variable genes
ccle.data$MUT <- ccle.data$MUT[rowSums(ccle.data$MUT)!=0, ]
ccle.data$CNV <- ccle.data$CNV[complete.cases(ccle.data$CNV), ]

#Log normalize transcriptomics data
ccle.data$GEXP <-  apply(ccle.data$GEXP,2,function(x) log2(x+1))

cell.freq <- table(unlist(sapply(ccle.data, colnames)))
# n <- length(ccle.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subet to common cells
ccle.order <- lapply(ccle.data, function(x){
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
lapply(ccle.order, dim)

MOFAobject.ccle <- create_mofa(ccle.order)
data_opts <- get_default_data_options(MOFAobject.ccle)
model_opts <- get_default_model_options(MOFAobject.ccle)

#7 factors
#model_opts$num_factors <- 7
train_opts <- get_default_training_options(MOFAobject.ccle)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ccle_7 <- prepare_mofa(MOFAobject.ccle,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)
MOFAobject.ccle_7 <- run_mofa(MOFAobject.ccle_7)
MOFAobject.ccle_15 <- MOFAobject.ccle_7
#10 fators
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject.ccle)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ccle_10 <- prepare_mofa(MOFAobject.ccle,
                                   data_options = data_opts,
                                   model_options = model_opts,
                                   training_options = train_opts)
MOFAobject.ccle_10 <- run_mofa(MOFAobject.ccle_10)

save(MOFAobject.ccle_7,MOFAobject.ccle_15, MOFAobject.ccle_10, file="~/Desktop/FIMM_Work/MOFA/CCLE_MOFA.RData", compress = "bzip2")

load("~/Desktop/FIMM_Work/MOFA/CCLE_MOFA.RData")
slotNames(MOFAobject.ccle)
names(MOFAobject.ccle@data)
dim(MOFAobject.ccle@data$Drugs$group1)
names(MOFAobject.ccle@expectations)
dim(MOFAobject.ccle@expectations$Z$group1)
dim(MOFAobject.ccle@expectations$W$GEXP)


plot_factor_cor(MOFAobject.ccle_7)


pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F15_Variance.pdf", height = 4, width = )
plot_variance_explained(MOFAobject.ccle_15, title="")
plot_variance_explained(MOFAobject.ccle_15, plot_total = T)[[2]]
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F7_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.ccle_7)
plot_variance_explained(MOFAobject.ccle_7, plot_total = T)[[2]]
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F10_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.ccle_10)
plot_variance_explained(MOFAobject.ccle_10, plot_total = T)[[2]]
dev.off()

#
pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F7_Factor_wts.pdf", height = 4, width = 6)
plot_factor(MOFAobject.ccle_7, factors = 2, color_by = "Factor2")
plot_weights(MOFAobject.ccle_7,view = "GEXP",factor = 2,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ccle_7,view = "METH",factor = 2,nfeatures = 10, scale = T )
plot_factor(MOFAobject.ccle_7, factors = 5, color_by = "Factor5")
plot_weights(MOFAobject.ccle_7,view = "GEXP",factor = 5,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ccle_7,view = "METH",factor = 5,nfeatures = 10, scale = T )
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F10_Factor_wts.pdf", height = 4, width = 6)
plot_factor(MOFAobject.ccle_10, factors = 2, color_by = "Factor2")
plot_weights(MOFAobject.ccle_10,view = "GEXP",factor = 2,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ccle_10,view = "METH",factor = 2,nfeatures = 10, scale = T )
plot_factor(MOFAobject.ccle_10, factors = 9, color_by = "Factor9")
plot_weights(MOFAobject.ccle_10,view = "GEXP",factor = 9,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ccle_10,view = "METH",factor = 9,nfeatures = 10, scale = T )
dev.off()


#
pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F7_Factor_top_wts.pdf", height = 5, width = 4)
plot_top_weights(MOFAobject.ccle_7, view = "GEXP",factor = 2,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ccle_7,view = "METH",factor = 2, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ccle_7, view = "GEXP",factor = 5,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ccle_7,view = "METH",factor = 5, nfeatures = 25,scale = T)
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F10_Factor_top_wts.pdf", height = 5, width = 4)
plot_top_weights(MOFAobject.ccle_10, view = "GEXP",factor = 2,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ccle_10,view = "METH",factor = 2, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ccle_10, view = "GEXP",factor = 9,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ccle_10,view = "METH",factor = 9, nfeatures = 25,scale = T)
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/CCLE_F15_Factor_top_wts.pdf", height = 5, width = 4)
for (i in 1:15){
  print(plot_top_weights(MOFAobject.ccle_15, view = "GEXP",factor = i,nfeatures = 25,scale = T))
  print(plot_top_weights(MOFAobject.ccle_15,view = "METH",factor = i, nfeatures = 25,scale = T))
}
dev.off()

##### MOFA GDSC #####

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

load("~/Desktop/FIMM_Work/CLIP_Datasets/CNV.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/GEXP.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/FUNC.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/METH.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/MUT.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/PEXP.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/PHOS.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/TAS.RData")

homo.genes <- read.table("~/Desktop/Data/Housekeeping/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

#
gdsc.data <- list(CNV=as.matrix(cnv.gdsc),
                  GEXP=as.matrix(gexp.gdsc),
                  METH=as.matrix(meth.gdsc),
                  MUT=as.matrix(mut.gdsc),
                  FUNC_CRISPR=as.matrix(func.gdsc))
#TAS=as.matrix(tas.gdsc))

#Subset to Breast cell lines and coding genes
load("~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis/BRCA_CSP_Gene_Mat_Prot.RData")
breast.cells <- colnames(brca.CSPgenes.Mat.Wide)[- 104]
gdsc.data <- lapply(gdsc.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
} )
lapply(gdsc.data, dim)

#Check density distribution
par(mfrow=c(2,4))
pDensity <- lapply(gdsc.data, function(x){
  plot(density(as.matrix(x), na.rm = T))
})

#Subset MUT data to variable genes
gdsc.data$MUT <- gdsc.data$MUT[rowSums(gdsc.data$MUT)!=0, ]
gdsc.data$CNV <- gdsc.data$CNV[complete.cases(gdsc.data$CNV), ]

cell.freq <- table(unlist(sapply(gdsc.data, colnames)))
# n <- length(gdsc.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subet to common cells
gdsc.order <- lapply(gdsc.data, function(x){
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
lapply(gdsc.order, dim)

MOFAobject.gdsc <- create_mofa(gdsc.order)
data_opts <- get_default_data_options(MOFAobject.gdsc)
model_opts <- get_default_model_options(MOFAobject.gdsc)

#7 factors
#model_opts$num_factors <- 7
train_opts <- get_default_training_options(MOFAobject.gdsc)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.gdsc_7 <- prepare_mofa(MOFAobject.gdsc,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)
MOFAobject.gdsc_7 <- run_mofa(MOFAobject.gdsc_7)
MOFAobject.gdsc_15 <- MOFAobject.gdsc_7

#10 fators
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject.gdsc)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.gdsc_10 <- prepare_mofa(MOFAobject.gdsc,
                                   data_options = data_opts,
                                   model_options = model_opts,
                                   training_options = train_opts)
MOFAobject.gdsc_10 <- run_mofa(MOFAobject.gdsc_10)

save(MOFAobject.gdsc_7, MOFAobject.gdsc_15,MOFAobject.gdsc_10, file="~/Desktop/FIMM_Work/MOFA/GDSC_MOFA.RData", compress = "bzip2")

load("~/Desktop/FIMM_Work/MOFA/GDSC_MOFA.RData")
slotNames(MOFAobject.gdsc)
names(MOFAobject.gdsc@data)
dim(MOFAobject.gdsc@data$Drugs$group1)
names(MOFAobject.gdsc@expectations)
dim(MOFAobject.gdsc@expectations$Z$group1)
dim(MOFAobject.gdsc@expectations$W$GEXP)

plot_factor_cor(MOFAobject.gdsc_7)

pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F7_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.gdsc_7)
plot_variance_explained(MOFAobject.gdsc_7, plot_total = T)[[2]]
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F10_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.gdsc_10)
plot_variance_explained(MOFAobject.gdsc_10, plot_total = T)[[2]]
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F15_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.gdsc_15)
plot_variance_explained(MOFAobject.gdsc_15, plot_total = T)[[2]]
dev.off()
#
pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F7_Factor_wts.pdf", height = 4, width = 6)
plot_factor(MOFAobject.gdsc_7, factors = 1, color_by = "Factor1")
plot_weights(MOFAobject.gdsc_7,view = "GEXP",factor = 1,nfeatures = 10, scale = T )
plot_weights(MOFAobject.gdsc_7,view = "METH",factor = 1,nfeatures = 10, scale = T )
plot_factor(MOFAobject.gdsc_7, factors = 6, color_by = "Factor6")
plot_weights(MOFAobject.gdsc_7,view = "GEXP",factor = 6,nfeatures = 10, scale = T )
plot_weights(MOFAobject.gdsc_7,view = "METH",factor = 6,nfeatures = 10, scale = T )
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F10_Factor_wts.pdf", height = 4, width = 6)
plot_factor(MOFAobject.gdsc_10, factors = 1, color_by = "Factor1")
plot_weights(MOFAobject.gdsc_10,view = "GEXP",factor = 1,nfeatures = 10, scale = T )
plot_weights(MOFAobject.gdsc_10,view = "METH",factor = 1,nfeatures = 10, scale = T )
plot_factor(MOFAobject.gdsc_10, factors = 6, color_by = "Factor6")
plot_weights(MOFAobject.gdsc_10,view = "GEXP",factor = 6,nfeatures = 10, scale = T )
plot_weights(MOFAobject.gdsc_10,view = "METH",factor = 6,nfeatures = 10, scale = T )
dev.off()

#
pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F7_Factor_top_wts.pdf", height = 5, width = 4)
plot_top_weights(MOFAobject.gdsc_7, view = "GEXP",factor = 1,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_7,view = "METH",factor = 1, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_7, view = "GEXP",factor = 6,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_7,view = "METH",factor = 6, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10, view = "GEXP",factor = 4,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10,view = "METH",factor = 4, nfeatures = 25,scale = T)
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F10_Factor_top_wts.pdf", height = 5, width = 4)
plot_top_weights(MOFAobject.gdsc_10, view = "GEXP",factor = 1,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10,view = "METH",factor = 1, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10, view = "GEXP",factor = 6,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10,view = "METH",factor = 6, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10, view = "GEXP",factor = 4,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.gdsc_10,view = "METH",factor = 4, nfeatures = 25,scale = T)
dev.off()


pdf("~/Desktop/FIMM_Work/MOFA/GDSC_F15_Factor_top_wts.pdf", height = 5, width = 4)
for (i in 1:15){
  print(plot_top_weights(MOFAobject.gdsc_15, view = "GEXP",factor = i,nfeatures = 25,scale = T))
  print(plot_top_weights(MOFAobject.gdsc_15,view = "METH",factor = i, nfeatures = 25,scale = T))
}
dev.off()


##### MOFA OHSU #####

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

load("~/Desktop/FIMM_Work/CLIP_Datasets/CNV.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/GEXP.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/FUNC.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/METH.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/MUT.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/PEXP.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/PHOS.RData")
load("~/Desktop/FIMM_Work/CLIP_Datasets/TAS.RData")

homo.genes <- read.table("~/Desktop/Data/Housekeeping/Homo_sapiens.gene_info", sep="\t", header=F, stringsAsFactors = F, quote="", fill=T)
homo.genes.coding.table <- homo.genes[homo.genes$V10 == "protein-coding", ]
homo.genes.coding.table$ENSEMBLID <- gsub("Ensembl:", "", sapply(strsplit(homo.genes.coding.table$V6, "\\|"), "[", 3))
homo.genes.coding <- homo.genes.coding.table$V3

#
ohsu.data <- list(CNV=as.matrix(cnv.ohsu),
                  GEXP=as.matrix(gexp.ohsu),
                  METH=as.matrix(meth.ohsu),
                  MUT=as.matrix(mut.ohsu))

#Subset to Breast cell lines and coding genes
load("~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis/BRCA_CSP_Gene_Mat_Prot.RData")
breast.cells <- colnames(brca.CSPgenes.Mat.Wide)[- 104]
ohsu.data <- lapply(ohsu.data, function(x){
  x <- x[, sort(intersect(colnames(x), breast.cells))]
  x <- x[sort(intersect(rownames(x), homo.genes.coding)), ]
} )
lapply(ohsu.data, dim)

#Check density distribution
par(mfrow=c(2,4))
pDensity <- lapply(ohsu.data, function(x){
  plot(density(as.matrix(x), na.rm = T))
})

#Subset MUT data to variable genes
ohsu.data$MUT <- ohsu.data$MUT[rowSums(ohsu.data$MUT)!=0, ]
ohsu.data$CNV <- ohsu.data$CNV[complete.cases(ohsu.data$CNV), ]

#Log normalize transcriptomics data
ohsu.data$GEXP <-  apply(ohsu.data$GEXP,2,function(x) log2(x+1))

cell.freq <- table(unlist(sapply(ohsu.data, colnames)))
# n <- length(ohsu.data)-1
# cell.freq <- cell.freq[cell.freq>n ]

#Arrange order and subet to common cells
ohsu.order <- lapply(ohsu.data, function(x){
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
lapply(ohsu.order, dim)

MOFAobject.ohsu <- create_mofa(ohsu.order)
data_opts <- get_default_data_options(MOFAobject.ohsu)
model_opts <- get_default_model_options(MOFAobject.ohsu)

#7 factors
#model_opts$num_factors <- 7
train_opts <- get_default_training_options(MOFAobject.ohsu)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ohsu_7 <- prepare_mofa(MOFAobject.ohsu,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)
MOFAobject.ohsu_7 <- run_mofa(MOFAobject.ohsu_7)
MOFAobject.ohsu_15 <- MOFAobject.ohsu_7

#10 fators
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject.ohsu)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
MOFAobject.ohsu_10 <- prepare_mofa(MOFAobject.ohsu,
                                   data_options = data_opts,
                                   model_options = model_opts,
                                   training_options = train_opts)
MOFAobject.ohsu_10 <- run_mofa(MOFAobject.ohsu_10)
save(MOFAobject.ohsu_7,MOFAobject.ohsu_15, MOFAobject.ohsu_10, file="~/Desktop/FIMM_Work/MOFA/OHSU_MOFA.RData", compress = "bzip2")

load("~/Desktop/FIMM_Work/MOFA/OHSU_MOFA.RData")
slotNames(MOFAobject.ohsu)
names(MOFAobject.ohsu@data)
dim(MOFAobject.ohsu@data$Drugs$group1)
names(MOFAobject.ohsu@expectations)
dim(MOFAobject.ohsu@expectations$Z$group1)
dim(MOFAobject.ohsu@expectations$W$GEXP)


plot_factor_cor(MOFAobject.ohsu_7)

pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F7_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.ohsu_7)
plot_variance_explained(MOFAobject.ohsu_7, plot_total = T)[[2]]
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F10_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.ohsu_10)
plot_variance_explained(MOFAobject.ohsu_10, plot_total = T)[[2]]
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F15_Variance.pdf", height = 4, width = 9)
plot_variance_explained(MOFAobject.ohsu_15)
plot_variance_explained(MOFAobject.ohsu_15, plot_total = T)[[2]]
dev.off()


#
pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F7_Factor_wts.pdf", height = 4, width = 6)
plot_factor(MOFAobject.ohsu_7, factors = 1, color_by = "Factor1")
plot_weights(MOFAobject.ohsu_7,view = "GEXP",factor = 1,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ohsu_7,view = "METH",factor = 1,nfeatures = 10, scale = T )
plot_factor(MOFAobject.ohsu_7, factors = 5, color_by = "Factor5")
plot_weights(MOFAobject.ohsu_7,view = "GEXP",factor = 5,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ohsu_7,view = "METH",factor = 5,nfeatures = 10, scale = T )
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F10_Factor_wts.pdf", height = 4, width = 6)
plot_factor(MOFAobject.ohsu_10, factors = 1, color_by = "Factor1")
plot_weights(MOFAobject.ohsu_10,view = "GEXP",factor = 1,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ohsu_10,view = "METH",factor = 1,nfeatures = 10, scale = T )
plot_factor(MOFAobject.ohsu_10, factors = 10, color_by = "Factor10")
plot_weights(MOFAobject.ohsu_10,view = "GEXP",factor = 10,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ohsu_10,view = "METH",factor = 10,nfeatures = 10, scale = T )
plot_factor(MOFAobject.ohsu_10, factors = 5, color_by = "Factor5")
plot_weights(MOFAobject.ohsu_10,view = "GEXP",factor = 5,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ohsu_10,view = "METH",factor = 5,nfeatures = 10, scale = T )
plot_factor(MOFAobject.ohsu_10, factors = 3, color_by = "Factor3")
plot_weights(MOFAobject.ohsu_10,view = "GEXP",factor = 3,nfeatures = 10, scale = T )
plot_weights(MOFAobject.ohsu_10,view = "METH",factor = 3,nfeatures = 10, scale = T )
dev.off()

#
pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F7_Factor_top_wts.pdf", height = 5, width = 4)
plot_top_weights(MOFAobject.ohsu_7, view = "GEXP",factor = 1,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_7,view = "METH",factor = 1, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_7, view = "GEXP",factor = 5,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_7,view = "METH",factor = 5, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_7, view = "GEXP",factor = 3,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_7,view = "METH",factor = 3, nfeatures = 25,scale = T)
dev.off()

pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F10_Factor_top_wts.pdf", height = 5, width = 4)
plot_top_weights(MOFAobject.ohsu_10, view = "GEXP",factor = 1,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10,view = "METH",factor = 1, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10, view = "GEXP",factor = 10,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10,view = "METH",factor = 10, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10, view = "GEXP",factor = 5,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10,view = "METH",factor = 5, nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10, view = "GEXP",factor = 3,nfeatures = 25,scale = T)
plot_top_weights(MOFAobject.ohsu_10,view = "METH",factor = 3, nfeatures = 25,scale = T)
dev.off()


pdf("~/Desktop/FIMM_Work/MOFA/OHSU_F15_Factor_top_wts.pdf", height = 5, width = 4)
for (i in 1:15){
  print(plot_top_weights(MOFAobject.ohsu_15, view = "GEXP",factor = i,nfeatures = 25,scale = T))
  print(plot_top_weights(MOFAobject.ohsu_15,view = "METH",factor = i, nfeatures = 25,scale = T))
}
dev.off()




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
