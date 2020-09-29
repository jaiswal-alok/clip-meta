##### load site-specific pre-processed data ######

load("clip-meta/data/CNV.RData")
load("clip-meta/data/GEXP.RData")
load("clip-meta/data/FUNC.RData")
load("clip-meta/data/METH.RData")
load("clip-meta/data/MUT.RData")
load("clip-meta/data/PEXP.RData")
load("clip-meta/data/PHOS.RData")
load("clip-meta/data/TAS.RData")

##### load functions and packages #####
source("clip-meta/R/reproducibility_analysis/reproducibility_analysis_functions.R")

######################################################################
#####     Assess reproducibility of METHYLATION profiles           ###
######################################################################
rm(list = ls())

load("clip-meta/data/METH.RData")

#
METH.COR.NCI60.OHSU <- CompareTwoStudies_Rcorr(meth.nci60, meth.ohsu, "spearman")
METH.COR.NCI60.GDSC  <- CompareTwoStudies_Rcorr(meth.nci60, meth.gdsc, "spearman")
METH.COR.NCI60.BROAD  <- CompareTwoStudies_Rcorr(meth.nci60, meth.broad, "spearman")
METH.COR.OHSU.GDSC  <- CompareTwoStudies_Rcorr(meth.ohsu, meth.gdsc, "spearman")
METH.COR.OHSU.BROAD  <- CompareTwoStudies_Rcorr(meth.ohsu, meth.broad, "spearman")
METH.COR.GDSC.BROAD  <- CompareTwoStudies_Rcorr(meth.gdsc, meth.broad, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("METH.COR", allvars)]
METH.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/meth_cor_long.rds", METH.COR.long,compress="bzip2")

######################################################################
#####     Assess reproducibility of MUTATION profiles              ###
######################################################################
load("clip-meta/data/MUT.RData")

#
# Matthew's Correlation coefficient = Pearson correlation coefficient estimated for two binary variables
MUT.COR.GDSC.CCLE <- CompareTwoStudies_Rcorr(mut.gdsc, mut.broad, "pearson")
MUT.COR.GDSC.NCI60 <- CompareTwoStudies_Rcorr(mut.gdsc, mut.nci60, "pearson")
MUT.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(mut.gdsc, mut.ohsu, "pearson")
MUT.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(mut.gdsc, mut.gcsi, "pearson")

MUT.COR.CCLE.NCI60 <- CompareTwoStudies_Rcorr(mut.broad, mut.nci60, "pearson")
MUT.COR.CCLE.OHSU <- CompareTwoStudies_Rcorr(mut.broad, mut.ohsu, "pearson")
MUT.COR.CCLE.gCSI <- CompareTwoStudies_Rcorr(mut.broad, mut.gcsi, "pearson")

MUT.COR.NCI60.OHSU <- CompareTwoStudies_Rcorr(mut.nci60, mut.ohsu, "pearson")
MUT.COR.NCI60.gCSI <- CompareTwoStudies_Rcorr(mut.nci60, mut.gcsi, "pearson")
MUT.COR.OHSU.gCSI <- CompareTwoStudies_Rcorr(mut.ohsu, mut.gcsi, "pearson")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("MUT.COR", allvars)]
MUT.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/mut_cor_long.rds", MUT.COR.long, compress="bzip2")

######################################################################
#####     Assess reproducibility of COPY NUMBER profiles           ###
######################################################################

load("clip-meta/data/CNV.RData")

CNV.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.gdsc, "spearman")
CNV.COR.UHN.CCLE <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.broad, "spearman")
CNV.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.nci60, "spearman")
CNV.COR.UHN.OHSU <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.ohsu, "spearman")
CNV.COR.UHN.gCSI <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.gcsi, "spearman")

CNV.COR.CCLE.GDSC <- CompareTwoStudies_Rcorr(cnv.broad, cnv.gdsc, "spearman")
CNV.COR.CCLE.NCI60 <- CompareTwoStudies_Rcorr(cnv.broad, cnv.nci60, "spearman")
CNV.COR.CCLE.OHSU <- CompareTwoStudies_Rcorr(cnv.broad, cnv.ohsu, "spearman")
CNV.COR.CCLE.gCSI <- CompareTwoStudies_Rcorr(cnv.broad, cnv.gcsi, "spearman")

CNV.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(cnv.nci60, cnv.gdsc, "spearman")
CNV.COR.NCI60.OHSU <- CompareTwoStudies_Rcorr(cnv.nci60, cnv.ohsu, "spearman")
CNV.COR.NCI60.gCSI <- CompareTwoStudies_Rcorr(cnv.nci60, cnv.gcsi, "spearman")

CNV.COR.OHSU.GDSC <- CompareTwoStudies_Rcorr(cnv.ohsu, cnv.gdsc, "spearman")
CNV.COR.OHSU.gCSI <- CompareTwoStudies_Rcorr(cnv.ohsu, cnv.gcsi, "spearman")
CNV.COR.gCSI.GDSC <- CompareTwoStudies_Rcorr(cnv.gcsi, cnv.gdsc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("CNV.COR", allvars)]
CNV.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/cnv_cor_long.rds", CNV.COR.long, compress="bzip2")




######################################################################
#####     Assess reproducibility of TRANSCRIPTOME profiles         ###
######################################################################

load("clip-meta/data/GEXP.RData")

#RNAseq vs. RNAseq
GEXP.COR.UHN.OHSU <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.ohsu, "spearman")
GEXP.COR.UHN.BROAD <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.broad, "spearman")
GEXP.COR.UHN.gCSI <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.gcsi, "spearman")
GEXP.COR.BROAD.OHSU <- CompareTwoStudies_Rcorr(gexp.broad, gexp.ohsu, "spearman")
GEXP.COR.BROAD.gCSI <- CompareTwoStudies_Rcorr(gexp.broad, gexp.gcsi, "spearman")
GEXP.COR.OHSU.gCSI <- CompareTwoStudies_Rcorr(gexp.ohsu, gexp.gcsi, "spearman")

#RNAseq vs. Array
GEXP.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.gdsc, "spearman")
GEXP.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.nci60, "spearman")
GEXP.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(gexp.broad, gexp.gdsc, "spearman")
GEXP.COR.BROAD.NCI60 <- CompareTwoStudies_Rcorr(gexp.broad, gexp.nci60, "spearman")
GEXP.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(gexp.gdsc, gexp.ohsu, "spearman")
GEXP.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(gexp.gdsc, gexp.gcsi, "spearman")
GEXP.COR.OHSU.NCI60 <- CompareTwoStudies_Rcorr(gexp.ohsu, gexp.nci60, "spearman")
GEXP.COR.NCI60.gCSI <- CompareTwoStudies_Rcorr(gexp.nci60, gexp.gcsi, "spearman")

#Array vs. Array
GEXP.COR.GDSC.NCI60 <- CompareTwoStudies_Rcorr(gexp.gdsc, gexp.nci60, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("GEXP.COR", allvars)]
GEXP.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/gexp_cor_long.rds", GEXP.COR.long,compress="bzip2")



######################################################################
#####     Assess reproducibility of PROTEOME profiles              ###
######################################################################

load("clip-meta/data/PEXP.RData")

# TMT labeled vs. TMT labeled
PEXP.COR.BROAD.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.ccle, pexp.mghcc_breast, "spearman")
PEXP.COR.GDSC.MGHCC_BREAST   <- CompareTwoStudies_Rcorr(pexp.gdsc, pexp.mghcc_breast, "spearman")
PEXP.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(pexp.ccle, pexp.gdsc, "spearman")

# Non-labeled vs. Non-labeled
PEXP.COR.UW_TNBC.NCI60 <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.nci60, "spearman")
PEXP.COR.NCI60.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.mipb_hgsoc, "spearman")

# TMT labeled  vs. Non-labeled
PEXP.COR.UW_TNBC.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.mghcc_breast, "spearman")
PEXP.COR.NCI60.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.mghcc_breast, "spearman")
PEXP.COR.UW_TNBC.BROAD <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.ccle, "spearman")
PEXP.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.gdsc, "spearman")
PEXP.COR.NCI60.BROAD <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.ccle, "spearman")
PEXP.COR.BROAD.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.ccle, pexp.mipb_hgsoc, "spearman")
PEXP.COR.GDSC.MPIB_HGSOC  <- CompareTwoStudies_Rcorr(pexp.gdsc, pexp.mipb_hgsoc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("PEXP.COR", allvars)]
PEXP.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/pexp_cor_long.rds", PEXP.COR.long,compress="bzip2")

######################################################################
#####     Assess reproducibility of PHOSPHOPROTEOME profiles       ###
######################################################################
load("clip-meta/data/PHOS.RData")

# RPPA vs. RPPA
PHOS.COR.UHN.OHSU <- CompareTwoStudies_Rcorr(phos.uhn, phos.ohsu, "spearman")
PHOS.COR.UHN.MCLP <- CompareTwoStudies_Rcorr(phos.uhn, phos.mclp, "spearman")
PHOS.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(phos.uhn, phos.nci60, "spearman")
PHOS.COR.UHN.BROAD <- CompareTwoStudies_Rcorr(phos.uhn, phos.ccle, "spearman")
PHOS.COR.OHSU.MCLP <- CompareTwoStudies_Rcorr(phos.ohsu, phos.mclp, "spearman")
PHOS.COR.OHSU.NCI60 <- CompareTwoStudies_Rcorr(phos.ohsu, phos.nci60, "spearman")
PHOS.COR.OHSU.BROAD <- CompareTwoStudies_Rcorr(phos.ohsu, phos.ccle, "spearman")
PHOS.COR.MCLP.NCI60 <- CompareTwoStudies_Rcorr(phos.mclp, phos.nci60, "spearman")
PHOS.COR.MCLP.BROAD <- CompareTwoStudies_Rcorr(phos.mclp, phos.ccle, "spearman")
PHOS.COR.BROAD.NCI60 <- CompareTwoStudies_Rcorr(phos.ccle, phos.nci60, "spearman")

# RPPA vs. MS-based
PHOS.COR.MCLP.GDSC <- CompareTwoStudies_Rcorr(phos.mclp, phos.gdsc, "spearman")
PHOS.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(phos.ccle, phos.gdsc, "spearman")
PHOS.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(phos.nci60, phos.gdsc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("PHOS.COR", allvars)]
PHOS.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/phos_cor_long.rds", PHOS.COR.long,compress="bzip2")



######################################################################
#####     Assess reproducibility of FUNCTIONAL profiles            ###
######################################################################
load("clip-meta/data/FUNC.RData")

# RNAi vs. RNAi
FUNC.COR.DRIVE.ACHILLES <- CompareTwoStudies_Rcorr(func.drive, func.broad_achilles, "spearman")
FUNC.COR.DRIVE.UHN <- CompareTwoStudies_Rcorr(func.drive, func.uhn, "spearman")
FUNC.COR.ACHILLES.UHN <- CompareTwoStudies_Rcorr(func.broad_achilles, func.uhn, "spearman")

# CRISPR vs. CRISPR
FUNC.COR.BROAD_AVANA.BROAD_GECKO <- CompareTwoStudies_Rcorr(func.broad_avana, func.broad_gecko, "spearman")
FUNC.COR.BROAD_AVANA.BROAD_AML <- CompareTwoStudies_Rcorr(func.broad_avana, func.broad_aml, "spearman")
FUNC.COR.BROAD_AVANA.GDSC <- CompareTwoStudies_Rcorr(func.broad_avana, func.gdsc, "spearman")
FUNC.COR.GDSC.BROAD_AML <- CompareTwoStudies_Rcorr(func.gdsc, func.broad_aml, "spearman")

# RNAi vs. CRISPR
FUNC.COR.DRIVE.BROAD_GECKO <- CompareTwoStudies_Rcorr(func.drive, func.broad_gecko, "spearman")
FUNC.COR.ACHILLES.BROAD_GECKO <- CompareTwoStudies_Rcorr(func.broad_achilles, gecko.broad.crispr, "spearman")
FUNC.COR.UHN.BROAD_GECKO <- CompareTwoStudies_Rcorr(func.uhn, func.broad_gecko, "spearman")
FUNC.COR.DRIVE.BROAD_AVANA <- CompareTwoStudies_Rcorr(func.drive, func.broad_avana, "spearman")
FUNC.COR.ACHILLES.BROAD_AVANA <- CompareTwoStudies_Rcorr(func.broad_achilles, func.broad_avana, "spearman")
FUNC.COR.UHN.BROAD_AVANA <- CompareTwoStudies_Rcorr(func.uhn, func.broad_avana, "spearman")
FUNC.COR.DRIVE.BROAD_AML <- CompareTwoStudies_Rcorr(func.drive, func.broad_aml, "spearman")
FUNC.COR.ACHILLES.BROAD_AML <- CompareTwoStudies_Rcorr(func.broad_achilles, func.broad_aml, "spearman")
FUNC.COR.DRIVE.GDSC <- CompareTwoStudies_Rcorr(func.drive, func.gdsc, "spearman")
FUNC.COR.ACHILLES.GDSC <- CompareTwoStudies_Rcorr(func.broad_achilles, func.gdsc, "spearman")
FUNC.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(func.uhn, func.gdsc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("FUNC.COR", allvars)]
FUNC.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/func_cor_long.rds", FUNC.COR.long,compress="bzip2")



######################################################################
#####     Assess reproducibility of DRUG SENSITIVITY profiles      ###
######################################################################

load("clip-meta/data/DSS.RData")

# CTG based vs. CTG based
DSS.COR.BROAD_CCLE.BROAD_CTRP <- CompareTwoStudies_Rcorr(dss.broad_ccle, dss.broad_ctrp, "spearman")
DSS.COR.BROAD_CCLE.FIMM <- CompareTwoStudies_Rcorr(dss.broad_ccle, dss.fimm, "spearman")
DSS.COR.BROAD_CCLE.gCSI <- CompareTwoStudies_Rcorr(dss.broad_ccle, dss.gcsi, "spearman")
DSS.COR.BROAD_CCLE.OHSU <- CompareTwoStudies_Rcorr(dss.broad_ccle, dss.ohsu, "spearman")
DSS.COR.BROAD_CTRP.FIMM <- CompareTwoStudies_Rcorr(dss.broad_ctrp, dss.fimm, "spearman")
DSS.COR.BROAD_CTRP.gCSI <- CompareTwoStudies_Rcorr(dss.broad_ctrp, dss.gcsi, "spearman")
DSS.COR.BROAD_CTRP.OHSU <- CompareTwoStudies_Rcorr(dss.broad_ctrp, dss.ohsu, "spearman")
DSS.COR.FIMM.gCSI <- CompareTwoStudies_Rcorr(dss.fimm, dss.gcsi, "spearman")
DSS.COR.FIMM.OHSU <- CompareTwoStudies_Rcorr(dss.fimm, dss.ohsu, "spearman")
DSS.COR.gCSI.OHSU <- CompareTwoStudies_Rcorr(dss.gcsi, dss.ohsu, "spearman")

# CTG based vs. Fluorophore based
DSS.COR.GDSC.BROAD_CCLE <- CompareTwoStudies_Rcorr(dss.gdsc, dss.broad_ccle, "spearman")
DSS.COR.GDSC.BROAD_CTRP <- CompareTwoStudies_Rcorr(dss.gdsc, dss.broad_ctrp, "spearman")
DSS.COR.GDSC.FIMM <- CompareTwoStudies_Rcorr(dss.gdsc, dss.fimm, "spearman")
DSS.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(dss.gdsc, dss.gcsi, "spearman")
DSS.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(dss.gdsc, dss.ohsu, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("DSS.COR", allvars)]
DSS.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/dss_cor_long.rds", DSS.COR.long, compress="bzip2")


######################################################################
#####     Assess reproducibility of TARGET ADDICTION profiles      ###
######################################################################
load("clip-meta/data/TAS.RData")

# CTG based vs. CTG based
TAS.COR.BROAD_CCLE.BROAD_CTRP <- CompareTwoStudies_Rcorr(tas.broad_ccle, tas.broad_ctrp, "spearman")
TAS.COR.BROAD_CCLE.FIMM <- CompareTwoStudies_Rcorr(tas.broad_ccle, tas.fimm, "spearman")
TAS.COR.BROAD_CCLE.gCSI <- CompareTwoStudies_Rcorr(tas.broad_ccle, tas.gcsi, "spearman")
TAS.COR.BROAD_CCLE.OHSU <- CompareTwoStudies_Rcorr(tas.broad_ccle, tas.ohsu, "spearman")
TAS.COR.BROAD_CTRP.FIMM <- CompareTwoStudies_Rcorr(tas.broad_ctrp, tas.fimm, "spearman")
TAS.COR.BROAD_CTRP.gCSI <- CompareTwoStudies_Rcorr(tas.broad_ctrp, tas.gcsi, "spearman")
TAS.COR.BROAD_CTRP.OHSU <- CompareTwoStudies_Rcorr(tas.broad_ctrp, tas.ohsu, "spearman")
TAS.COR.FIMM.gCSI <- CompareTwoStudies_Rcorr(tas.fimm, tas.gcsi, "spearman")
TAS.COR.FIMM.OHSU <- CompareTwoStudies_Rcorr(tas.fimm, tas.ohsu, "spearman")
TAS.COR.gCSI.OHSU <- CompareTwoStudies_Rcorr(tas.gcsi, tas.ohsu, "spearman")


# CTG based vs. Fluorophore based
TAS.COR.GDSC.BROAD_CCLE <- CompareTwoStudies_Rcorr(tas.gdsc, tas.broad_ccle, "spearman")
TAS.COR.GDSC.BROAD_CTRP <- CompareTwoStudies_Rcorr(tas.gdsc, tas.broad_ctrp, "spearman")
TAS.COR.GDSC.FIMM <- CompareTwoStudies_Rcorr(tas.gdsc, tas.fimm, "spearman")
TAS.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(tas.gdsc, tas.gcsi, "spearman")
TAS.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(tas.gdsc, tas.ohsu, "spearman")


#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("TAS.COR", allvars)]
TAS.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/tas_cor_long.rds", TAS.COR.long, compress="bzip2")

######################################################################
#####     Plot Reproduciblity measures                             ###
######################################################################

#load correlation data frame for each modality
listF <- list.files("clip-meta/processed_files/reprod_analysis/")
listF <- listF[grep("_long.rds", listF)]
cor.long.list <- lapply(listF, function(fl){
  pathF <- paste("clip-meta/processed_files/reprod_analysis/", fl, sep = "")
  f.rds <- readRDS(pathF)
  return(f.rds)
})

#
all.modal.cor.long <- do.call(rbind, cor.long.list)
all.modal.cor.long$ID[all.modal.cor.long$ID == 1] <- "Identical"
all.modal.cor.long$ID[all.modal.cor.long$ID == 0] <- "Non-identical"

all.modal.cor.long$DataNumber = 1
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "MUT"] = 2
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "CNV"] = 3
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "GEXP"] = 4
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "PEXP"] = 5
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "PHOS"] = 6
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "FUNC"] = 7
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "DSS"] = 8
all.modal.cor.long$DataNumber[all.modal.cor.long$Datatype == "TAS"] = 9

# Identical vs. Non-identical cell lines for all modalities
library(ggplot2)
colors.plot <- c('#999999','#E69F00')
p <- ggplot(all.modal.cor.long, aes(factor(DataNumber), Cor, fill=colors.plot))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p <- p + geom_violin(aes(fill = factor(ID)), scale = "width", trim=TRUE, draw_quantiles =  0.5, width=0.6) 
print(p)



######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     METHYLATION profiles                                     ###
######################################################################

# load correlation data
METH.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/meth_cor_long.rds")
METH.COR.long <- METH.COR.long[METH.COR.long$ID ==1, ] #subset to identical cell lines

# add platform info
#GDSC and NCI60 = ARRAY 450k
#OHSU = ARRAY 27k
#BROAD = Bisulfite sequencing
METH.COR.long$DATA1_TYPE <- "ARRAY_450"
METH.COR.long$DATA2_TYPE <- "ARRAY_450"
METH.COR.long$DATA1_TYPE[grep("BROAD", METH.COR.long$Dataset1)] <- "BISULFITE"
METH.COR.long$DATA2_TYPE[grep("BROAD", METH.COR.long$Dataset2)] <- "BISULFITE"
METH.COR.long$DATA1_TYPE[grep("OHSU", METH.COR.long$Dataset1)] <- "ARRAY_27"
METH.COR.long$DATA2_TYPE[grep("OHSU", METH.COR.long$Dataset2)] <- "ARRAY_27"
METH.COR.long$DATATYPE_CLASS <- paste(METH.COR.long$DATA1_TYPE, METH.COR.long$DATA2_TYPE, sep="_" )
METH.COR.long$DATATYPE_SAME <- METH.COR.long$DATA1_TYPE == METH.COR.long$DATA2_TYPE

ARRAY_450.IDx <- grep("ARRAY_450", METH.COR.long$DATATYPE_CLASS)
ARRAY_27.IDx <- grep("ARRAY_27", METH.COR.long$DATATYPE_CLASS)
BISULFITE.IDx <- grep("BISULFITE", METH.COR.long$DATATYPE_CLASS)

# Statistical comparisons
# 450/450 vs. 450/bisulfite
x1 <-  grep("ARRAY_450_ARRAY_450", METH.COR.long$DATATYPE_CLASS)
x2 <-  grep("ARRAY_450_BISULFITE", METH.COR.long$DATATYPE_CLASS)
x <- wilcox.test(METH.COR.long$Cor[x1], METH.COR.long$Cor[x2])
x$p.value

# array/array vs. array/bisulfite
y1 <-  grep("ARRAY_27_ARRAY_450|ARRAY_450_ARRAY_27", METH.COR.long$DATATYPE_CLASS)
y2 <-  grep("ARRAY_27_BISULFITE", METH.COR.long$DATATYPE_CLASS)
y <- wilcox.test(METH.COR.long$Cor[y1], METH.COR.long$Cor[y2])
y$p.value

# Corrplots
MakeCorrPlot(METH.COR.long, order= c("BROAD", "OHSU", "NCI60", "SANGER") )

######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     FUNCTIONAL profiles                                      ###
######################################################################
# load correlation data
FUNC.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/func_cor_long.rds")
FUNC.COR.long <- FUNC.COR.long[FUNC.COR.long$ID ==1, ] #subset to identical cell lines

# add platform info
#BROAD_GECKO,BROAD_AVANA,BROAD_AML, GDSC = CRISPR
#ACHILLES, DRIVE, UHN = ARNAi
FUNC.COR.long$DATA1_TYPE <- "RNAi"
FUNC.COR.long$DATA2_TYPE <- "RNAi"
FUNC.COR.long$DATA1_TYPE[grep("BROAD_GECKO", FUNC.COR.long$Dataset1)] <- "CRISPR"
FUNC.COR.long$DATA2_TYPE[grep("BROAD_GECKO", FUNC.COR.long$Dataset2)] <- "CRISPR"
FUNC.COR.long$DATA1_TYPE[grep("BROAD_AVANA", FUNC.COR.long$Dataset1)] <- "CRISPR"
FUNC.COR.long$DATA2_TYPE[grep("BROAD_AVANA", FUNC.COR.long$Dataset2)] <- "CRISPR"
FUNC.COR.long$DATA1_TYPE[grep("BROAD_AML", FUNC.COR.long$Dataset1)] <- "CRISPR"
FUNC.COR.long$DATA2_TYPE[grep("BROAD_AML", FUNC.COR.long$Dataset2)] <- "CRISPR"
FUNC.COR.long$DATA1_TYPE[grep("GDSC", FUNC.COR.long$Dataset1)] <- "CRISPR"
FUNC.COR.long$DATA2_TYPE[grep("GDSC", FUNC.COR.long$Dataset2)] <- "CRISPR"
FUNC.COR.long$DATATYPE_CLASS <- paste(FUNC.COR.long$DATA1_TYPE, FUNC.COR.long$DATA2_TYPE, sep="_" )
table(FUNC.COR.long$DATATYPE_CLASS )


#Statistical comparisons
#select index
allRNAi_IDx <- grep("RNAi", FUNC.COR.long$DATATYPE_CLASS) #All RNAi
allCRISPR_IDx <- grep("CRISPR", FUNC.COR.long$DATATYPE_CLASS) #All CRISPR

RNAi_RNAi_IDx <- setdiff(allRNAi_IDx,allCRISPR_IDx) 
Crspr_Crspr_IDx <- setdiff(allCRISPR_IDx,allRNAi_IDx) 
RNAi_Crspr_IDx <- intersect(allRNAi_IDx,allCRISPR_IDx) 

#RNAi/RNAi vs. CRISPR/CRISPR
wilcox.test(FUNC.COR.long$Cor[RNAi_RNAi_IDx],FUNC.COR.long$Cor[Crspr_Crspr_IDx])

##RNAi/CRISPR vs. CRISPR/CRISPR
wilcox.test(FUNC.COR.long$Cor[RNAi_Crspr_IDx],FUNC.COR.long$Cor[Crspr_Crspr_IDx])

##RNAi/RNAi vs. RNAi/CRISPR 
wilcox.test(FUNC.COR.long$Cor[RNAi_RNAi_IDx],FUNC.COR.long$Cor[RNAi_Crspr_IDx])

#Make Corrplot
MakeCorrPlot(FUNC.COR.long, order= c("ACHILLES", "DRIVE", "UHN", "BROAD_AVANA", "BROAD_GECKO", "BROAD_AML", "GDSC") )



######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     PROTEOME profiles                                        ###
######################################################################
# load correlation data
PEXP.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/pexp_cor_long.rds")
PEXP.COR.long <- PEXP.COR.long[PEXP.COR.long$ID == 1, ] #subset to identical cell lines

# add platform info
#GDSC, BROAD, MGHCC_BREAST  = TMT-labelled
#NCI60, UW_TNBC,  MPIB_HGSOC = Non-labelled
PEXP.COR.long$DATATYPE_CLASS <- paste(PEXP.COR.long$Dataset1, PEXP.COR.long$Dataset2, sep="_vs_" )
PEXP.COR.long$DATA1_TYPE <- "NL"
PEXP.COR.long$DATA2_TYPE <- "NL"
PEXP.COR.long$DATA1_TYPE[grep("GDSC", PEXP.COR.long$Dataset1)] <- "TMT"
PEXP.COR.long$DATA2_TYPE[grep("GDSC", PEXP.COR.long$Dataset2)] <- "TMT"
PEXP.COR.long$DATA1_TYPE[grep("BROAD", PEXP.COR.long$Dataset1)] <- "TMT"
PEXP.COR.long$DATA2_TYPE[grep("BROAD", PEXP.COR.long$Dataset2)] <- "TMT"
PEXP.COR.long$DATA1_TYPE[grep("MGHCC_BREAST", PEXP.COR.long$Dataset1)] <- "TMT"
PEXP.COR.long$DATA2_TYPE[grep("MGHCC_BREAST", PEXP.COR.long$Dataset2)] <- "TMT"
PEXP.COR.long$DATATYPE_CLASS <- paste(PEXP.COR.long$DATA1_TYPE, PEXP.COR.long$DATA2_TYPE, sep="_" )
table(PEXP.COR.long$DATATYPE_CLASS )

#Statistical comparisons
tmt_tmt_IDx <- grep("TMT_TMT", PEXP.COR.long$DATATYPE_CLASS) 
tmt_nl_IDx <- grep("NL_TMT|TMT_NL", PEXP.COR.long$DATATYPE_CLASS) 
nl_nl_IDx <- grep("NL_NL", PEXP.COR.long$DATATYPE_CLASS)

#TMT/TMT vs. TMT/NL
wilcox.test(PEXP.COR.long$Cor[tmt_tmt_IDx],PEXP.COR.long$Cor[tmt_nl_IDx])$p.value
#TMT/TMT vs. NL/NL
wilcox.test(PEXP.COR.long$Cor[tmt_tmt_IDx],PEXP.COR.long$Cor[nl_nl_IDx])$p.value
#NL/NL vs. TML/NL
wilcox.test(PEXP.COR.long$Cor[tmt_nl_IDx],PEXP.COR.long$Cor[nl_nl_IDx])$p.value

#Make Violin Plot
PEXP.COR.long$DATATYPE_NUM <- 1
PEXP.COR.long$DATATYPE_NUM[PEXP.COR.long$DATATYPE_CLASS == "NL_TMT"] <- 2
PEXP.COR.long$DATATYPE_NUM[PEXP.COR.long$DATATYPE_CLASS == "TMT_NL"] <- 2
PEXP.COR.long$DATATYPE_NUM[PEXP.COR.long$DATATYPE_CLASS == "TMT_TMT"] <- 3
p <- ggplot(PEXP.COR.long, aes(factor(DATATYPE_NUM), Cor))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p <- p + geom_violin(aes(fill = factor(DATATYPE_NUM)), scale = "width", trim=FALSE, colour = "gray30") 
cols <- c("1" = "powderblue", "2" = "seagreen3","3" = "lightslateblue")
p <- p + scale_fill_manual(values = cols)
p + stat_summary(fun = "mean", geom = "point", shape = 24, size = 3, color = "gray4",fill="white")

#Make Corrplot
MakeCorrPlot(PEXP.COR.long, order= c("BROAD", "MGHCC_BREAST", "SANGER", "NCI60", "UW_TNBC", "MPIB_HGSOC") )


######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     PROTEOME profiles - Coefficient of variation CCV)        ###
######################################################################

# CV correlation between UW_TNBC and HAAS.BRCA
load("clip-meta/data/PEXP.RData")
housekeep.genes <- read.table("clip-meta/data/housekeeping-genes.txt", sep="\t", header = T, stringsAsFactors = F)


##
ccle.prot <- read.csv("~/Downloads/summed_sn_non_normalized.csv", stringsAsFactors = F)
ccle.prot <- ccle.prot[, c(2,47:466)]
colnames(ccle.prot) <- sapply(strsplit(colnames(ccle.prot), "\\_"), "[[", 1)
ccle.prot[1:5, 1:10]
ccle.prot[,-1] <- apply(ccle.prot[,-1],2,as.numeric)

ccle.prot.list <- split(ccle.prot[,-1], ccle.prot$Gene)
ccle.prot.mat <- lapply(ccle.prot.list, function(x) apply(x,2, function(y) mean(y, na.rm=T)))
ccle.prot.mat <- do.call(rbind, ccle.prot.mat)

#Check correlation of cell line duplicates
dupProt <- ccle.prot.mat[, grep("SW948|CAL120|HCT15", colnames(ccle.prot.mat))]
rcorr(dupProt, type = "spearman")

#Remove duplcates
ccle.prot.mat <- ccle.prot.mat[, setdiff(colnames(ccle.prot.mat), c("SW948.1", "CAL120.1", "HCT15.1"))]
ccle.prot.mat <- ccle.prot.mat[sort(rownames(ccle.prot.mat)), sort(colnames(ccle.prot.mat))]

#
PlotCV(pexp.uw_tnbc, pexp.nci60, housekeep.genes, "UW_TNBC (CV)", "NCI60 (CV)")
PlotCV(pexp.mipb_hgsoc, pexp.nci60, housekeep.genes, "MIPB_HGSOC (CV)", "NCI60 (CV)")

PlotCV(pexp.nci60, pexp.ccle, housekeep.genes, "NCI60 (CV)", "BROAD (CV)")
#PlotCV(pexp.nci60, ccle.prot.mat, housekeep.genes, "NCI60 (CV)", "BROAD (CV)")
PlotCV(pexp.nci60, pexp.mghcc_breast, housekeep.genes, "NCI60 (CV)", "MGHCC_BREAST (CV)")
PlotCV(pexp.nci60, pexp.gdsc, housekeep.genes, "NCI60 (CV)", "GDSC (CV)")

PlotCV(pexp.mipb_hgsoc, pexp.ccle, housekeep.genes, "MIPB_HGSOC (CV)", "BROAD (CV)")
#PlotCV(pexp.mipb_hgsoc, ccle.prot.mat, housekeep.genes, "MIPB_HGSOC (CV)", "BROAD (CV)")

PlotCV(pexp.uw_tnbc, pexp.mghcc_breast, housekeep.genes, "UW_TNBC (CV)", "MGHCC_BREAST (CV)")
PlotCV(pexp.uw_tnbc, pexp.ccle, housekeep.genes, "UW_TNBC (CV)", "BROAD (CV)")
#PlotCV(pexp.uw_tnbc, ccle.prot.mat, housekeep.genes, "UW_TNBC (CV)", "BROAD (CV)")

PlotCV(pexp.mghcc_breast, pexp.ccle, housekeep.genes, "MGHCC_BREAST (CV)", "BROAD (CV)")
PlotCV(pexp.mghcc_breast, ccle.prot.mat, housekeep.genes, "MGHCC_BREAST (CV)", "BROAD (CV)")
#PlotCV(pexp.gdsc, ccle.prot.mat, housekeep.genes, "GDSC (CV)", "BROAD (CV)")
PlotCV(pexp.gdsc, pexp.ccle, housekeep.genes, "GDSC (CV)", "BROAD (CV)")




######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     TRANSCRIPTOME profiles                                   ###
######################################################################
# load correlation data
GEXP.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/gexp_cor_long.rds")
GEXP.COR.long <- GEXP.COR.long[GEXP.COR.long$ID == 1, ] #subset to identical cell lines

# add platform info
#GDSC, NCI60  = Array
#BROAD, UHN, gCSI, OHSU = RNAseq
GEXP.COR.long$DATA1_TYPE <- "RNASEQ"
GEXP.COR.long$DATA2_TYPE <- "RNASEQ"
GEXP.COR.long$DATA1_TYPE[grep("GDSC|NCI60", GEXP.COR.long$Dataset1)] <- "ARRAY"
GEXP.COR.long$DATA2_TYPE[grep("GDSC|NCI60", GEXP.COR.long$Dataset2)] <- "ARRAY"
GEXP.COR.long$DATATYPE_CLASS <- paste(GEXP.COR.long$DATA1_TYPE, GEXP.COR.long$DATA2_TYPE, sep="_" )
GEXP.COR.long$DATATYPE_SAME <- GEXP.COR.long$DATA1_TYPE == GEXP.COR.long$DATA2_TYPE
table(GEXP.COR.long$DATATYPE_CLASS )

#Statistical comparisons
arr_arr_IDx <- grep("ARRAY_ARRAY", GEXP.COR.long$DATATYPE_CLASS) 
rnaseq_rnseq_IDx <- grep("RNASEQ_RNASEQ", GEXP.COR.long$DATATYPE_CLASS) 
arr_rnaseq_IDx <- grep("ARRAY_RNASEQ|RNASEQ_ARRAY", GEXP.COR.long$DATATYPE_CLASS)

#Array/Array vs. RNAseq/RNAseq
wilcox.test(GEXP.COR.long$Cor[arr_arr_IDx],GEXP.COR.long$Cor[rnaseq_rnseq_IDx])$p.value
#RNAseq/RNAseq vs. Array/RNAseq
wilcox.test(GEXP.COR.long$Cor[rnaseq_rnseq_IDx],GEXP.COR.long$Cor[arr_rnaseq_IDx])$p.value
#Array/Array vs. Array/RNAseq
wilcox.test(GEXP.COR.long$Cor[arr_arr_IDx],GEXP.COR.long$Cor[arr_rnaseq_IDx])$p.value

#Make Corrplot
MakeCorrPlot(GEXP.COR.long, order= c("BROAD", "UHN", "OHSU", "gCSI", "GDSC", "NCI60") )


#Make Violin Plot
GEXP.COR.long$DATATYPE_NUM <- 1
GEXP.COR.long$DATATYPE_NUM[GEXP.COR.long$DATATYPE_CLASS == "ARRAY_RNASEQ"] <- 2
GEXP.COR.long$DATATYPE_NUM[GEXP.COR.long$DATATYPE_CLASS == "RNASEQ_ARRAY"] <- 2
GEXP.COR.long$DATATYPE_NUM[GEXP.COR.long$DATATYPE_CLASS == "RNASEQ_RNASEQ"] <- 3
p <- ggplot(GEXP.COR.long, aes(factor(DATATYPE_NUM), Cor))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p <- p + geom_violin(aes(fill = factor(DATATYPE_NUM)), scale = "width", trim=FALSE, colour = "gray30") 
cols <- c("1" = "powderblue", "2" = "seagreen3","3" = "lightslateblue")
p <- p + scale_fill_manual(values = cols)
p + stat_summary(fun = "mean", geom = "point", shape = 24, size = 3, color = "gray4",fill="white")



######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     DRUG SENSITIVITY profiles                                ###
######################################################################
# load correlation data
DSS.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/dss_cor_long.rds")
DSS.COR.long <- DSS.COR.long[DSS.COR.long$ID == 1, ] #subset to identical cell lines

# add platform info
#GDSC  = Fluorophore
#BROAD, UHN, gCSI, OHSU = CTG
DSS.COR.long$DATA1_TYPE <- "CTG"
DSS.COR.long$DATA2_TYPE <- "CTG"
DSS.COR.long$DATA1_TYPE[grep("GDSC", DSS.COR.long$Dataset1)] <- "FLUOR"
DSS.COR.long$DATA2_TYPE[grep("GDSC", DSS.COR.long$Dataset2)] <- "FLUOR"
DSS.COR.long$DATATYPE_CLASS <- paste(DSS.COR.long$DATA1_TYPE, DSS.COR.long$DATA2_TYPE, sep="_" )
DSS.COR.long$DATATYPE_SAME <- DSS.COR.long$DATA1_TYPE == DSS.COR.long$DATA2_TYPE
table(DSS.COR.long$DATATYPE_CLASS )

#Statistical comparisons
ctg_ctg_IDx <- grep("CTG_CTG", DSS.COR.long$DATATYPE_CLASS) 
ctg_fluor_IDx <- grep("FLUOR_CTG", DSS.COR.long$DATATYPE_CLASS) 

##CTG/CTG vs CTG/FLUOR
wilcox.test(DSS.COR.long$Cor[ctg_ctg_IDx],DSS.COR.long$Cor[ctg_fluor_IDx])$p.value

#Make density plot
alpha_value = 0.6
c1 <- adjustcolor("mediumorchid1", alpha.f=alpha_value) 
alpha_value = 1
c2 <- adjustcolor("mediumorchid4", alpha.f=alpha_value)

plot(density(DSS.COR.long$Cor[ctg_ctg_IDx]),  col=c2, las=1, main="", xlim=c(0,1), xlab="Correlation")
polygon(density(DSS.COR.long$Cor[ctg_ctg_IDx]), col=c2, border=c2)
lines(density(DSS.COR.long$Cor[ctg_fluor_IDx]), col=c1)
polygon(density(DSS.COR.long$Cor[ctg_fluor_IDx]), col=c1, border=c1)
legend(0.01, 3.8, legend=c("CTG/CTG", "Syto60/CTG"),
       pch=15, cex=1,col=c("mediumorchid4", "mediumorchid1"),
       box.col="white", title="", 
       bty="n", y.intersp=0.5, x.intersp=0.4)

#Make Corrplot
MakeCorrPlot(DSS.COR.long, order= c("BROAD_CCLE", "BROAD_CTRP", "FIMM", "gCSI", "GDSC") )



######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     TRAGET ADDICTION profiles                                ###
######################################################################
# load correlation data
TAS.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/tas_cor_long.rds")
TAS.COR.long <- TAS.COR.long[TAS.COR.long$ID == 1, ] #subset to identical cell lines

# add platform info
#GDSC  = Fluorophore
#BROAD, UHN, gCSI, OHSU = CTG
TAS.COR.long$DATA1_TYPE <- "CTG"
TAS.COR.long$DATA2_TYPE <- "CTG"
TAS.COR.long$DATA1_TYPE[grep("GDSC", TAS.COR.long$Dataset1)] <- "FLUOR"
TAS.COR.long$DATA2_TYPE[grep("GDSC", TAS.COR.long$Dataset2)] <- "FLUOR"
TAS.COR.long$DATATYPE_CLASS <- paste(TAS.COR.long$DATA1_TYPE, TAS.COR.long$DATA2_TYPE, sep="_" )
TAS.COR.long$DATATYPE_SAME <- TAS.COR.long$DATA1_TYPE == TAS.COR.long$DATA2_TYPE
table(TAS.COR.long$DATATYPE_CLASS )

#Statistical comparisons
ctg_ctg_IDx <- grep("CTG_CTG", TAS.COR.long$DATATYPE_CLASS) 
ctg_fluor_IDx <- grep("FLUOR_CTG", TAS.COR.long$DATATYPE_CLASS) 

##CTG/CTG vs CTG/FLUOR
wilcox.test(TAS.COR.long$Cor[ctg_ctg_IDx],TAS.COR.long$Cor[ctg_fluor_IDx])$p.value

#Make density plot
alpha_value = 0.3
c1 <- adjustcolor("red1", alpha.f=alpha_value) 
alpha_value = 0.8
c2 <- adjustcolor("red1", alpha.f=alpha_value)

plot(density(TAS.COR.long$Cor[ctg_ctg_IDx]),  col=c2, las=1, main="", xlim=c(0,1), ylim=c(0,2.5), xlab="Correlation")
polygon(density(TAS.COR.long$Cor[ctg_ctg_IDx]), col=c2, border=c2)
lines(density(TAS.COR.long$Cor[ctg_fluor_IDx]), col=c1)
polygon(density(TAS.COR.long$Cor[ctg_fluor_IDx]), col=c1, border=c1)
legend(0.01, 3.8, legend=c("CTG/CTG", "Syto60/CTG"),
       pch=15, cex=1,col=c("mediumorchid4", "mediumorchid1"),
       box.col="white", title="", 
       bty="n", y.intersp=0.5, x.intersp=0.4)

#Make Corrplot
MakeCorrPlot(TAS.COR.long, order= c("BROAD_CCLE", "BROAD_CTRP", "FIMM", "gCSI", "GDSC") )





######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     PHOSPHOPROTEOME profiles                                 ###
######################################################################
# load correlation data
PHOS.COR.long <- readRDS(file="clip-meta/processed_files/reprod_analysis/phos_cor_long.rds")
PHOS.COR.long <- PHOS.COR.long[PHOS.COR.long$ID == 1, ] #subset to identical cell lines

# add platform info
#GDSC  = MS-based 
#MCLP, UHN, BROAD, OHSU, NCI60 = RPPA
PHOS.COR.long$DATA1_TYPE <- "RPPA"
PHOS.COR.long$DATA2_TYPE <- "RPPA"
PHOS.COR.long$DATA1_TYPE[grep("GDSC", PHOS.COR.long$Dataset1)] <- "MS-based"
PHOS.COR.long$DATA2_TYPE[grep("GDSC", PHOS.COR.long$Dataset2)] <- "MS-based"
PHOS.COR.long$DATATYPE_CLASS <- paste(PHOS.COR.long$DATA1_TYPE, PHOS.COR.long$DATA2_TYPE, sep="_" )
PHOS.COR.long$DATATYPE_SAME <- PHOS.COR.long$DATA1_TYPE == PHOS.COR.long$DATA2_TYPE
table(PHOS.COR.long$DATATYPE_CLASS )

#Statistical comparisons
rppa_rppa_IDx <- grep("RPPA_RPPA", PHOS.COR.long$DATATYPE_CLASS) 
rppa_ms_IDx <- grep("RPPA_MS-based", PHOS.COR.long$DATATYPE_CLASS) 

##CTG/CTG vs CTG/FLUOR
wilcox.test(PHOS.COR.long$Cor[rppa_rppa_IDx],PHOS.COR.long$Cor[rppa_ms_IDx])$p.value

#Make density plot
alpha_value = 0.8
c1 <- adjustcolor("lightseagreen", alpha.f=alpha_value) 
alpha_value = 0.8
c2 <- adjustcolor("indianred1", alpha.f=alpha_value)

#Make density plot
plot(density(PHOS.COR.long$Cor[rppa_ms_IDx]),  col=c1, las=1, main="", xlim=c(0,1), xlab="Correlation")
polygon(density(PHOS.COR.long$Cor[rppa_ms_IDx]), col=c1, border=c1)
lines(density(PHOS.COR.long$Cor[rppa_rppa_IDx]), col=c2)
polygon(density(PHOS.COR.long$Cor[rppa_rppa_IDx]), col=c2, border=c2)
legend(0.5, 3.9, legend=c("MS-based/RPPA", "RPPA/RPPA"),
       pch=15, cex=1,col=c(c1, c2),
       box.col="white", title="", 
       bty="n", y.intersp=0.3, x.intersp=0.4)



######################################################################
#####     Assess reproducibility of PROTEOME profiles              ###
#####     Normalized and un-normalized                             ###
######################################################################


load("~/Desktop/FIMM_Work/CLIP_Datasets/PEXP.RData")


# Non-normalized CCLE data
ccle.prot <- read.csv("~/Downloads/summed_sn_non_normalized.csv", stringsAsFactors = F)

#ccle.prot <- readxl::read_xlsx("~/Downloads/Table_S2_Protein_Quant_Normalized.xlsx", sheet=2)
ccle.prot <- as.data.frame(ccle.prot, stringsAsFactors=F)
ccle.prot <- ccle.prot[, 1:426]
ccle.prot <- ccle.prot[, c(2,47:426)]

#Normalize by bridge
bridgeIDx <- grep("bridge", colnames(ccle.prot))
ccle.prot.norm <- c()
for (i in bridgeIDx){
  n2 = i-1
  n1 = i - 9
  ccle.prot1 <- ccle.prot[, c(1,n1:n2) ]
  ccle.prot1 <- ccle.prot1[order(ccle.prot1$Gene.Symbol), ]
  ccle.prot1 <- ccle.prot1[ccle.prot1$Gene.Symbol %in% rownames(pexp.ccle), ]
  ccle.prot1[,-1] <- apply(ccle.prot1[,-1],2,function(x) log2(x+1))
  colnames(ccle.prot1) <- sapply(strsplit(colnames(ccle.prot1), "\\_"), "[[", 1)
  ccle.prot1.list <- split(ccle.prot1[,-1], ccle.prot1$Gene)
  ccle.prot.mat <- lapply(ccle.prot1.list, function(x) apply(x,2, function(y) mean(y, na.rm=T)))
  ccle.prot.mat <- do.call(rbind, ccle.prot.mat)
  ccle.prot.norm <- cbind(ccle.prot.norm, ccle.prot.mat)
}
ccle.prot.norm <- as.data.frame(ccle.prot.norm, stringsAsFactors =F)


# Before
pdf("~/Desktop/FIMM_Work/CLIP_Review/PEXP_before.pdf", height = 6, width = 8)
par(mfrow=c(2,3))
plot(density(as.matrix(pexp.ccle), na.rm=T), main="CCLE", xlab="Bridge normalized", las=1)
plot(density(as.matrix(pexp.mghcc_breast), na.rm=T), main="MGHCC_BREAST", xlab="Bridge normalized", las=1)
plot(density(as.matrix(pexp.gdsc), na.rm=T), main="GDSC", cex.lab=0.7,las=1,
     xlab="Sum of column-normalized (gene) \nTMT spectrum intensities \nfollowed by row-mean scaling")
plot(density(as.matrix(pexp.uw_tnbc), na.rm=T), main="UW_TNBC", xlab="iBAQ", las=1)
plot(density(as.matrix(pexp.nci60), na.rm=T), main="NCI60", xlab="LFQ", las=1)
plot(density(as.matrix(pexp.mipb_hgsoc), na.rm=T), main="MIPB_HGSOC", xlab="LFQ", las=1)
dev.off()

# Normalize datasets
#After
pdf("~/Desktop/FIMM_Work/CLIP_Review/PEXP_after.pdf", height = 6, width = 8)
par(mfrow=c(2,3))
plot(density(as.matrix(ccle.prot.norm), na.rm=T), main="CCLE" , xlab="Non-normalized", las=1)
pexp.mghcc_breast1 <- apply(pexp.mghcc_breast,2,log2)
plot(density(as.matrix(pexp.mghcc_breast1), na.rm=T), main="MGHCC_BREAST", xlab="Bridge normalized + log2 transform", las=1)
pexp.uw_tnbc1 <- apply(pexp.uw_tnbc,2,log2)
plot(density(as.matrix(pexp.uw_tnbc1), na.rm=T), main="UW_TNBC", xlab="iBAQ + log2 transform", las=1)
pexp.nci601 <- apply(pexp.nci60,2,log2)
plot(density(as.matrix(pexp.nci601), na.rm=T), main="NCI60", xlab="LFQ + log2 transform", las=1)
dev.off()


#Assess reproducibility
#TMT labeled vs. TMT labeled
PEXP.COR.BROAD.MGHCC_BREAST <- CompareTwoStudies_Rcorr(ccle.prot.norm, pexp.mghcc_breast1, "spearman")
PEXP.COR.GDSC.MGHCC_BREAST   <- CompareTwoStudies_Rcorr(pexp.gdsc, pexp.mghcc_breast1, "spearman")
PEXP.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(ccle.prot.norm, pexp.gdsc, "spearman")

# Non-labeled vs. Non-labeled
PEXP.COR.UW_TNBC.NCI60 <- CompareTwoStudies_Rcorr(pexp.uw_tnbc1, pexp.nci601, "spearman")
PEXP.COR.NCI60.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.nci601, pexp.mipb_hgsoc, "spearman")

# TMT labeled  vs. Non-labeled
PEXP.COR.UW_TNBC.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.uw_tnbc1, pexp.mghcc_breast1, "spearman")
PEXP.COR.NCI60.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.nci601, pexp.mghcc_breast1, "spearman")
PEXP.COR.UW_TNBC.BROAD <- CompareTwoStudies_Rcorr(pexp.uw_tnbc1, ccle.prot.norm, "spearman")
PEXP.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(pexp.nci601, pexp.gdsc, "spearman")
PEXP.COR.NCI60.BROAD <- CompareTwoStudies_Rcorr(pexp.nci601, ccle.prot.norm, "spearman")
PEXP.COR.BROAD.MPIB_HGSOC <- CompareTwoStudies_Rcorr(ccle.prot.norm, pexp.mipb_hgsoc, "spearman")
PEXP.COR.GDSC.MPIB_HGSOC  <- CompareTwoStudies_Rcorr(pexp.gdsc, pexp.mipb_hgsoc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("PEXP.COR", allvars)]
PEXP.COR.long <- CollapseCorMatrices(corvars)




#### Corrplot GEXP vs. PEXP ######

#
GEXP.PEXP.COR.CCLE.UC_TNBC <- CompareTwoStudies_Rcorr(ccle.rnaseq, TNBC.prot, "spearman")
GEXP.PEXP.COR.CCLE.CCLE <- CompareTwoStudies_Rcorr(ccle.rnaseq, ccle.prot.mat, "spearman")
GEXP.PEXP.COR.CCLE.GDSC <- CompareTwoStudies_Rcorr(ccle.rnaseq, COLREC.prot, "spearman")
GEXP.PEXP.COR.CCLE.MGHCC_BREAST <- CompareTwoStudies_Rcorr(ccle.rnaseq, HAAS.BRCA.prot, "spearman")
GEXP.PEXP.COR.CCLE.NCI60 <- CompareTwoStudies_Rcorr(ccle.rnaseq, NCI60.prot, "spearman")
GEXP.PEXP.COR.CCLE.MPIB_HGSOC <- CompareTwoStudies_Rcorr(ccle.rnaseq, OVCA.prot, "spearman")

GEXP.PEXP.COR.GDSC.UC_TNBC <- CompareTwoStudies_Rcorr(gdsc.array, TNBC.prot, "spearman")
GEXP.PEXP.COR.GDSC.CCLE <- CompareTwoStudies_Rcorr(gdsc.array, ccle.prot.mat, "spearman")
GEXP.PEXP.COR.GDSC.GDSC <- CompareTwoStudies_Rcorr(gdsc.array, COLREC.prot, "spearman")
GEXP.PEXP.COR.GDSC.MGHCC_BREAST <- CompareTwoStudies_Rcorr(gdsc.array, HAAS.BRCA.prot, "spearman")
GEXP.PEXP.COR.GDSC.NCI60 <- CompareTwoStudies_Rcorr(gdsc.array, NCI60.prot, "spearman")
GEXP.PEXP.COR.GDSC.MPIB_HGSOC <- CompareTwoStudies_Rcorr(gdsc.array, OVCA.prot, "spearman")

GEXP.PEXP.COR.NCI60.UC_TNBC <- CompareTwoStudies_Rcorr(nci60.array, TNBC.prot, "spearman")
GEXP.PEXP.COR.NCI60.CCLE <- CompareTwoStudies_Rcorr(nci60.array, ccle.prot.mat, "spearman")
GEXP.PEXP.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(nci60.array, COLREC.prot, "spearman")
GEXP.PEXP.COR.NCI60.MGHCC_BREAST <- CompareTwoStudies_Rcorr(nci60.array, HAAS.BRCA.prot, "spearman")
GEXP.PEXP.COR.NCI60.NCI60 <- CompareTwoStudies_Rcorr(nci60.array, NCI60.prot, "spearman")
GEXP.PEXP.COR.NCI60.MPIB_HGSOC <- CompareTwoStudies_Rcorr(nci60.array, OVCA.prot, "spearman")


GEXP.PEXP.COR.UHN.UC_TNBC <- CompareTwoStudies_Rcorr(bfg.rnaseq, TNBC.prot, "spearman")
GEXP.PEXP.COR.UHN.CCLE <- CompareTwoStudies_Rcorr(bfg.rnaseq, ccle.prot.mat, "spearman")
GEXP.PEXP.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(bfg.rnaseq, COLREC.prot, "spearman")
GEXP.PEXP.COR.UHN.MGHCC_BREAST <- CompareTwoStudies_Rcorr(bfg.rnaseq, HAAS.BRCA.prot, "spearman")
GEXP.PEXP.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(bfg.rnaseq, NCI60.prot, "spearman")
GEXP.PEXP.COR.UHN.MPIB_HGSOC <- CompareTwoStudies_Rcorr(bfg.rnaseq, OVCA.prot, "spearman")


GEXP.PEXP.COR.OHSU.UC_TNBC <- CompareTwoStudies_Rcorr(gray.rnaseq, TNBC.prot, "spearman")
GEXP.PEXP.COR.OHSU.CCLE <- CompareTwoStudies_Rcorr(gray.rnaseq, ccle.prot.mat, "spearman")
GEXP.PEXP.COR.OHSU.GDSC <- CompareTwoStudies_Rcorr(gray.rnaseq, COLREC.prot, "spearman")
GEXP.PEXP.COR.OHSU.MGHCC_BREAST <- CompareTwoStudies_Rcorr(gray.rnaseq, HAAS.BRCA.prot, "spearman")
GEXP.PEXP.COR.OHSU.NCI60 <- CompareTwoStudies_Rcorr(gray.rnaseq, NCI60.prot, "spearman")
GEXP.PEXP.COR.OHSU.MPIB_HGSOC <- CompareTwoStudies_Rcorr(gray.rnaseq, OVCA.prot, "spearman")


GEXP.PEXP.COR.gCSI.UC_TNBC <- CompareTwoStudies_Rcorr(klijn.rnaseq, TNBC.prot, "spearman")
GEXP.PEXP.COR.gCSI.CCLE <- CompareTwoStudies_Rcorr(klijn.rnaseq, ccle.prot.mat, "spearman")
GEXP.PEXP.COR.gCSI.GDSC <- CompareTwoStudies_Rcorr(klijn.rnaseq, COLREC.prot, "spearman")
GEXP.PEXP.COR.gCSI.MGHCC_BREAST <- CompareTwoStudies_Rcorr(klijn.rnaseq, HAAS.BRCA.prot, "spearman")
GEXP.PEXP.COR.gCSI.NCI60 <- CompareTwoStudies_Rcorr(klijn.rnaseq, NCI60.prot, "spearman")
GEXP.PEXP.COR.gCSI.MPIB_HGSOC <- CompareTwoStudies_Rcorr(klijn.rnaseq, OVCA.prot, "spearman")


#Process all cor matrices
x <- ls()
x <- x[grep("GEXP.PEXP.COR", x)]

compare.views.cor.list.long <- list()
for (n in x){
  n.df <- get(n)
  compare.views.cor.list.long[[n]] <- n.df
}

compare.views.cor.list.filtered <- list()
for (n in 1:length(compare.views.cor.list.long)){
  n.df <- compare.views.cor.list.long[[n]]
  n.df$ID <- 0
  n.df$ID[which(n.df$Cell1 == n.df$Cell2)] <- 1
  n.df <- n.df[n.df$P < 0.05, ]
  #n.df <- n.df[n.df$n >= 10, ]
  if(nrow(n.df) == 0) next()
  listName <- names(compare.views.cor.list.long)[n]
  compare.views.cor.list.filtered[[listName]] <- cbind(n.df,listName)
}
do.call(rbind, lapply(compare.views.cor.list.filtered, dim))
compare.views.cor.list.filtered <- do.call(rbind, compare.views.cor.list.filtered)
compare.views.cor.list.filtered <- as.data.frame(compare.views.cor.list.filtered, stringsAsFactors = F)
colnames(compare.views.cor.list.filtered)[7] <- "Comparison"
compare.views.cor.list.filtered$Comparison <- as.character(compare.views.cor.list.filtered$Comparison)
compare.views.cor.list.filtered$Cell1 <- as.character(compare.views.cor.list.filtered$Cell1)
compare.views.cor.list.filtered$Cell2 <- as.character(compare.views.cor.list.filtered$Cell2)
#compare.views.cor.list.filtered$Datatype <-  sapply(strsplit(compare.views.cor.list.filtered$Comparison, "\\."), "[[", 1)
compare.views.cor.list.filtered$Dataset1 <-  sapply(strsplit(compare.views.cor.list.filtered$Comparison, "\\."), "[[", 4)
compare.views.cor.list.filtered$Dataset2 <-  sapply(strsplit(compare.views.cor.list.filtered$Comparison, "\\."), "[[", 5)
rownames(compare.views.cor.list.filtered) <- NULL
compare.views.cor.list.filtered$Cor <- as.numeric(compare.views.cor.list.filtered$Cor)
compare.views.cor.list.filtered$P <- as.numeric(compare.views.cor.list.filtered$P)
compare.views.cor.list.filtered <- compare.views.cor.list.filtered[!is.na(compare.views.cor.list.filtered$Cor), ]
save(file="~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis//GEXP_PEXP_Correlation.RData", 
     compare.views.cor.list.long, compare.views.cor.list.filtered,compress="bzip2")

load("~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis//GEXP_PEXP_Correlation.RData")

#
#Subset Identical and non-identical
compare.cor.views.ID <- compare.views.cor.list.filtered[compare.views.cor.list.filtered$ID == 1, ]
compare.cor.views.nonID <- compare.views.cor.list.filtered[compare.views.cor.list.filtered$ID == 0, ]


Cor.GEXP.PEXP.COMPARE <- compare.cor.views.ID
Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS <- paste(Cor.GEXP.PEXP.COMPARE$Dataset1, Cor.GEXP.PEXP.COMPARE$Dataset2, sep="_vs_" )
Cor.GEXP.PEXP.COMPARE$DATA1_TYPE <- "RNAseq"
Cor.GEXP.PEXP.COMPARE$DATA2_TYPE <- "NL"
Cor.GEXP.PEXP.COMPARE$DATA1_TYPE[grep("GDSC|NCI60", Cor.GEXP.PEXP.COMPARE$Dataset1)] <- "Array"
Cor.GEXP.PEXP.COMPARE$DATA2_TYPE[grep("GDSC|CCLE|MGHCC_BREAST", Cor.GEXP.PEXP.COMPARE$Dataset2)] <- "TMT"
Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS <- paste(Cor.GEXP.PEXP.COMPARE$DATA1_TYPE, Cor.GEXP.PEXP.COMPARE$DATA2_TYPE, sep="_" )
table(Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS )
summary(Cor.GEXP.PEXP.COMPARE$Cor~Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS)
summary(Cor.GEXP.PEXP.COMPARE$Cor~Cor.GEXP.PEXP.COMPARE$Dataset2)

#Make density plot
p1 <- grep("Array_NL", Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS) #Array vs NL
p2 <- grep("Array_TMT", Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS) #Array vs TMT
p3 <- grep("RNAseq_NL", Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS) #RNAseq vs NL
p4 <- grep("RNAseq_TMT", Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS) #RNAseq vs TMT

pcor1 <- Cor.GEXP.PEXP.COMPARE$Cor[p1]
pcor2 <- Cor.GEXP.PEXP.COMPARE$Cor[p2]
pcor3 <- Cor.GEXP.PEXP.COMPARE$Cor[p3]
pcor4 <- Cor.GEXP.PEXP.COMPARE$Cor[p4]
wilcox.test(pcor1,pcor2)
wilcox.test(pcor1,pcor2)$p.value
#p-value =3.114575e-66
t.test(pcor1,pcor2)
#Mean 0.4023826 0.2453164

wilcox.test(pcor1,pcor3)
wilcox.test(pcor1,pcor3)$p.value
#p-value =8.372372e-53
t.test(pcor1,pcor3)
#Mean 0.4023826 0.5425047 

wilcox.test(pcor1,pcor4)
wilcox.test(pcor1,pcor4)$p.value
#p-value =2.45304e-47
t.test(pcor1,pcor4)
#Mean 0.4023826 0.3116343 

wilcox.test(pcor2,pcor3)
wilcox.test(pcor2,pcor3)$p.value
#p-value = 3.822005e-89
t.test(pcor2,pcor3)
#Mean 0.2453164 0.5425047 

wilcox.test(pcor2,pcor4)
wilcox.test(pcor2,pcor4)$p.value
#p-value = 8.850125e-46
t.test(pcor2,pcor4)
#Mean 0.2453164 0.3116343

wilcox.test(pcor3,pcor4)
wilcox.test(pcor3,pcor4)$p.value
#p-value = 4.119578e-105
t.test(pcor3,pcor4)
#Mean 0.5425047 0.3116343 

Cor.GEXP.PEXP.COMPARE$DATATYPE_NUM <- 1
Cor.GEXP.PEXP.COMPARE$DATATYPE_NUM[Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS == "RNAseq_TMT"] <- 2
Cor.GEXP.PEXP.COMPARE$DATATYPE_NUM[Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS == "Array_NL"] <- 3
Cor.GEXP.PEXP.COMPARE$DATATYPE_NUM[Cor.GEXP.PEXP.COMPARE$DATATYPE_CLASS == "Array_TMT"] <- 4
p <- ggplot(Cor.GEXP.PEXP.COMPARE, aes(factor(DATATYPE_NUM), Cor))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p <- p + geom_violin(aes(fill = factor(DATATYPE_NUM)), scale = "width", trim=FALSE, colour = "gray30") 
cols <- c("1" = "powderblue", "2" = "seagreen3",
          "3" = "lightslateblue", "4" = "mediumpurple4")
p <- p + scale_fill_manual(values = cols)
#p <-  p + geom_jitter(height = 0, width = 0.3, alpha=0.2, shape=19)
p + stat_summary(fun.y = "mean", geom = "point", shape = 24, size = 3, color = "gray4",fill="white")
quartz.save("/Volumes/Alok_New/FIMM/Breast_Genomics_Analysis_Tmp//Analysis/For_Paper/ViolinPlot_Correlation_GEXP_PEXP.pdf", 
            type = "pdf", device = dev.cur(), dpi = 100, width=6, height=4)

quartz.save("~/Desktop/FIMM_Work/Misc_Projects/Breast_Genomics_Analysis//Analysis/For_Paper/ViolinPlot_Correlation_GEXP_PEXP.pdf", 
            type = "pdf", device = dev.cur(), dpi = 100, width=6, height=4)


#Corrplot
Cor.GEXP.PEXP.COMPARE.MAT <- dcast(Cor.GEXP.PEXP.COMPARE[,c(8,9,3)],Dataset1~Dataset2, fill=0, fun.aggregate = mean)
Cor.GEXP.PEXP.COMPARE.MAT[,-1] <- apply(Cor.GEXP.PEXP.COMPARE.MAT[,-1],2, function(x){
  x[x==0] = NA
  return(x)})

rownames(Cor.GEXP.PEXP.COMPARE.MAT) <- Cor.GEXP.PEXP.COMPARE.MAT$Dataset1
Cor.GEXP.PEXP.COMPARE.MAT <- Cor.GEXP.PEXP.COMPARE.MAT[,-1]
M <- as.matrix(Cor.GEXP.PEXP.COMPARE.MAT)
colnames(M)[1] =  c("BROAD")
rownames(M)[1] =  c("BROAD")
M <- M[c(1,2,5,6,3,4), c(1,2,3,6,5,4)]
r <- range(M, na.rm = T)
bk1 = seq(0,0.1,length=10)
bk2 = seq(0.2,0.7,length=30)
hmcols1 <- colorRampPalette(c("gray"))(length(bk1))
hmcols2 <- colorRampPalette(c("white","red"))(length(bk2))
bk <- c(bk1, bk2)
hmcols <- c(hmcols1,hmcols2)
pheatmap(M, cluster_rows = F, cluster_cols = F, 
         display_numbers=T, fontsize_number=12,
         border_color = NA, color=hmcols, breaks=bk)
quartz.save("/Volumes/Alok_New/FIMM/Breast_Genomics_Analysis_Tmp//Analysis/For_Paper/CorrPlot_Correlation_GEXP.PEXP_HM.pdf", 
            type = "pdf", device = dev.cur(), dpi = 100, width=4, height=4)





