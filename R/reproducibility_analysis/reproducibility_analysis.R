##### load functions and packages #####
source("clip-meta/R/reproducibility_analysis/reproducibility_analysis_functions.R")

##### load modality-specific pre-processed data ######
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


######################################################################
#####     Assess reproducibility of METHYLATION profiles           ###
######################################################################

# load Methylation profiles
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


# Compute correlation
METH.COR.NCI60.OHSU <- CompareTwoStudies_Rcorr(meth.nci60, meth.ohsu, "spearman")
METH.COR.NCI60.GDSC  <- CompareTwoStudies_Rcorr(meth.nci60, meth.gdsc, "spearman")
METH.COR.NCI60.BROAD  <- CompareTwoStudies_Rcorr(meth.nci60, meth.broad, "spearman")
METH.COR.OHSU.GDSC  <- CompareTwoStudies_Rcorr(meth.ohsu, meth.sanger, "spearman")
METH.COR.OHSU.BROAD  <- CompareTwoStudies_Rcorr(meth.ohsu, meth.broad, "spearman")
METH.COR.GDSC.BROAD  <- CompareTwoStudies_Rcorr(meth.gdsc, meth.broad, "spearman")

# save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("METH.COR", allvars)]
METH.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/meth_cor_long.rds", METH.COR.long,compress="bzip2")

######################################################################
#####     Assess reproducibility of MUTATION profiles              ###
######################################################################


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


#
# Matthew's Correlation coefficient = Pearson correlation coefficient estimated for two binary variables
MUT.COR.GDSC.CCLE <- CompareTwoStudies_Rcorr(mut.sanger, mut.broad, "pearson")
MUT.COR.GDSC.NCI60 <- CompareTwoStudies_Rcorr(mut.sanger, mut.nci60, "pearson")
MUT.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(mut.sanger, mut.ohsu, "pearson")
MUT.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(mut.sanger, mut.gcsi, "pearson")

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


CNV.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.sanger, "spearman")
CNV.COR.UHN.CCLE <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.broad, "spearman")
CNV.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.nci60, "spearman")
CNV.COR.UHN.OHSU <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.ohsu, "spearman")
CNV.COR.UHN.gCSI <- CompareTwoStudies_Rcorr(cnv.uhn, cnv.gcsi, "spearman")

CNV.COR.CCLE.GDSC <- CompareTwoStudies_Rcorr(cnv.broad, cnv.sanger, "spearman")
CNV.COR.CCLE.NCI60 <- CompareTwoStudies_Rcorr(cnv.broad, cnv.nci60, "spearman")
CNV.COR.CCLE.OHSU <- CompareTwoStudies_Rcorr(cnv.broad, cnv.ohsu, "spearman")
CNV.COR.CCLE.gCSI <- CompareTwoStudies_Rcorr(cnv.broad, cnv.gcsi, "spearman")

CNV.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(cnv.nci60, cnv.sanger, "spearman")
CNV.COR.NCI60.OHSU <- CompareTwoStudies_Rcorr(cnv.nci60, cnv.ohsu, "spearman")
CNV.COR.NCI60.gCSI <- CompareTwoStudies_Rcorr(cnv.nci60, cnv.gcsi, "spearman")

CNV.COR.OHSU.GDSC <- CompareTwoStudies_Rcorr(cnv.ohsu, cnv.sanger, "spearman")
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

#RNAseq vs. RNAseq
GEXP.COR.UHN.OHSU <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.ohsu, "spearman")
GEXP.COR.UHN.BROAD <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.broad, "spearman")
GEXP.COR.UHN.gCSI <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.gcsi, "spearman")
GEXP.COR.BROAD.OHSU <- CompareTwoStudies_Rcorr(gexp.broad, gexp.ohsu, "spearman")
GEXP.COR.BROAD.gCSI <- CompareTwoStudies_Rcorr(gexp.broad, gexp.gcsi, "spearman")
GEXP.COR.OHSU.gCSI <- CompareTwoStudies_Rcorr(gexp.ohsu, gexp.gcsi, "spearman")

#RNAseq vs. Array
GEXP.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.sanger, "spearman")
GEXP.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(gexp.uhn, gexp.nci60, "spearman")
GEXP.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(gexp.broad, gexp.gdsc, "spearman")
GEXP.COR.BROAD.NCI60 <- CompareTwoStudies_Rcorr(gexp.broad, gexp.nci60, "spearman")
GEXP.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(gexp.sanger, gexp.ohsu, "spearman")
GEXP.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(gexp.sanger, gexp.gcsi, "spearman")
GEXP.COR.OHSU.NCI60 <- CompareTwoStudies_Rcorr(gexp.ohsu, gexp.nci60, "spearman")
GEXP.COR.NCI60.gCSI <- CompareTwoStudies_Rcorr(gexp.nci60, gexp.gcsi, "spearman")

#Array vs. Array
GEXP.COR.GDSC.NCI60 <- CompareTwoStudies_Rcorr(gexp.sanger, gexp.nci60, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("GEXP.COR", allvars)]
GEXP.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/gexp_cor_long.rds", GEXP.COR.long,compress="bzip2")



######################################################################
#####     Assess reproducibility of PROTEOME profiles              ###
######################################################################

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

# TMT labeled vs. TMT labeled
PEXP.COR.BROAD.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.broad, pexp.mghcc_breast, "spearman")
PEXP.COR.GDSC.MGHCC_BREAST   <- CompareTwoStudies_Rcorr(pexp.sanger, pexp.mghcc_breast, "spearman")
PEXP.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(pexp.broad, pexp.sanger, "spearman")

# Non-labeled vs. Non-labeled
PEXP.COR.UW_TNBC.NCI60 <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.nci60, "spearman")
PEXP.COR.NCI60.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.mpib_hgsoc, "spearman")

# TMT labeled  vs. Non-labeled
PEXP.COR.UW_TNBC.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.mghcc_breast, "spearman")
PEXP.COR.NCI60.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.mghcc_breast, "spearman")
PEXP.COR.UW_TNBC.BROAD <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.broad, "spearman")
PEXP.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.sanger, "spearman")
PEXP.COR.NCI60.BROAD <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.broad, "spearman")
PEXP.COR.BROAD.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.broad, pexp.mpib_hgsoc, "spearman")
PEXP.COR.GDSC.MPIB_HGSOC  <- CompareTwoStudies_Rcorr(pexp.sanger, pexp.mpib_hgsoc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("PEXP.COR", allvars)]
PEXP.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/pexp_cor_long.rds", PEXP.COR.long,compress="bzip2")

######################################################################
#####     Assess reproducibility of PHOSPHOPROTEOME profiles       ###
######################################################################
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


# RPPA vs. RPPA
PHOS.COR.UHN.OHSU <- CompareTwoStudies_Rcorr(phos.uhn, phos.ohsu, "spearman")
PHOS.COR.UHN.MCLP <- CompareTwoStudies_Rcorr(phos.uhn, phos.mclp, "spearman")
PHOS.COR.UHN.NCI60 <- CompareTwoStudies_Rcorr(phos.uhn, phos.nci60, "spearman")
PHOS.COR.UHN.BROAD <- CompareTwoStudies_Rcorr(phos.uhn, phos.ccle, "spearman")
PHOS.COR.OHSU.MCLP <- CompareTwoStudies_Rcorr(phos.ohsu, phos.mclp, "spearman")
PHOS.COR.OHSU.NCI60 <- CompareTwoStudies_Rcorr(phos.ohsu, phos.nci60, "spearman")
PHOS.COR.OHSU.BROAD <- CompareTwoStudies_Rcorr(phos.ohsu, phos.broad, "spearman")
PHOS.COR.MCLP.NCI60 <- CompareTwoStudies_Rcorr(phos.mclp, phos.nci60, "spearman")
PHOS.COR.MCLP.BROAD <- CompareTwoStudies_Rcorr(phos.mclp, phos.broad, "spearman")
PHOS.COR.BROAD.NCI60 <- CompareTwoStudies_Rcorr(phos.broad, phos.nci60, "spearman")

# RPPA vs. MS-based
PHOS.COR.MCLP.GDSC <- CompareTwoStudies_Rcorr(phos.mclp, phos.sanger, "spearman")
PHOS.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(phos.broad, phos.sanger, "spearman")
PHOS.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(phos.nci60, phos.sanger, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("PHOS.COR", allvars)]
PHOS.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/phos_cor_long.rds", PHOS.COR.long,compress="bzip2")



######################################################################
#####     Assess reproducibility of FUNCTIONAL profiles            ###
######################################################################

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


# RNAi vs. RNAi
FUNC.COR.DRIVE.ACHILLES <- CompareTwoStudies_Rcorr(func.drive, func.broad_achilles, "spearman")
FUNC.COR.DRIVE.UHN <- CompareTwoStudies_Rcorr(func.drive, func.uhn, "spearman")
FUNC.COR.ACHILLES.UHN <- CompareTwoStudies_Rcorr(func.broad_achilles, func.uhn, "spearman")

# CRISPR vs. CRISPR
FUNC.COR.BROAD_AVANA.BROAD_GECKO <- CompareTwoStudies_Rcorr(func.broad_avana, func.broad_gecko, "spearman")
FUNC.COR.BROAD_AVANA.BROAD_AML <- CompareTwoStudies_Rcorr(func.broad_avana, func.broad_aml, "spearman")
FUNC.COR.BROAD_AVANA.GDSC <- CompareTwoStudies_Rcorr(func.broad_avana, func.sanger, "spearman")
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
FUNC.COR.ACHILLES.GDSC <- CompareTwoStudies_Rcorr(func.broad_achilles, func.sanger, "spearman")
FUNC.COR.UHN.GDSC <- CompareTwoStudies_Rcorr(func.uhn, func.sanger, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("FUNC.COR", allvars)]
FUNC.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/func_cor_long.rds", FUNC.COR.long,compress="bzip2")



######################################################################
#####     Assess reproducibility of DRUG SENSITIVITY profiles      ###
######################################################################

# load Drug sensitivity data
path = "clip-meta/data/"
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
DSS.COR.GDSC.BROAD_CCLE <- CompareTwoStudies_Rcorr(dss.sanger, dss.broad_ccle, "spearman")
DSS.COR.GDSC.BROAD_CTRP <- CompareTwoStudies_Rcorr(dss.sanger, dss.broad_ctrp, "spearman")
DSS.COR.GDSC.FIMM <- CompareTwoStudies_Rcorr(dss.sanger, dss.fimm, "spearman")
DSS.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(dss.sanger, dss.gcsi, "spearman")
DSS.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(dss.sanger, dss.ohsu, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("DSS.COR", allvars)]
DSS.COR.long <- CollapseCorMatrices(corvars)
saveRDS(file="clip-meta/processed_files/reprod_analysis/dss_cor_long.rds", DSS.COR.long, compress="bzip2")


######################################################################
#####     Assess reproducibility of TARGET ADDICTION profiles      ###
######################################################################
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
TAS.COR.GDSC.BROAD_CCLE <- CompareTwoStudies_Rcorr(tas.sanger, tas.broad_ccle, "spearman")
TAS.COR.GDSC.BROAD_CTRP <- CompareTwoStudies_Rcorr(tas.sanger, tas.broad_ctrp, "spearman")
TAS.COR.GDSC.FIMM <- CompareTwoStudies_Rcorr(tas.sanger, tas.fimm, "spearman")
TAS.COR.GDSC.gCSI <- CompareTwoStudies_Rcorr(tas.sanger, tas.gcsi, "spearman")
TAS.COR.GDSC.OHSU <- CompareTwoStudies_Rcorr(tas.sanger, tas.ohsu, "spearman")


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
MakeCorrPlot(PEXP.COR.long, order= c("BROAD", "MGHCC_BREAST", "GDSC", "NCI60", "UW_TNBC", "MPIB_HGSOC") )


######################################################################
#####     Reproduciblity measures by platform/technique            ###
#####     PROTEOME profiles - Coefficient of variation CCV)        ###
######################################################################
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

# CV correlation between UW_TNBC and HAAS.BRCA
# Read House keeping genes table
# Obtained from Publication: https://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00089-9#secsect0045
housekeep.genes <- read.table("clip-meta/data/housekeeping-genes.txt", sep="\t", header = T, stringsAsFactors = F)

#
PlotCV(pexp.uw_tnbc, pexp.nci60, housekeep.genes, "UW_TNBC (CV)", "NCI60 (CV)")
PlotCV(pexp.mpib_hgsoc, pexp.nci60, housekeep.genes, "MPIB_HGSOC (CV)", "NCI60 (CV)")

PlotCV(pexp.nci60, pexp.broad, housekeep.genes, "NCI60 (CV)", "BROAD (CV)")
PlotCV(pexp.nci60, pexp.mghcc_breast, housekeep.genes, "NCI60 (CV)", "MGHCC_BREAST (CV)")
PlotCV(pexp.nci60, pexp.sanger, housekeep.genes, "NCI60 (CV)", "GDSC (CV)")

PlotCV(pexp.uw_tnbc, pexp.mghcc_breast, housekeep.genes, "UW_TNBC (CV)", "MGHCC_BREAST (CV)")
PlotCV(pexp.mpib_hgsoc, pexp.broad, housekeep.genes, "MPIB_HGSOC (CV)", "BROAD (CV)")
PlotCV(pexp.uw_tnbc, pexp.broad, housekeep.genes, "UW_TNBC (CV)", "BROAD (CV)")

PlotCV(pexp.mghcc_breast, pexp.broad, housekeep.genes, "MGHCC_BREAST (CV)", "BROAD (CV)")
PlotCV(pexp.sanger, pexp.broad, housekeep.genes, "GDSC (CV)", "BROAD (CV)")




######################################################################
#####     Assess reproducibility of PROTEOME profiles              ###
#####     Non-normalized BROAD PEXP profiles                       ###
######################################################################


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


# read Non-normalized BROAD PEXP data
# Downloaded from  https://gygi.med.harvard.edu/publications/ccle
PEXP_NONORM.list <- h5dump("clip-meta/data/PEXP_BROAD_NONNORMALIZED.h5",load=TRUE)
pexp.broad.nonnorm <- PEXP_NONORM.list$BROAD$mat
colnames(pexp.broad.nonnorm) <- PEXP_NONORM.list$BROAD$colnames
pexp.broad.nonnorm <- as.data.frame(pexp.broad.nonnorm, stringsAsFactors = F)
pexp.broad.nonnorm[,-1] <- apply(pexp.broad.nonnorm[,-1],2,as.numeric)

#Remove bridge columns 
bridgeIDx <- grep("bridge", colnames(pexp.broad.nonnorm))
pexp.broad.nonnorm <- pexp.broad.nonnorm[, -bridgeIDx]
pexp.broad.nonnorm <- pexp.broad.nonnorm[order(pexp.broad.nonnorm$Gene.Symbol), ]

#log transform
pexp.broad.nonnorm[,-1] <- apply(pexp.broad.nonnorm[,-1],2,function(x) log2(x+1))
colnames(pexp.broad.nonnorm) <- sapply(strsplit(colnames(pexp.broad.nonnorm), "\\_"), "[[", 1)

#Average proteins with multiple peptides
pexp.broad.nonnorm.list <- split(pexp.broad.nonnorm[,-1], pexp.broad.nonnorm$Gene)
pexp.broad.nonnorm.list <- lapply(pexp.broad.nonnorm.list, function(x) apply(x,2, function(y) mean(y, na.rm=T)))
pexp.broad.nonnorm <- do.call(rbind, pexp.broad.nonnorm.list)
rownames(pexp.broad.nonnorm) <- gsub("#", "", rownames(pexp.broad.nonnorm))
pexp.broad.nonnorm <- as.data.frame(pexp.broad.nonnorm, stringsAsFactors = F)

#Assess reproducibility with before normalize data
#TMT labeled vs. TMT labeled
PEXP.COR.BROAD.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.broad.nonnorm, pexp.mghcc_breast, "spearman")
PEXP.COR.GDSC.MGHCC_BREAST   <- CompareTwoStudies_Rcorr(pexp.sanger, pexp.mghcc_breast, "spearman")
PEXP.COR.BROAD.GDSC <- CompareTwoStudies_Rcorr(pexp.broad.nonnorm, pexp.sanger, "spearman")

# Non-labeled vs. Non-labeled
PEXP.COR.UW_TNBC.NCI60 <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.nci60, "spearman")
PEXP.COR.NCI60.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.mpib_hgsoc, "spearman")

# TMT labeled  vs. Non-labeled
PEXP.COR.UW_TNBC.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.mghcc_breast, "spearman")
PEXP.COR.NCI60.MGHCC_BREAST <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.mghcc_breast, "spearman")
PEXP.COR.UW_TNBC.BROAD <- CompareTwoStudies_Rcorr(pexp.uw_tnbc, pexp.broad.nonnorm, "spearman")
PEXP.COR.NCI60.GDSC <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.sanger, "spearman")
PEXP.COR.NCI60.BROAD <- CompareTwoStudies_Rcorr(pexp.nci60, pexp.broad.nonnorm, "spearman")
PEXP.COR.BROAD.MPIB_HGSOC <- CompareTwoStudies_Rcorr(pexp.broad.nonnorm, pexp.mpib_hgsoc, "spearman")
PEXP.COR.GDSC.MPIB_HGSOC  <- CompareTwoStudies_Rcorr(pexp.sanger, pexp.mpib_hgsoc, "spearman")

#save correlation matrices 
allvars <- ls()
corvars <- allvars[grep("PEXP.COR", allvars)]
PEXP.COR.long <- CollapseCorMatrices(corvars)


#Subset to identical cell lines
PEXP.COR.long <- PEXP.COR.long[PEXP.COR.long$ID == 1, ] #subset to identical cell lines

#Change identifier
PEXP.COR.long$Dataset1[PEXP.COR.long$Dataset1 == "GDSC"] = "SANGER"
PEXP.COR.long$Dataset2[PEXP.COR.long$Dataset2 == "GDSC"] = "SANGER"

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


#Make Corrplot
MakeCorrPlot(PEXP.COR.long, order= c("BROAD", "MGHCC_BREAST", "SANGER", "NCI60", "UW_TNBC", "MPIB_HGSOC") )





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



