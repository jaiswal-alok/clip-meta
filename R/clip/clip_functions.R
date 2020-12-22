QuantNormScale <- function(dataInput){
  # dataInput should be a matrix or a data frame, rows = Genes. columns = Cell lines
  #Quantile normalization for making two distributions identical in statistical properties for all cell lines

  library(preprocessCore)

  print(dim(dataInput))
  dataInput <- apply(dataInput,2, function(x){
    nan_vec <- is.nan(x)
    x[nan_vec] <- NA
    return(x)
  } )
  dataInput.Quant <- normalize.quantiles(as.matrix(dataInput))

  rownames(dataInput.Quant) <- rownames(dataInput)
  colnames(dataInput.Quant) <- colnames(dataInput)
  
  #mean scaling of each gene on quantile normalized data 
  dataInput.Quant.scaled <- apply(t(dataInput.Quant),2,function(e) scale(e, center = T))
  dataInput.Quant.scaled <- t(dataInput.Quant.scaled)
  
  rownames(dataInput.Quant.scaled) <- rownames(dataInput)
  colnames(dataInput.Quant.scaled) <- colnames(dataInput)
  return(dataInput.Quant.scaled)
}

ResidueAvg <- function(input){
  #input should be a matrix/dataframe, rows = genes. columns = cell lines
  
  #Average for genes with multiple phosphosites
  inputDF <- as.data.frame(cbind(GeneName=rownames(input), input), stringsAsFactors = F)
  inputDF$GeneName <- as.character(inputDF$GeneName)
  inputDF$GeneName <- sapply(strsplit(inputDF$GeneName, "\\_"), "[[", 1)
  inputDF.list <- split(inputDF, inputDF$GeneName)
  inputDF.list.mean <- do.call(rbind, lapply(inputDF.list, function(x){ y <- colMeans(x[,-1], na.rm = T)}))
  outputDF <- as.data.frame(inputDF.list.mean, stringsAsFactors = F)
  return(outputDF)
  
}

BinarizeCNV <- function(CNVData, threshDEL, threshAMP){
  #CNVData  = copy number data matrix/dataframe
  #threshDEL = threshold for calling deletion
  #threshAMP =  threshold for calling amplification
  
  
  # DELETION 
  CNVData_DEL <- apply(CNVData, 2, function(x){
    y <- x
    #for binary thresholding
    y[x <= -threshDEL] <- 1
    y[x > -threshDEL] <- 0
    return(y)
  })
  rownames(CNVData_DEL) <- rownames(CNVData)
  
  # AMPLIFICATION 
  CNVData_AMP <- apply(CNVData, 2, function(x){
    y <- x
    #For binary thresholding
    y[x >= threshAMP] <- 1
    y[x < threshAMP] <- 0
    return(y)
  })
  rownames(CNVData_AMP) <- rownames(CNVData)
  
  return(list(DEL=CNVData_DEL, AMP=CNVData_AMP))
}

CalculatePropScores <- function(binCNVdata){
  #binCNVdata = binarized copy number data, output of BinarizeCNV
  
  #Proportion of CNV alteration as the measure for context specific property
  geneSums <- rowSums(binCNVdata, na.rm = T)
  geneSums.Freq <- table(geneSums)
  
  geneProps.Mat <- c()
  for (s in 1:length(geneSums.Freq)){
    
    sName <- as.numeric(names(geneSums.Freq)[s])
    sIdx <- which(geneSums == sName)
    #print(c(s, sName))
    prop.value <- sName/ncol(binCNVdata)
    
    df.sIdx <- binCNVdata[sIdx, , drop=F]
    if(nrow(df.sIdx) == 1) df.sIdx[1,][df.sIdx[1,] == 1] <- prop.value
    if(nrow(df.sIdx) > 1) {
      df.sIdx <- apply(df.sIdx,2,function(E){
        E1 = E
        E1[E == 1] <- prop.value
        E1[E == 0] <- 0
        return(E1)
      })
    }
    
    geneProps.Mat <- rbind(geneProps.Mat, df.sIdx)
    
  }
  
  geneProps.Mat <- geneProps.Mat[sort(rownames(geneProps.Mat)), sort(colnames(geneProps.Mat))]
  
  return(geneProps.Mat)
}

IntegrateCCS <- function(listDataStudy, Modality, cellNames, SingleSelectThreshold, SelectIndx, SelectThreshold){
  # listDataStudy = Quantile normalized and scaled datasets representing CCS scores from various studies
  # Modality = Modality name
  # cellNames =  Cell lines on which meta-analysis is performed
  # SingleSelectThreshold  = Threshold for selection of genes from a cellline profiled in only 1 study
  # SelectThreshold =  percentage of false prediction (pfp) threshold for selection of CCS genes from Rank product analysis
  
  cellIntList <- list()
  cellIntResults <- list()
  
  # for each cell line
  for (c in cellNames){
    print(c)
    cellStudyList <- list()
    
    # for each cell line extract data from each study
    for (i in 1:length(listDataStudy)){
      dataMat <- listDataStudy[[i]]
      dataName <- names(listDataStudy)[i]
      cellName <- paste("^", c, "$", sep="")
      cellName.Idx <- grep(cellName, colnames(dataMat))
      if(length(cellName.Idx) == 0) next()
      if(length(cellName.Idx) > 0)
      {
        dataMat.Cell <- as.data.frame(cbind(GeneID = rownames(dataMat),dataMat[,cellName.Idx]), stringsAsFactors = F)
        colnames(dataMat.Cell)[2] <- dataName
        dataMat.Cell[[2]] <- as.numeric(dataMat.Cell[[2]])
        cellStudyList[[dataName]] <- dataMat.Cell
      }
    }
    
    # merge CCS/PS scores for a cell line from various sites
    cellIntMat <- Reduce(function(...) merge(...,by = "GeneID", all=T), cellStudyList)
    cellIntMat <- unique(cellIntMat)
    
    if(is.null(cellIntMat)) next()
    cellIntList[[c]] <- cellIntMat
    
    # Meta-analysis of CCS/PS scores
    if(Modality %in% c("CNV_DEL", "CNV_AMP", "MUT")) {
      cellIntResults[[c]] <- PerformIntBin(cellIntMat, Modality)
    } else {
      cellIntResults[[c]] <- PerformIntCont(cellIntMat, 
                                           SelectIndx, 
                                           SelectThreshold,
                                           SingleSelectThreshold, 
                                           Modality)
    }
  }
  
  return(list(cellIntList= cellIntList, cellIntResults = cellIntResults))
}

PerformIntBin <- function(cellIntMat, Modality){
  #cellIntMat = dataframe of geneXstudy for each cell line
  #Modality = any binary data modalities i.e. MUT, CNv_AMP, CNV_DEL
  
  # Integration of PS scores for binary modalities
  
  data.sub = cellIntMat
  rownames(data.sub) <- data.sub$GeneID
  
  #For cell lines with data from 1 study
  if(ncol(data.sub) == 2){
    data = data.sub
    Y = X = data[,2]
    X[X==0] <- NA
    Y[Y==0] <- NA
    Y[X>0.1] <- 0
    Y[X<=0.1] <- 1 # Proportion scores <= 0.1 defined as CCS 
    data[,2] = Y
    colnames(data) = c("GeneID", "Sum")
    data <- data[which(data$Sum == 1),]
    colnames(data) <- paste(Modality, colnames(data), sep="_")
  }
  
  #For cell lines with data from >= 2 studies
  if(ncol(data.sub) > 2){
    
    #Remove Genes with no alterations
    data = data.sub[rowSums(data.sub[,-1], na.rm = T) != 0, ]
    rownames(data) <- data$GeneID
    
    #Mark genes with proportion scores <= 10% in a dataset as CCS genes
    if(dim(data)[1] == 0) data = data
    if(dim(data)[1] > 0) data[,-1] <- apply(data[,-1],2, function(X){
      Y <- X
      X[X==0] <- NA
      Y[Y==0] <- NA
      Y[X>0.1] <- 0
      Y[X<=0.1] <- 1
      return(Y)
    })
    data$Sum <- rowSums(data[,-1], na.rm = T)
    
    #Select genes with CCS evidence in >= 2 datasets
    data <- data[data$Sum >= 1,]
    data <- data[,c("GeneID", "Sum")]
    colnames(data) <- paste(Modality, colnames(data), sep="_")
  }
  
  return(list(Table1 = data, Table1.genes= data[,1]))
}

PerformIntCont <- function(cellIntMat, SelectIndx, SelectThreshold, SingleSelectThreshold, Modality){
  
  #cellIntMat = dataframe of geneXstudy for each cell line
  #Modality = any continuous data modalities i.e. GEXP. METH. FUNC, TAS and PHOS
  #SelectThreshold = pfp threshold for selection of CCS genes after Rank integration
  #SingleSelectThreshold = threshold for selection of top N CCS genes for cell lines profiled in only 1 study
  
  
  # Non-parametric integration of CCS scores from continuous modalities
  
  library(RankProd)
  
  rownames(cellIntMat) <- cellIntMat$GeneID
  
  #For cell lines with data from 1 study
  if(ncol(cellIntMat) == 2){
    
    # Table 1: list genes that are CCS_DOWN (except for FUNC modality, where it's CCS_UP)
    Table1 <- cellIntMat
    Table1 <- Table1[complete.cases(Table1), ]
    Table1 <- Table1[order(Table1[[2]], decreasing = F), ]
    Table1.genes <- Table1$GeneID[1:SingleSelectThreshold]
    rownames(Table1) <- NULL
    colnames(Table1) = c("GeneID", "Score")
    colnames(Table1) <- paste(Modality, colnames(Table1), sep="_")
    
    # Table 2: list genes that are CCS_UP (except for FUNC modality, where it's CCS_DOWN)
    Table2 <- cellIntMat
    Table2 <- Table2[complete.cases(Table2), ]
    Table2 <- Table2[order(Table2[[2]], decreasing = T), ]
    Table2.genes <- Table2$GeneID[1:SingleSelectThreshold]  #Select top genes above the specified threshold
    rownames(Table2) <- NULL
    colnames(Table2) = c("GeneID", "Score")
    colnames(Table2) <- paste(Modality, colnames(Table2), sep="_")
    
  }
  
  #For cell lines with data from >= 2 studies
  if(ncol(cellIntMat) > 2){
    cellIntMatF <- cellIntMat
    cellIntMatF[,-1] <- apply(cellIntMatF[,-1],2, function(x) rank(x, na.last = "keep"))
    cellIntMatF <- as.data.frame(cellIntMatF[,-1,drop =F])
    
    dType.data.cl <- rep(1, ncol(cellIntMatF)) #single class rank product integration
    cellIntMat.RP <- RankProducts(cellIntMatF, dType.data.cl, calculateProduct = T, 
                            logged=FALSE,na.rm=FALSE, plot=FALSE, rand=123) #Rank product integration of ranked CCS scores
    cellIntMat.RP.Top <- topGene(cellIntMat.RP, cutoff=100, method="pfp", logged=FALSE, gene.names=rownames(cellIntMatF))
    cellIntMat.RP.Top <- lapply(cellIntMat.RP.Top, function(x) as.data.frame(x))
    
    # Table 1: list genes that are CCS_DOWN (except for FUNC modality, where it's CCS_UP)
    Table1 <- cellIntMat.RP.Top$Table1
    Table1 <- Table1[complete.cases(Table1), ]
    Table1.genes <- rownames(Table1)[Table1[[SelectIndx]] < SelectThreshold]
    colnames(Table1) <- paste(Modality, colnames(Table1), sep="_")
    Table1 <- cbind(GeneID = rownames(Table1), Table1)
    rownames(Table1) <- NULL
    Table1=Table1[,c(1,3,5,6)]
    
    # Table 2: list genes that are CCS_UP (except for FUNC modality, where it's CCS_DOWN)
    Table2 <- cellIntMat.RP.Top$Table2
    Table2 <- Table2[complete.cases(Table2), ]
    Table2.genes <- rownames(Table2)[Table2[[SelectIndx]] < SelectThreshold]
    colnames(Table2) <- paste(Modality, colnames(Table2), sep="_")
    Table2 <- cbind(GeneID = rownames(Table2), Table2)
    rownames(Table2) <- NULL
    Table2 = Table2[,c(1,3,5,6)]
    
    #Remove Genes with lower PFP value if Table1 and Table2 overlap
    commonT1T2genes <- intersect(Table1.genes, Table2.genes)
    for (gene in commonT1T2genes){
      #Fetch PFP values
      pfp1 <- Table1[[3]][Table1$GeneID == gene]
      pfp2 <- Table2[[3]][Table2$GeneID == gene]
      if (pfp1 < pfp2) Table2.genes = setdiff(Table2.genes, gene)
      if (pfp2 < pfp1) Table1.genes = setdiff(Table1.genes, gene)
    }
    
  }
  
  return(list(Table1=Table1, Table2 = Table2, Table1.genes= Table1.genes, Table2.genes=Table2.genes))
}

`%ni%` <- Negate(`%in%`) 

RearrangeCCS <- function(x){
  
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
}

ProcessCCSlist <- function(CCSMat.list, keepMode){
  
  #CCSMat.list = list of CCS genes after integration, output of IntegrateCCS
  #keepMode = Obligatorily include a specific modality to identify rCCS genes
  
  # select rCCS genes
  rCCSmat <- lapply(CCSMat.list, function(XDF){
    if (keepMode %in% "essential genes"){ # for Genetic interaction analysis, obligate criteria for rCCS genes
      XDF <- XDF[which(XDF$FUNC_UP == 1 | XDF$TAS_UP == 1) , ]
    } else { XDF }
    
    XDF.Genes <- XDF$GeneID[XDF$Sum >=2] # Select genes with CCS evidence >= 2 modalities
    if(length(XDF.Genes)==0) XDF.Genes <- NA
    XDF.Genes <- cbind("CellName", XDF.Genes)
    colnames(XDF.Genes) <- c("CellName", "GeneID")
    return(XDF.Genes)
  })
  for (i in 1:length(rCCSmat)){
    c = names(rCCSmat)[i]
    rCCSmat[[i]][,1] <- c 
    
  }

  # convert list to matric
  rCCSmat <- do.call(rbind, rCCSmat)
  rCCSmat <- as.data.frame(rCCSmat, stringsAsFactors = F)
  rCCSmat$Var <- 1
  rCCSmat <- unique(rCCSmat)
  rCCSmat <- rCCSmat[complete.cases(rCCSmat), ]
  rCCSmat.Wide <- dcast(rCCSmat, GeneID~CellName, fill=0)
  rCCSmat.Wide[,-1] <- apply(rCCSmat.Wide[,-1],2,as.numeric)
  rCCSmat.Wide$Sum <- rowSums(rCCSmat.Wide[,-1])
  return(rCCSmat.Wide)
}

PerformFisherExact <- function(data, group1.cells, group2.cells, side){
  
  #data = rCCS matrix/dataframe, rows = genes and cols = cell lines
  #side = alternative hypothesis for Fisher's Exact test
  
  group1.data <- data[ ,na.omit(match(group1.cells, colnames(data)))]
  group2.data <- data[ ,na.omit(match(group2.cells, colnames(data)))]
  print(ncol(group1.data))
  print(ncol(group2.data))
  
  all.data <- cbind(group1.data, group2.data)
  count <- rowSums(all.data)
  
  #keep genes present in atleast 1 cell line
  all.data <-  all.data[which(count > 0), ]
  
  #Fisher's exact test between two groups
  group <- rep(NA, ncol(all.data))
  n1 <- ncol(group1.data)
  n2 <- n1+1
  n3 <- ncol(all.data)
  group[1:n1] <- 1
  group[n2:n3] <- 0
  group <- as.factor(group)
  
  fisher.exact.pval <- apply(all.data,1, function(x){
    x <- as.factor(as.numeric(x))
    x.fisher <- fisher.test(x,group, alternative = side)
    x.fisher$p.value
  })
  return(sort(fisher.exact.pval))
}

PerformLIMMA  <- function(data, group1.cells, group2.cells, rtrvIDx){
  #data = rCCS matrix/dataframe, rows = genes and cols = cell lines
  #rtrvIDx = index to choose from resulting DE table 
  
  library(limma)
  
  group1.data <- data[ ,na.omit(match(group1.cells, colnames(data)))]
  group2.data <- data[ ,na.omit(match(group2.cells, colnames(data)))]
  
  nc1 <- ncol(group1.data)
  nc2 <- ncol(group2.data)
  evalData <- cbind(group1.data, group2.data)

  dsgn <- c(rep(0, nc1), rep(1, nc2))
  
  #Remove Genes with NA in one group
  fit <- lmFit(evalData, design = model.matrix(~dsgn))
  genes.onegroup.na <- is.na(fit$coefficients[,2])
  evalData <- evalData[!genes.onegroup.na, ]
  
  #Refit linear model
  fit <- lmFit(evalData, design = model.matrix(~dsgn))
  fit <- eBayes(fit, trend=TRUE)
  fitTable <- topTable(fit, coef=2, number = Inf)
  fitTable <- cbind(Genes=rownames(fitTable), fitTable)
  
  colnames(fitTable)[1] <- "GeneID"
  
  fitTable.Sub <- fitTable[, c(1,rtrvIDx)]

  return(fitTable.Sub)
}

PerformRP <- function(data){
  
  #data = geneXStudy data frame of OES scores
  #input for rank product integration
  library(RankProd)
  dType.data.cl <- rep(1, ncol(data)) # define single class
  data.RP <- RankProducts(data, dType.data.cl, calculateProduct = T, 
                          logged=FALSE,na.rm=FALSE, plot=FALSE, rand=123) #Rank product
  data.RP.Top <- topGene(data.RP, cutoff=100, method="pfp", logged=FALSE, gene.names=rownames(data))
  data.RP.Top <- lapply(data.RP.Top, function(x) as.data.frame(x))
  
  return(data.RP.Top)
}

PerformDE <- function(brca.subtype, modal.list, modal_name,  ref.drivers ){
  #brca.subtype = data frame of breast cancer subtype info for cell lines in each column
  #modal.list = list of data sets from several studies for  each modality 
  #modal_name = data modality name
  #ref.drivers = known gold standard driver genes
  
  drivers.list <- list()
  
  nc2 <- length(colnames(brca.subtype))
  for (c in c(2:nc2)){
    
    #Select cells belonging to subtypes
    subtype.compare.name <- colnames(brca.subtype)[c]
    subtype.pos.cells <- brca.subtype$CellName[which(brca.subtype[,c]==1)]
    subtype.neg.cells <- brca.subtype$CellName[brca.subtype[,c]==0]
  
    
    # LIMMA analysis (get fold change)
    dE.LIMMA.FC <- lapply(modal.list, function(XDF){
      Pos.cells <- intersect(subtype.pos.cells, colnames(XDF))
      Neg.cells <- intersect(subtype.neg.cells, colnames(XDF))
      dFC <- PerformLIMMA(XDF, Pos.cells, Neg.cells, 2)
      return(dFC)
    })
    dE.LIMMA.FC <- Reduce(function(...) merge(..., by = "GeneID", all=T), dE.LIMMA.FC)
    colnames(dE.LIMMA.FC)[-1] <- names(modal.list)
    dE.LIMMA.FC[,-1] <- apply(dE.LIMMA.FC[,-1, drop=F],2,as.numeric)
    rownames(dE.LIMMA.FC) <- dE.LIMMA.FC$GeneID
    #Estimate average fold change across datasets for direction of change
    library(psych)
    fc.mean <- apply(dE.LIMMA.FC[,-1],1, function(x) mean(x, na.rm=T))
    
    # LIMMA analysis (get B stats)
    dE.LIMMA.B <- lapply(modal.list, function(XDF){
      Pos.cells <- intersect(subtype.pos.cells, colnames(XDF))
      Neg.cells <- intersect(subtype.neg.cells, colnames(XDF))
      dB <- PerformLIMMA(XDF, Pos.cells, Neg.cells, 7)
      return(dB)
    })
    dE.LIMMA.B <- Reduce(function(...) merge(..., by = "GeneID", all=T), dE.LIMMA.B)
    colnames(dE.LIMMA.B)[-1] <- names(modal.list)
    dE.LIMMA.B[,-1] <- apply(dE.LIMMA.B[,-1, drop=F],2,as.numeric)
    rownames(dE.LIMMA.B) <- dE.LIMMA.B$GeneID
    #RankProduct of B stats
    dE.Limma.B.RP <- dE.LIMMA.B
    dE.Limma.B.RP[,-1] <- apply(dE.Limma.B.RP[,-1],2, function(x) rank(x, na.last = "keep"))
    dE.Limma.B.RP <- PerformRP(dE.Limma.B.RP[,-1])
    
    # Extract statistically significant genes up-regulated/amplified/hypermethylated in Subtype of interest 
    fdr.thresh = 0.05
    # "UP" driver genes
    dE.genes.UP <- rownames(dE.Limma.B.RP$Table2)[dE.Limma.B.RP$Table2$pfp <= fdr.thresh]
    dE.genes.UP <- intersect(dE.genes.UP, ref.drivers)
    dE.genes.UP <- fc.mean[dE.genes.UP] 
    if (modal_name %in% c("GEXP", "PEXP", "CNV", "METH")) dE.genes.UP <- dE.genes.UP[dE.genes.UP < 0]
    if (modal_name %in% c("FUNC")) dE.genes.UP <- dE.genes.UP[dE.genes.UP > 0]
    
    fc_thresh = 0 # threshold for average fold change in expression
    dE.genes.UP <- dE.genes.UP[abs(dE.genes.UP) > fc_thresh ]
    dE.genes.UP <- names(dE.genes.UP)
    
    modality <- paste("DRIVERS", subtype.compare.name, modal_name, sep="_")
    drivers.list[[modality]] <- dE.genes.UP
  }
  return(drivers.list)
}

ExtractTopFeatures <- function(mofaweightdf){
  unique.views <- unique(mofaweightdf$view)
  factor.vec <- unique(mofaweightdf$factor)
  
  mofa.top.weights <- c()
  for (v in unique.views){
    mofa.weights.sub <- mofaweightdf[mofaweightdf$view == v, ]
    for (f in factor.vec){
      mofa.weights.sub2 <- mofa.weights.sub[mofa.weights.sub$factor == f, ]
      mofa.weights.sub2 <- mofa.weights.sub2[order(mofa.weights.sub2$value, decreasing = T), ]
      mofa.weights.top10.pos <- mofa.weights.sub2[1:10, ]
      mofa.weights.top10.pos$Direction = "Positive"
      mofa.weights.sub2 <- mofa.weights.sub2[order(mofa.weights.sub2$value, decreasing = F), ]
      mofa.weights.top10.neg <- mofa.weights.sub2[1:10, ]
      mofa.weights.top10.neg$Direction = "Negative"
      mofa.weights.top20 <- rbind(mofa.weights.top10.pos, mofa.weights.top10.neg)
      mofa.top.weights <- rbind(mofa.top.weights, mofa.weights.top20)
    }
  }
  mofa.top.weights <- as.data.frame(mofa.top.weights, stringsAsFactors =F)
  mofa.top.weights$FeatureID <- sapply(strsplit(as.character(mofa.top.weights$feature), "\\_"), "[[", 1)
  return(mofa.top.weights)
}


