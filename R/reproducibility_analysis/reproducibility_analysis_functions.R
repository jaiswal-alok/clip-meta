CompareTwoStudies_Rcorr <- function(data1, data2, useCond){
  
  library(Hmisc)
  library(reshape2)
  
  # useCond = spearman
  
  #get common cells and genes 
  cells.d1.d2 <-  intersect(colnames(data1), colnames(data2))
  genes.d1.d2 <-  intersect(rownames(data1), rownames(data2))
  
  if (length(cells.d1.d2) == 0){
    cat("No overlap","\n")
  }
  if (length(cells.d1.d2) == 1){
    
    #format data matrices
    d1.mat <- data1[genes.d1.d2, cells.d1.d2]
    d2.mat <- data2[genes.d1.d2, cells.d1.d2]
    d1.d2 <- cbind(d1.mat,d2.mat)
    d1.d2 <- as.matrix(d1.d2)
    
    #compute correlation
    d1.d2.cor.list <- rcorr(d1.d2, type=useMethod)
    d1.d2.cor.list <- lapply(d1.d2.cor.list, function(X){
      d1.d2.CorMat <- X[1,2]
      return(d1.d2.CorMat)
    }) 
    d1.d2.cor.df <- do.call(cbind, d1.d2.cor.list)
    rownames(d1.d2.cor.df) <- cell
    
    #extract sample sizes for inter-dataset correlations
    rcorMat <- rbind(rcorMat, d1.d2.cor.df)
    
    cat("Mean=",mean(rcorMat[,1], na.rm=T),"\n")
    cat("Median=", median(rcorMat[,1], na.rm=T),"\n")
    cat("Dimension=", length(genes.d1.d2),"\t", length(cells.d1.d2), "\n")
    return(rcorMat)
  } 
  if (length(cells.d1.d2) >= 2){
    
    #format data matrices
    d1.mat <- data1[genes.d1.d2, cells.d1.d2]
    d1.mat <- as.data.frame(d1.mat)
    
    d2.mat <- data2[genes.d1.d2, cells.d1.d2]
    d2.mat <- as.data.frame(d2.mat)
    
    #compute correlation
    cormat.d1.d2 <- rcorr(as.matrix(d1.mat), as.matrix(d2.mat), type=useCond)
    
    #extract inter-dataset correlation values
    n1 = ncol(d1.mat)
    n2 = n1+1
    n3 = n1*2
    
    y1 <- cormat.d1.d2$r[1:n1,n2:n3]
    y1.long <- melt(y1)

    #extract sample sizes for inter-dataset correlations
    y2 <- cormat.d1.d2$n[1:n1,n2:n3]
    y2.long <- melt(y2)
    
    #extract P-values for inter-dataset correlations
    y3 <- cormat.d1.d2$P[1:n1,n2:n3]
    y3.long <- melt(y3)
    
    #Bind computed measures
    corLong <- cbind(y1.long, y2.long[,3], y3.long[,3])
    colnames(corLong) <- c("Cell1", "Cell2", "Cor", "n", "P")
    # Cell1 %in% Dataset1
    # Cell2 %in% Dataset2
    # Cor = computed correlation
    # n = number of measurements used for computing correlatin
    # P = P-value

    #correlation stats for identical cell lines
    mean.y1 <- mean(diag(y1), na.rm=T)
    median.y1 <- median(diag(y1), na.rm=T)    

    cat("# common genes  =", dim(d1.mat)[1],"\n")
    cat("# common cell lines  =", dim(d1.mat)[2],"\n")
    cat("Mean correlation =",mean.y1,"\n")
    
    return(corLong)
  }
  
}

CollapseCorMatrices <- function(cormatenv){
  
  collapseCorlist <- list()
  for (n in cormatenv){
    n.df <- get(n)
    collapseCorlist[[n]] <- n.df
  }
  
  collapseCorlist.sig <- list()
  for (n in 1:length(collapseCorlist)){
    n.df <- collapseCorlist[[n]]
    
    n.df$ID <- 0 # mark non-identical cell lines 
    n.df$ID[which(n.df$Cell1 == n.df$Cell2)] <- 1 # mark identical cell lines 
    
    #subset to statistically significant and comparisons with # gene sets >=10
    n.df <- n.df[n.df$P < 0.05, ]
    n.df <- n.df[n.df$n >= 10, ]
    
    if(nrow(n.df) == 0) next()
    listName <- names(collapseCorlist)[n]
    collapseCorlist.sig[[listName]] <- cbind(n.df,listName)
  }
  
  collapseCorlist.sig <- do.call(rbind, collapseCorlist.sig)
  collapseCorlist.sig <- as.data.frame(collapseCorlist.sig, stringsAsFactors = F)
  colnames(collapseCorlist.sig)[7] <- "Comparison"
  collapseCorlist.sig$Comparison <- as.character(collapseCorlist.sig$Comparison)
  collapseCorlist.sig$Cell1 <- as.character(collapseCorlist.sig$Cell1)
  collapseCorlist.sig$Cell2 <- as.character(collapseCorlist.sig$Cell2)
  collapseCorlist.sig$Datatype <-  sapply(strsplit(collapseCorlist.sig$Comparison, "\\."), "[[", 1)
  collapseCorlist.sig$Dataset1 <-  sapply(strsplit(collapseCorlist.sig$Comparison, "\\."), "[[", 3)
  collapseCorlist.sig$Dataset2 <-  sapply(strsplit(collapseCorlist.sig$Comparison, "\\."), "[[", 4)
  rownames(collapseCorlist.sig) <- NULL
  collapseCorlist.sig$Cor <- as.numeric(collapseCorlist.sig$Cor)
  collapseCorlist.sig$P <- as.numeric(collapseCorlist.sig$P)
  collapseCorlist.sig <- collapseCorlist.sig[!is.na(collapseCorlist.sig$Cor), ]
  return(collapseCorlist.sig)
}

MakeCorrPlot <- function(cor.long.df, order){
  library(corrplot)
  library(reshape2)
  COR.MAT <- dcast(cor.long.df[,c(9,10,3)],Dataset1~Dataset2, fill=0, fun.aggregate = mean)
  COR.MAT <- melt(COR.MAT)
  
  x1 <- COR.MAT
  colnames(x1) <- c("D1", "D2", "Cor")
  
  x2 <- COR.MAT[, c(2,1,3)]
  colnames(x2) <- c("D1", "D2", "Cor")
  
  COR.MAT <- rbind(x1,x2)
  COR.MAT <- COR.MAT[COR.MAT$Cor != 0,]
  COR.MAT <- dcast(COR.MAT,D1~D2, fill=-0.000000001, fun.aggregate = mean)
  rownames(COR.MAT) <- COR.MAT$D1
  COR.MAT <- COR.MAT[,-1]
  COR.MAT <- COR.MAT[sort(rownames(COR.MAT)), sort(colnames(COR.MAT))]
  COR.MAT <- COR.MAT[order, order]
  
  #
  M <- as.matrix(COR.MAT)
  p.mat = 1-M
  M[M == -0.000000001] = NA
  col3 <- colorRampPalette(c("blue","white","red")) 
  print(corrplot(M, type = "lower", col=col3(50),method = "color",
                 p.mat = p.mat, sig.level = 0.98, insig = "blank",
                 tl.col = "black", tl.srt = 90, tl.cex=0.8,
                 na.label="square", na.label.col="red4",outline="red4", 
                 number.cex = 1, addCoef.col = "black"))
}

PlotCV <- function(data1, data2, housekeep.genes, label1, label2){
  
  library(Hmisc)
  library(reshape2)
  library(raster)
  
  #get common cells and genes 
  cells.d1.d2 <-  intersect(colnames(data1), colnames(data2))
  genes.d1.d2 <-  intersect(rownames(data1), rownames(data2))
  
  if (length(cells.d1.d2) == 0){
    cat("No overlap","\n")
  }
  
  if (length(cells.d1.d2) >= 2){
    
    #format data matrices
    d1.mat <- data1[genes.d1.d2, cells.d1.d2]
    d1.mat <- as.data.frame(d1.mat)
    
    d2.mat <- data2[genes.d1.d2, cells.d1.d2]
    d2.mat <- as.data.frame(d2.mat)
    
    d1.CV <- apply(d1.mat,1,function(x) raster::cv(x, na.rm=T))
    d1.CV <- cbind(GeneName = names(d1.CV), CV = d1.CV)
    d1.CV <- as.data.frame(d1.CV)
    
    d2.CV <- apply(d2.mat,1,function(x) raster::cv(x[!is.nan(x)], na.rm=T))
    d2.CV <- cbind(GeneName = names(d2.CV), CV = d2.CV)
    d1.CV <- as.data.frame(d1.CV)
 
    PEXP.CV <- merge(d1.CV, d2.CV, by = "GeneName", all=T)
    PEXP.CV.df <- PEXP.CV[complete.cases(PEXP.CV), ]
    PEXP.CV.df <- as.data.frame(PEXP.CV.df, stringsAsFactors = F)
    PEXP.CV.df[,-1] <- apply(PEXP.CV.df[,-1],2,as.numeric)
    PEXP.CV.df$HK <-  0
    PEXP.CV.df$HK[PEXP.CV.df$GeneName %in% housekeep.genes$Gene] <- 1
    
    # plot
    alpha_value = 0.3
    PEXP.CV.df$Col = adjustcolor("gray48", alpha.f=alpha_value)
    alpha_value = 0.25
    PEXP.CV.df$Col[PEXP.CV.df$HK == 1]  = adjustcolor("red3", alpha.f=alpha_value)
    
    colVec <- PEXP.CV.df$Col
    xlimvar <- round(range(as.matrix(PEXP.CV.df[,2]))[2]+5)
    ylimvar <- round(range(as.matrix(PEXP.CV.df[,3]))[2]+5)
    print(plot(PEXP.CV.df[,2], PEXP.CV.df[,3], pch=19, col=colVec, xlim=c(0,xlimvar), ylim=c(0,ylimvar),
         cex=1, las=1, main="", xlab = label1, ylab=label2))
    
    corCV <- cor(PEXP.CV.df[,2], PEXP.CV.df[,3], method="spearman", use="pairwise.complete.obs")
    
    cat("# common genes  =", dim(d1.mat)[1],"\n")
    cat("# common cell lines  =", dim(d1.mat)[2],"\n")
    cat("Correlation of CV = ", corCV, "\n")
  

    
  }
  
}

