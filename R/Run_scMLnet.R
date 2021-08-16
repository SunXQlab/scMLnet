#Part1:getHighExpGene
#output:GeneList
getHighExpGene <- function(GCMat,barCluTable,CluNum.1,CluNum.2,pval,logfc,cores)
{
  
  #Function
  getNormData <- function(data.use){
    
    LogNorm <- function(data, scale_factor, display_progress = TRUE) {
      .Call('_Seurat_LogNorm', PACKAGE = 'Seurat', data, scale_factor, display_progress)
    }
    
    LogNormalize <- function(data, scale.factor = 1e4, display.progress = TRUE)
    {
      if (class(x = data) == "data.frame")
      {
        data <- as.matrix(x = data)
      }
      if (class(x = data) != "dgCMatrix") {
        data <- as(object = data, Class = "dgCMatrix")
      }
      # call Rcpp function to normalize
      if (display.progress) {
        cat("Performing log-normalization\n", file = stderr())
      }
      norm.data <- LogNorm(data, scale_factor = scale.factor, display_progress = display.progress)
      colnames(x = norm.data) <- colnames(x = data)
      rownames(x = norm.data) <- rownames(x = data)
      return(norm.data)
    }
    
    data.test <- LogNormalize(data.use)
    
    return(data.test)
    
  }
  
  ChooseGene <- function(data.use,cells.1,cells.2,logfc){3
    
    #Function
    calpct <- function(data.use,cells,thresh.min=0){
      
      if(length(genes.use) > 20000)
      {
        num=ceiling(length(genes.use)/5000)
        
        data.temp <- c()
        for(i in seq(1,num,1))
        {
          start=1+5000*(i-1)
          if(i == num)
          {
            end=length(genes.use)
          }
          else
          {
            end=5000*i
          }
          
          tempGeneUse=genes.use[start:end]
          
          data.tempX <- round(
            x = apply(
              X = data.use[tempGeneUse, cells, drop = F],
              MARGIN = 1,
              FUN = function(x) {
                return(sum(x > thresh.min) / length(x = x))
              }
            ),
            digits = 3
          )#round
          data.temp <- c(data.temp,data.tempX)
        }#for
      }else{
        data.temp <- round(
          x = apply(
            X = data.use[, cells, drop = F],
            MARGIN = 1,
            FUN = function(x) {
              return(sum(x > thresh.min) / length(x = x))
            }
          ),
          digits = 3
        )#round
      }
      
      return(data.temp)
    }
    
    genes.use <- rownames(data.use)
    print(paste0("gene:",length(genes.use)))
    
    #choose Gene based on percent expressed
    min.pct <- 0.05
    min.diff.pct <- -Inf
    
    ##calculate pct
    data.temp1 <- calpct(data.use,cells.1)
    data.temp2 <- calpct(data.use,cells.2)
    data.alpha <- cbind(data.temp1, data.temp2)
    colnames(x = data.alpha) <- c("pct.1","pct.2")
    
    ##max pct
    alpha.max <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.max) <- rownames(x = data.alpha)
    genes.use <- names(x = which(x = alpha.max > min.pct))
    
    ##diff pct
    alpha.diff <- alpha.max - apply(X = data.alpha, MARGIN = 1, FUN = min)
    genes.use <- names(
      x = which(x = alpha.max > min.pct & alpha.diff > min.diff.pct)
    )
    
    #choose Gene based on average difference
    logfc.threshold <- logfc
    print("logfc.threshold:")
    print(logfc.threshold)
    pseudocount.use <- 1
    
    ##diff expr
    data.1 <- apply(X = data.use[genes.use, cells.1, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
    data.2 <- apply(X = data.use[genes.use, cells.2, drop = F], MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
    total.diff <- (data.1 - data.2)
    genes.diff <- names(x = which(x = total.diff > logfc.threshold))
    
    #filter
    genes.use <- intersect(x = genes.use, y = genes.diff)
    
    #output final gene Table(gene,logfc)
    geneTable <- data.frame(gene = genes.use,logfc = total.diff[genes.use])
    
    return(geneTable)
  }
  
  do.ttest <- function(data.test,genes.use,cells.1,cells.2,cores){
    
    #Function
    getpval <- function(gene){
      x <- as.vector(data.test[gene, cells.1])
      y <- as.vector(data.test[gene, cells.2])
      result <- t.test(x, y)
      return(result$p.value)
    }
    
    if(cores == 1){
      
      print("do T-test")
      p_val <- lapply(genes.use,getpval)
      p_val <- unlist(p_val)
      
    }else{
      
      print("T-test in parallel")
      cl <- makeCluster(cores)
      clusterEvalQ(cl, library(stats))
      clusterEvalQ(cl, library(Matrix))
      clusterExport(cl, varlist=c("data.test", "cells.1", "cells.2"), envir=environment())
      p_val <- parLapply(cl, genes.use, getpval)
      p_val <- unlist(p_val)
      stopCluster(cl)
      
    }
    
    gc()
    return(p_val)
  }
  
  #Main
  stri <- paste0("get High Exp Gene in ",CluNum.1)
  print(stri)
  barListResult <- getBarList(CluNum.1,CluNum.2,barCluTable)
  cells.1 <- barListResult[[1]]
  cells.2 <- barListResult[[2]]
  stri <- paste0(CluNum.1,":",length(cells.1))
  print(stri)
  
  data.use <- GCMat
  ChoGeneTable <- ChooseGene(data.use,cells.1,cells.2,logfc)
  genes.use <- rownames(ChoGeneTable)
  
  #normalize
  data.test <- getNormData(GCMat)
  data.test <- data.test[genes.use,c(cells.1,cells.2)]
  
  #do t test
  p_val <- do.ttest(data.test,genes.use,cells.1,cells.2,cores)
  
  #get final Gene
  GenePval <- data.frame(Gene = genes.use,LogFC = ChoGeneTable$logfc ,Pval = p_val,row.names = genes.use)
  GenePval1 <- GenePval[which(GenePval$Pval < pval),]
  
  HighExpGene <- GenePval1$Gene
  HighExpGene <- as.vector(HighExpGene)
  HighExpGeneNum <- length(HighExpGene)
  
  gc()
  stri <- paste("find high gene num:",HighExpGeneNum,sep="")
  print(stri)
  print("-----------------------")
  return(HighExpGene)
}
#Part1.end

#Part2.getLigRec
getLigRec <- function(LigRecTable,LigCluGene,RecCluGene)
{
  #read Ligand Receptor
  LigRecTable <- read.table(LigRecTable,header = TRUE,sep="\t")
  LigGene <- as.vector(unique(LigRecTable$Ligand))
  RecGene <- as.vector(unique(LigRecTable$Receptor))
  TotLigRec <- as.vector(unique(LigRecTable$Key))
  
  #get cluster final Lig_Rec result
  LRList <- c()
  LigHighGene <- intersect(LigGene,LigCluGene)
  RecHighGene <- intersect(RecGene,RecCluGene)
  
  for(x in LigHighGene)
  {
    for(y in RecHighGene)
    {
      stri=paste(x,"_",y,sep="")
      if(length(intersect(stri,TotLigRec)) == 1)
      {
        LRList <- c(LRList,stri)
      }
    }
  }
  
  return(LRList)
  
}

#Part2.end

#Part3.TFTargetFisher
TFTargetFisher <- function(TFTableFile,Tup,TotalGene)
{
  TFTable <- read.table(TFTableFile,sep="\t",header = TRUE)
  TFGene <- unique(as.vector(TFTable$TF))
  
  tempTable <- NULL
  TFTarGeneTable <- NULL
  for(v in seq(1,length(TFGene),1))
  {
    ThisTF <- TFGene[v]
    TarGenes <- as.vector(TFTable[which(TFTable$TF == ThisTF),]$Target)  #This TF all target genes
    
    TarGeneTup <- intersect(TarGenes,Tup)  #This TF's Up target Gene
    
    #get a,b,c,d
    a <- TarGeneTup
    b <- setdiff(Tup,a)
    c <- setdiff(TarGenes,a)
    
    abc <- union(a,b)
    abc <- union(abc,c)
    
    d <- setdiff(TotalGene,abc)
    
    ab <- union(a,b)
    cd <- union(c,d)
    ac <- union(a,c)
    
    #get p
    p <- (choose(length(ab),length(a)) * choose(length(cd),length(c))) / choose(length(TotalGene),length(ac))
    temp <- c(ThisTF,p)
    tempTable <- rbind(tempTable,temp)
    
    if(length(TarGeneTup) > 0)
    {
      for(t in seq(1,length(TarGeneTup),1))
      {
        key <- paste(ThisTF,TarGeneTup[t],sep="_")
        temp2 <- c(ThisTF,TarGeneTup[t],key)
        TFTarGeneTable <- rbind(TFTarGeneTable,temp2)
      }
    }
    
  }
  
  tempTable <- as.data.frame(tempTable)
  colnames(tempTable) <- c("TF","pval")
  FinalTable <- tempTable[which(as.vector(tempTable$pval) < 0.05),]
  colnames(FinalTable) <- c("TF","pval")
  
  #TF-TargetGene,all TargetGene in Up Gene List
  colnames(TFTarGeneTable) <- c("TF","TargetGene","key")
  TFTarGeneTable <- as.data.frame(TFTarGeneTable)
  
  #result:TFTarGeneTable(a TF have TargetGene up),TF p value < 0.05
  result <- list(TFTarGeneTable,FinalTable)
  
  return(result)
}
#Part3.end

#Part4.Receptor-TF Fisher Test
RecTFFisher <- function(RecTFTableFile,TFP)
{
  
  #main
  RecTable <- read.table(RecTFTableFile,sep="\t",header = TRUE)
  RecAll <- as.vector(unique(RecTable$Receptor))
  TFAll <- as.vector(unique(RecTable$TF))
  
  tempTable <- NULL
  for(v in seq(1,length(RecAll),1))
  {
    #for every receptor
    ThisRec <- RecAll[v]
    ThisRecTF <- unique(as.vector(RecTable[which(RecTable$Receptor == ThisRec),]$TF))  #This Receptor's all TF
    
    #get a b c d
    a <- intersect(TFP,ThisRecTF)
    b <- setdiff(TFP,a)
    c <- setdiff(ThisRecTF,a)
    
    ab <- union(a,b)
    ac <- union(a,c)
    
    abc <- union(a,b)
    abc <- union(abc,c)
    d <- setdiff(TFAll,abc)
    
    cd <- union(c,d)
    
    #get p value
    p <- (choose(length(ab),length(a)) * choose(length(cd),length(c))) / choose(length(TFAll),length(ac))
    temp <- c(ThisRec,p)
    tempTable <- rbind(tempTable,temp)
  }
  
  #get final table
  tempTable <- as.data.frame(tempTable)
  colnames(tempTable) <- c("Receptor","pval")
  FinalTable <- tempTable[which(as.vector(tempTable$pval) < 0.05),]
  colnames(FinalTable) <- c("Receptor","pval")
  
  #get XZRec of TF list
  XZRecTFList <- c()
  XZRec <- as.vector(FinalTable$Receptor)
  for(v in seq(1,length(XZRec),1))
  {
    #for every XZ receptor
    ThisRec <- XZRec[v]
    ThisRecTF <- unique(as.vector(RecTable[which(RecTable$Receptor == ThisRec),]$TF))  #This Receptor's all TF
    XZTF <- intersect(ThisRecTF,TFP)
    
    if(length(XZTF) > 0)
    {
      for(i in seq(1,length(XZTF),1))
      {
        thisXZTF <- paste(ThisRec,XZTF[i],sep="_")
        XZRecTFList <- c(XZRecTFList,thisXZTF)
      }
    }
    
  }
  
  result <- list(FinalTable,XZRecTFList)
  return(result)
  
}
#Part4.end

#Part5.get correlation
getCorR <- function(GCMat,keys,cells)
{
  
  GCMat <- GCMat[,cells]
  
  #check key
  keys <- unique(keys)
  key_df <- data.frame(key=keys)
  key_df$gene.1 <- unlist(lapply(keys,function(key){strsplit(key,"_")[[1]][1]}))
  key_df$gene.2 <- unlist(lapply(keys,function(key){strsplit(key,"_")[[1]][2]}))
  
  key_df <- key_df[ key_df$gene.1 %in% rownames(GCMat) & key_df$gene.2 %in% rownames(GCMat),]
  key_df <- key_df[ Matrix::rowMeans(GCMat[key_df$gene.1,])>0 & Matrix::rowMeans(GCMat[key_df$gene.2,])>0 ,]
  
  #cal cor
  cor_df <- matrix(nrow = nrow(key_df), ncol = 2)
  for(i in 1:nrow(key_df)){
    
    key_info <- key_df[i,]
    
    cor <- cor.test(GCMat[key_info$gene.1,],GCMat[key_info$gene.2,],method = c("kendall"))
    pval <- signif(cor$p.value,digits = 3)
    R <- signif(cor$estimate,digits = 2)
    
    cor_df[i,] <- c(R,pval)
    
  }
  
  rownames(cor_df) <- key_df$key
  colnames(cor_df) <- c("R","pval")
  cor_df <- cor_df[ cor_df[,1] > 0 &  cor_df[,2] < 0.05,]
  
  return(rownames(cor_df))
  
}
#Part5 end

#function1:get Clu barcode
getBarList <- function(Aclu,Bclu,barcluTable)
{
  AllCluster <- unique(barcluTable$Cluster)
  
  AcluBar <- as.vector(barcluTable$Barcode[which(barcluTable$Cluster == Aclu)])
  BcluBar <- as.vector(barcluTable$Barcode[which(barcluTable$Cluster == Bclu)])
  
  result <- list(AcluBar,BcluBar)
  return(result)
}
#function1:end

#function2:get Pair List
getPairGeneList <- function(PairList,AB)
{
  GeneList <- c()
  for(v in seq(1,length(PairList),1))
  {
    stri <- PairList[v]
    if(AB == "A")
    {
      thisGene <- unlist(strsplit(stri,split="_"))[1]
    }else if(AB == "B")
    {
      thisGene <- unlist(strsplit(stri,split="_"))[2]
    }else{
      print("type error!")
      thisGene <- ""
    }
    GeneList <- c(GeneList,thisGene)
  }
  GeneList <- unique(GeneList)
  return(GeneList)
}
#function2:end

#function3:search in pair
SearchInPair <- function(PairList,SearchList,SearAB)
{
  FPairList <- c()
  if(SearAB == "A")
  {
    for(v in seq(1,length(PairList),1))
    {
      if(unlist(strsplit(PairList[v],"_"))[1] %in% SearchList)
      {
        FPairList <- c(FPairList,PairList[v])
      }
    }
  }else if(SearAB == "B")
  {
    for(v in seq(1,length(PairList),1))
    {
      if(unlist(strsplit(PairList[v],"_"))[2] %in% SearchList)
      {
        FPairList <- c(FPairList,PairList[v])
      }
    }
  }
  return(FPairList)
}
#function3:end

#function4:pair to table
PairToTable <- function(PairList,LName,RName)
{
  Table <- NULL
  for(v in seq(1,length(PairList),1))
  {
    temp <- unlist(strsplit(PairList[v],"_"))
    tempStri <- c(temp[1],temp[2],PairList[v])
    Table <- rbind(Table,tempStri)
  }
  Table <- as.data.frame(Table)
  colnames(Table) <- c(LName,RName,"key")
  return(Table)
  
}
#function4:end

###################

#' @title Generate Multi-layer Signal Networks
#'
#' @description
#' This function constructs the Ligand_Receptor, Receptor_TF, TF_TarGene networks between the central cell and neighboring cells according to scRNA-Seq expression matrix and barcode table
#'
#' @param GCMat scRNA-seq data. The gene expression matrix（raw) with rows as genes (gene symbols) and columns as cells.
#' @param BarCluFile The annotation results for clustering. The first column is barcode and the second is cell type.
#' @param RecClu character: The central cell
#' @param LigClu character: The neighboring cell
#' @param pval Screening threshold for getting high expressed gene in cluster. The default setting is 0.05.
#' @param logfc Screening threshold for getting high expressed gene in cluster. The default setting is 0.15.
#' @param LigRecLib The file path of Ligand-Receptor interactions. The default setting is 'LigRec.txt' file in the /databases/ folder.
#' @param TFTarLib The file path of TF-TarFget gene interactions. The default setting is 'TFTargetGene.txt' file in the /databases/ folder.
#' @param RecTFLib The file path of Receptor-TF interactions. The default setting is 'RecTF.txt' file in the /databases/ folder.
#'
#' @import Seurat 
#' @import parallel 
#' @importFrom methods as
#' @importFrom stats cor.test t.test
#' @importFrom utils read.table
#'
#' @export
#'
#' @return
#' A list consists of the Ligand_Receptor, Receptor_TF and TF_TarGene signaling subnetwork.
#' The signaling subnetwork is returned as a dataframe object, including three columns:
#' the first column is molecule A, and the second column is molecule B.
#' There is an interaction between the molecules A and B, which correspond to the ligand, receptor,
#' transcription factor or target gene. The third column is used to visualize the signaling subnetwork.
#'
RunMLnet <- function(GCMat, BarCluFile, RecClu, LigClu,
                     pval = 0.05, logfc = 0.15, 
                     LigRecLib = "./database/LigRec.txt",
                     TFTarLib = "./database/TFTargetGene.txt",
                     RecTFLib = "./database/RecTF.txt")
{
  
  #database
  LigRecFile <- LigRecLib
  TFTableFile <- TFTarLib
  RecTFTableFile <- RecTFLib
  
  #get cores
  cores_num <- detectCores()
  cores = ifelse(cores_num>1,cores_num-1,1)
  
  #remove cell not in GCMat
  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  allCell <- colnames(GCMat)
  tableCell <- rownames(BarCluFile)
  print("check table cell:")
  BarCluTable <- BarCluTable[which(BarCluTable$Barcode %in% allCell),]
  print(dim(BarCluTable))
  
  #view some infomation
  print("Rec Cluster:")
  print(RecClu)
  print("Lig Cluster:")
  print(LigClu)
  print("p val:")
  print(pval)
  print("logfc:")
  print(logfc)
  
  #getHighExpGene
  RecClus <- getHighExpGene(GCMat,BarCluTable,RecClu,LigClu,pval,logfc,cores)
  LigClus <- getHighExpGene(GCMat,BarCluTable,LigClu,RecClu,pval,logfc,cores)
  
  #net - Lig_Rec List
  LigRecList <- getLigRec(LigRecFile,LigClus,RecClus)
  stri <- paste("Lig_Rec Num:",length(LigRecList),sep="")
  print(stri)
  
  #net - TF_Target List
  TotalGene <- rownames(GCMat)
  TFTargetTableResult <- TFTargetFisher(TFTableFile,RecClus,TotalGene)
  TFTargetTable <- TFTargetTableResult[[1]]
  TFPval <- TFTargetTableResult[[2]]
  TFP <- as.vector(TFPval$TF)
  TFTarOutSta <- TFTargetTable[which(TFTargetTable$TF %in% TFP),]
  stri <- paste("TF_Target Num:",nrow(TFTarOutSta),sep="")
  print(stri)
  
  #net - Receptor_TF List
  RecTFFisherResult <- RecTFFisher(RecTFTableFile,TFP)
  RecTFFisherTable <- RecTFFisherResult[[1]]
  XZRecXZTFList <- RecTFFisherResult[[2]]
  RecofRecTF <- as.vector(RecTFFisherTable$Receptor)
  stri <- paste("XZRec_XZTF Num:",length(XZRecXZTFList),sep="")
  print(stri)
  
  #get Hign Exp XZ Rec_XZ TF List
  RecofLigRecList <- getPairGeneList(LigRecList,"B")
  RecCommLigTF <- intersect(RecofLigRecList,RecofRecTF)
  stri <- paste("Rec common in LigRec and RecTF:",sep="")
  print(stri)
  if(length(RecCommLigTF) == 0){
    stop("error:LigRec and RecTF don't have the same Rec!")
  }
  
  #update Lig_Rec List
  LigRecList2 <- SearchInPair(LigRecList,RecCommLigTF,"B")
  
  #update Rec_TF List
  HignXZRecXZTF <- SearchInPair(XZRecXZTFList,RecCommLigTF,"A")
  
  #get TF Of HignXZRecXZTF
  TFOfHignXZRecXZTF <- getPairGeneList(HignXZRecXZTF,"B")
  TFCommRecTar <- intersect(TFOfHignXZRecXZTF,TFTargetTable$TF)
  stri <- paste("TF common in RecTF and TFTar:",sep="")
  print(stri)
  if(length(TFCommRecTar) == 0){
    stop("error:RecTF and TFTar don't have the same TF!")
  }
  
  #update TF_Tar List
  TFTarOutSta2 <- TFTargetTable[which(TFTargetTable$TF %in% TFCommRecTar),]
  TFTarF <- as.vector(TFTarOutSta2$key)
  
  #get RecClu barcode
  barListResult <- getBarList(RecClu,LigClu,BarCluTable)
  cells <- barListResult[[1]]
  
  #calculate Cor between RecTF
  print("calculate Cor between RecTF")
  CorRecTFResult <- getCorR(GCMat,HignXZRecXZTF,cells)
  if(length(CorRecTFResult) == 0){
    stop("error:There was no significant correlation between receptors and TF")
  }
  TFofcorRecTF <- getPairGeneList(CorRecTFResult,"B")
  RecofcorRecTF <- getPairGeneList(CorRecTFResult,"A")
  
  #calculate Cor between TFTar
  print("calculate Cor between TFTar")
  CorTFTarResult <- getCorR(GCMat,TFTarF,cells)
  if(length(CorTFTarResult) == 0){
    stop("error:There was no significant correlation between TF and TarGene")
  }
  TFofcorTFTar <- getPairGeneList(CorTFTarResult,"A")
  
  #Final net
  TFofRecTFTar <- intersect(TFofcorRecTF,TFofcorTFTar)
  FinalRecTF <- SearchInPair(CorRecTFResult,TFofRecTFTar,"B")
  FinalTFTar <- SearchInPair(CorTFTarResult,TFofRecTFTar,"A")
  RecoffinalRecTF <- getPairGeneList(FinalRecTF,"A")
  RecofLigRec <- getPairGeneList(LigRecList2,"B")
  RecofLigTF <- intersect(RecofLigRec,RecofcorRecTF)
  FinalLigRec <- SearchInPair(LigRecList2,RecofLigTF,"B")
  
  result <- list("LigRec" = FinalLigRec,
                 "RecTF" = FinalRecTF,
                 "TFTar" = FinalTFTar)
  return(result)
}

