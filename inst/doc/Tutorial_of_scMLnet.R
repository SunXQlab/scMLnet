## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  library(Seurat)
#  library(Matrix)
#  library(parallel)
#  library(scMLnet)

## ----eval=FALSE---------------------------------------------------------------
#  # import sample data
#  GCMat <- readRDS("./example/data.Rdata")
#  GCMat<- as(GCMat,"dgCMatrix")
#  
#  # import sample annotation
#  BarCluFile <- "./example/barcodetype.txt"
#  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  types <- unique(BarCluTable$Cluster)
#  
#  LigClu <- "B cells"       #types[4]
#  RecClu <- "Secretory"     #types[8]

## ----eval=FALSE---------------------------------------------------------------
#  pval <- 0.05
#  logfc <- 0.15
#  cores <- NULL
#  LigRecLib <- "./database/LigRec.txt"
#  TFTarLib <- "./database/TFTargetGene.txt"
#  RecTFLib <- "./database/RecTF.txt"

## ----eval=FALSE---------------------------------------------------------------
#  cores <- detectCores()-1
#  netList <- RunMLnet(GCMat, BarCluFile, RecClu, LigClu,
#                      pval, logfc, cores,
#                      LigRecLib, TFTarLib, RecTFLib)

## ----eval=FALSE---------------------------------------------------------------
#  workdir <- "sample"
#  DrawMLnet(netList,LigClu,RecClu,workdir,plotMLnet = F)

## ----eval=FALSE---------------------------------------------------------------
#  workdir <- "sample"
#  PyHome <- "D:/Miniconda3/envs/R36/python.exe" #for Window
#  DrawMLnet(netList,LigClu,RecClu,workdir,PyHome,plotMLnet = T)

## ----eval=FALSE---------------------------------------------------------------
#  
#  # import data
#  GCMat <- readRDS("./example/data.Rdata")
#  GCMat<- as(GCMat,"dgCMatrix")
#  
#  # import annotation
#  BarCluFile <- "./example/barcodetype.txt"
#  BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
#  
#  ## get LigClu
#  LigClus <- unique(BarCluTable$Cluster)
#  LigClus <- LigClus[-grep("Secretory|Doublets",LigClus)]
#  
#  ## get cores
#  cores <- detectCores()-1
#  
#  ## creat MLnet
#  netList <- list()
#  for(ligclu in LigClus){
#  
#    #sender cell and receiver cell
#    LigClu <- ligclu
#    RecClu <- "Secretory"
#  
#    name <- paste(strsplit(LigClu,split = "\\W")[[1]][1],RecClu,sep = "_")
#  
#    #main
#    netList[[name]] <- RunMLnet(GCMat,BarCluFile,RecClu,LigClu,cores = cores)
#  
#  }
#  
#  ## save output and plot MLnet
#  workdir <- "microenvironment"
#  for (name in names(netList)) {
#  
#    #scMLnet output
#    MLnetList <- netList[[name]]
#    print(paste0(name,":"))
#  
#    #sender cell and receiver cell
#    LigClu <- strsplit(name,"_")[[1]][1]
#    RecClu <- strsplit(name,"_")[[1]][2]
#  
#    #main
#    PyHome <- "D:/Miniconda3/envs/R36/python.exe" #for Window
#    DrawMLnet(MLnetList,LigClu,RecClu,workdir,PyHome,plotMLnet = T)
#  
#  }
#  

