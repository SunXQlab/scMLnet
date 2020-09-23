# Tutorial of scMLnet
Compiled: September 15, 2020

For this tutorial, we will be using scMLnet to construct the multi-layer signaling network between B cells and Secretory cells from scRNA-Seq data of BALF in COVID-19 patients. The expression matrix and annotation of clstuers can be found in the `/data` folder and the prior information about interactions in the `/database` folder.

# Preparation

## Required packages

We start by loading all required packages. The `Seurat` package is used for normalizing the raw scRNA-Seq data, the `Matrix` package is used for transformation between matrix and sparse matrix and the `parallel` package is used for parallel computation of t.test for screening potentially highly expressed genes.

        library(Seurat)
        library(Matrix)
        library(parallel)
        library(scMLnet)

## Input data

We then read a raw scRNA-Seq data with rows as genes (gene symbols) and columns as cells and the gene expression matrix is required as a sparse matrix. Annotation of cell type consists of two columns: the first column is barcode and the second is cell type. The column number of the gene expression matrix should be consistent with the row number of the annotation table. Numbers can also be used to replace specific cell types in BarCluFile, but with less biological explanation of cell-cell communication between groups.

        # import sample data
        GCMat <- readRDS("./data/data.Rdata")
        GCMat<- as(GCMat,"dgCMatrix")
        
        # import sample annotation
        BarCluFile <- "./data/barcodetype.txt"
        BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)

We next define the receiver cell and sender cell that we want to explorer the cell-cell communication between them. In this example, we focus on the inter-/intracellular signaling network between B cells as senders and Secretory cells as receivers (**NOTE**: make sure the `LigClu` and `RecClu` parameters are the values in the annotation table).

        types <- unique(BarCluTable$Cluster)
        
        LigClu <- "B cells"       #types[4]
        RecClu <- "Secretory"     #types[8]

## Default parameters

Default parameter settings are as follows. Users can provide their own databases which must including three colums: molecule A, molecule B and key (connecting A with B by underlined). Molecules A and B need to have clear identities (i.e., Ligand, Receptor, TF and Target gene), and there should be interactions between them (i.e., Ligand_Receptor, Receptor_TF and TF_Gene interactions).

According to the difference (subjected to `pct` parameter) and ratio (subjected to `logfc` parameter) of gene expression percentage between the sender cells and receiver cells, we firstly define specific highly expressed genes in each of the two types of cells. The highly expressed genes in sender cells are considered as potential ligands and the highly expressed genes in receiver cells are considered as potential receptors and target genes. We screen the potential Ligand-Receptor interactions by searching the Ligands-Receptor pairs in database. We then screen the potential Receptor-TF, TF-Target gene interactions by Fisher's Exact Test (subjected to `pval` parameter).

        pval <- 0.05
        logfc <- 0.15
        LigRecLib <- "./database/LigRec.txt"
        TFTarLib <- "./database/TFTargetGene.txt"
        RecTFLib <- "./database/RecTF.txt"

# Construction of Multi-layer Signaling Networks

Running time of getting highly expressed genes depends on the t-test of gene expression in two types of cells. Doing t-test in parallel can improve the performance of scMLnet (**Precondition**: the parallel package has been installed). The default cores in parallel is set to the number of available logical cores minus one.

        netList <- RunMLnet(GCMat, BarCluFile, RecClu, LigClu, 
                            pval, logfc, 
                            LigRecLib, TFTarLib, RecTFLib)

# Save and Visualization of Multi-layer Signaling Networks

The output `netList` is a list consisiting of gene pairs connecting each upstream layer and downstream layer (i.e., Ligand_Receptor, Receptor_TF and TF_Gene subnetworks). The signaling subnetwork is returned as a dataframe object as the same structure as sample databases. 

        workdir <- "sample"
        DrawMLnet(netList,LigClu,RecClu,workdir,plotMLnet = F)


If you want to visualize the tabular results of multi-layer signaling network, make sure the correct installation of following three python packages: networkx, matplotlib and pymnet and the correct setting of PyHome parameter in `DrawMLnet` funtion. The tabular and graphic output of the constructed multi-layer network will be saved in the `/output/sample/B_Hepatocytes` folder.

        workdir <- "sample"
        PyHome <- "D:/Miniconda3/envs/R36/python.exe" #for Window
        DrawMLnet(netList,LigClu,RecClu,workdir,PyHome,plotMLnet = T)


<div align=center>
<img src="https://github.com/YUZIXD/scMLnet/blob/master/figure/demo2.png" width="400" height="350" alt="demo2" />
</div>

---

# Construction of Multi-cellular Multi-layer Signaling Networks

In the above section, we only focus one pairs of cell type every time, however we can always remain the same receiver cells and only change sender cells so as to construct the multi-cellular multi-layer signaling networks of receiver cells (the central cells) affected by sender cells (neighbor cells). 

        # import data
        GCMat <- readRDS("./data/data.Rdata")
        GCMat<- as(GCMat,"dgCMatrix")
        
        # import annotation
        BarCluFile <- "./data/barcodetype.txt"
        BarCluTable <- read.table(BarCluFile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
        
        ## get LigClu
        LigClus <- unique(BarCluTable$Cluster)
        LigClus <- LigClus[-grep("Secretory|Doublets",LigClus)]
        
        ## creat MLnet
        netList <- list()
        for(ligclu in LigClus){
          
          #sender cell and receiver cell
          LigClu <- ligclu
          RecClu <- "Secretory"
          
          name <- paste(strsplit(LigClu,split = "\\W")[[1]][1],RecClu,sep = "_")
          
          #main
          netList[[name]] <- RunMLnet(GCMat,BarCluFile,RecClu,LigClu)
          
        }
        
        ## save output and plot MLnet
        workdir <- "multi-cellular"
        for (name in names(netList)) {
          
          #scMLnet output
          MLnetList <- netList[[name]]
          print(paste0(name,":"))
          
          #sender cell and receiver cell
          LigClu <- strsplit(name,"_")[[1]][1]
          RecClu <- strsplit(name,"_")[[1]][2]
          
          #main
          PyHome <- "D:/Miniconda3/envs/R36/python.exe" #for Window
          DrawMLnet(MLnetList,LigClu,RecClu,workdir,PyHome,plotMLnet = T)
          
        }

The demo of multi-cellular-mediated ACE2 regulation based on the scRNA-seq data of bronchoalveolar lavage fluid (BALF) in COVID-19 patients:

<div align=center>
<img src="https://github.com/YUZIXD/scMLnet/blob/master/figure/demo1.png" width="600" height="350" alt="demo1" />
</div>





















## Introduction

scMLnet is an R package developed to construct inter-/intracellular multilayer singaling network based on single-cell RNA-seq expression data. scMLnet constructs the multilayer network by integrating intercellular pathways (ligand-receptor interactions) and intracellular subnetworks (receptor-TF pathways and TF-target gene interactions) based on cell-type specific gene expression, prior network information and statistical inference. scMLnet can also visualize the constructed inter-/intracellular signaling pathways between the central cell and neighboring cells. scMLnet is implemented using R
(version 3.6.0) and Python (version 3.7).



The main steps of the scMLnet algorithm include:

* **Step1 Constructing Ligand-Receptor subnetwork**: defines potential Ligand-Receptor subnetworks from scRNA-Seq data and Ligand-Receptor database by getting highly expressed genes (HEGs). HEGs in Type A (sender cells) are considered as potential ligands and HEGs in Type B (receiver cells) as potential receptors.
* **Step2 Constructing TF-Target gene subnetwork**: defines potient TF-Target gene subnetworks from scRNA-Seq data and TF-Target gene database by getting HEG and Fisher’s exact test. HEGs in Type B are considered as potential target genes. Activated TFs can be inferred from the TF-Target gene subnetwork.
* **Step3 Constructing Receptor-TF subnetwork**: defines potient Receptor-TF subnetworks from activated TFs and Receptor-TF database by Fisher’s exact test. Activated receptors can be inferred from the TF-Target gene subnetwork.
* **Step4 constructing multi-layer signaling network**: defines multi-layer signaling network by overlapping the Ligand-Receptor, TF-Target gene, Receptor-TF subnetworks according to common receptors and TFs.

<div align=center>
<img src="./figure/illustration.png" width="800" height="564" alt="illustration" />
</div>


## Installation

### 1.Preparation 

The python package pymnet module should be installed for visualizing the multi-layer singal network, you can download the source files directly from <a href="http://bitbucket.org/bolozna/multilayer-networks-library" target="_blank">bitbucket</a> or <a href="https://github.com/bolozna/Multilayer-networks-library" target="_blank">github</a>, unzip and install:

       python setup.py install

and then copy the /pymnet/sampling folder to the installation path (%PYTHONHOME%/Lib/site-packages/pymnet-0.1-py3.7.egg/pymnet). Alternatively, you can simply copy the /pymnet/ folder to your python library path (%PYTHONHOME%/Lib/site-packages). **NOTE**: The following python packages should be installed first: networkx, matplotlib.

The following R packages should be installed for creating the multi-layer singal network, please install before the installation of scMLnet:
    
* Seurat
* Matrix
* parallel   
    
### 2.Installation

install scMLnet from github:

       install.packages("devtools")
       library(devtools)
       devtools::install_github("YUZIXD/scMLnet")
       library(scMLnet)
    
or install scMLnet module from the source code:

       install.packages("path/to/download/scMLnet_0.1.0.tar.gz", repos = NULL, type = "source")
       library(scMLnet)


## Input and output

### 1.Input

Before using scMLnet, scRNA-seq data should be processed and clustered to identify cell types for dissecting cell type-specific gene expressions by employing existing methods or tools (e.g., Seurat). scMLnet requires the following information as input:

(1) scRNA-Seq expression matrix (a Sparse matrix, where rows represent genes, columns represent cells);

(2) clustering results containing two columns: cell’s barcode and cluster identities;

(3) two cluster identities of sending cells and receiving cells.

### 2.Output

The output of scMLnet has two forms:

(1) tabular information of the constructed multilayer network, containing gene pairs connecting each upstream layer and downstream layer (i.e., Ligand_Receptor, Receptor_TF and TF_Gene subnetworks);

(2) visualization of the constructed multilayer network using the pymnet library.

<div align=center>
<img src="./figure/demo2.png" width="400" height="350" alt="demo" />
</div>

## Working directory structure

The Working Directory requires the following files and directories:

File|Description
---|---
/data|Input directory including scRNA-Seq data and clustering results
/database|Prior information about interactions between ligands, receptors, TFs and target genes
/database/LigRec.txt|The Ligand-Receptor interactions including three colums: Ligand, Receptor and Key (links connecting ligands with receptors by underlined)
/database/RecTF.txt|The Receptor-TF interactions including three colums: Receptor, TF and Key (links connecting receptors with TFs by underlined)
/database/TFTargetGene.txt|The TF-Target gene interactions including three colums: TF, Target gene and Key (links connecting TFs with Target genes by underlined)
/code/DrawNetNew.py|The python script for visualization
/output|Output directory, tabular and graphic results of the multi-layer signaling network will be saved in this folder


## Demonstration

A demonstration of using scMLnet to construct the multi-layer signaling network between B cells and Secretory cells from scRNA-Seq data of COVID-19 patients BALF can be found at following vignette. The expression matrix and annotation of clstuers can be found in the  /data folder and the prior information about interactions in the /database folder.

* Vignette: <a href="./vignettes/Tutorial_of_scMLnet.html" target="_blank">Tutorial of scMLnet</a>



