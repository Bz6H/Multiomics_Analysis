#This script is used to annotate and analyse a clean subset of clusters that are 
#Found by using both ATAC and RNA data as parameters for clustering, from the 
#Multiomic_script.
#The UMAP produced in this workflow will be used to find PAX8 and its targets'
#expression pattern along the developmental trajectory

library(ArchR)
library(presto)
library(dplyr)
library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(viridis)
library(Hmisc)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(pheatmap)

set.seed(1)                          # Makes results reproducible
addArchRThreads(threads = 4)          # Adjust to your Mac's CPU count
addArchRGenome('hg38')               # Human genome

library(future)
plan('multisession', workers = 4)  

proj_sub <- loadArchRProject(path =  "/Users/yuzhihuang/MultiomeAnalysis/ArchRSubset_cleaner")

names(proj_sub@reducedDims)
names(proj_sub@embeddings)
names(proj_sub@cellColData)



plotEmbedding(proj_sub,
              name      = 'Clusters_test_Combined',  #Produced from addClusters
              embedding = 'UMAP_Combined',            #Produced from addUMAP
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')
