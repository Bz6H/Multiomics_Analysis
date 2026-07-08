#Dont download the data using Safari as they automatically decompress the tsv.gz
#files and truncate them, use chrome, or untick the automatic opening option in 
# general > safari preference setting.
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


set.seed(1)                          # Makes results reproducible
addArchRThreads(threads = 4)          # Adjust to your Mac's CPU count
addArchRGenome('hg38')               # Human genome

library(future)
plan('multisession', workers = 4)    # Parallelise where possible

#Check current working directory and files
getwd()
list.files(recursive = TRUE)

# Point to fragment files
day10_frags <- "3_Day_10/d10_atac_fragments.tsv.gz"
names(day10_frags) <- "3_Day_10"

day12_frags <- "3_Day_12/d12_atac_fragments.tsv.gz"
names(day12_frags) <- "3_Day_12"

day14_frags <- "3_Day_14/d14_atac_fragments.tsv.gz"
names(day14_frags) <- "3_Day_14"

#Create Arrow files in ArchR
#This step takes long
createArrowFiles(day10_frags, force = TRUE)
createArrowFiles(day12_frags, force = TRUE)
createArrowFiles(day14_frags, force = TRUE)

#Step 2.3, create ArchR project combining 3 arrow files, and integrate RNA
#expression data
proj <- ArchRProject(
  ArrowFiles = c("3_Day_10.arrow", "3_Day_12.arrow", "3_Day_14.arrow"),
  outputDirectory = "ArchRProject",
  copyArrows = TRUE
)
#check the project
proj

#Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRProject", load = FALSE)

#Now loading RNA expression matrix data into the ArchR project:
arrowlist <- list.files(pattern = '*.arrow')
sample_names <- c("3_Day_10","3_Day_12","3_Day_14")
rna_files <- c('3_Day_10/d10_filtered_feature_bc_matrix.h5',
               '3_Day_12/d12_filtered_feature_bc_matrix.h5',
               '3_Day_14/d14_filtered_feature_bc_matrix.h5')
seRNA <- import10xFeatureMatrix(input = rna_files, names = sample_names)

# --- Add RNA to the project, takes a while ---
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA)

# --- Check how many cells you have before filtering ---
proj  # prints cell count
#Cell count = 3473 in this batch 2.
#Cell count = 50060 in batch 3 at this step

#2.4 — Filter Low-Quality Cells
#This step removes cells that are likely to be empty droplets, doublets, dead 
#cells, or poor-quality captures. The thresholds below match those in the 
#original code but you may need to adjust them based on your
#data quality metrics.

proj <- proj[
  !is.na(proj@cellColData$Gex_nGenes) &   # Must have RNA data
    proj$TSSEnrichment > 6           &       # Good ATAC signal at TSSs
    proj$nFrags        > 2500        &       # Enough ATAC fragments
    proj$Gex_nGenes    > 1000        &       # Enough detected genes
    proj$Gex_nUMI      > 5000        &       # Enough RNA molecules
    proj$Gex_nGenes    < 10000               # Remove likely doublets
]
proj
#Cell count = 2857 in this case
#Cell count = 10645
saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRProject", load = FALSE)

#SECTION 3
#Dimensionality Reduction & Clustering
#LSI, Harmony batch correction, UMAP, and cell type annotation

#This section takes you from a filtered project to labelled cell clusters on 
#a UMAP. The key steps are: 
#(1) reduce dimensions separately for ATAC and RNA using LSI
#(2) apply Harmony to remove any batch effects between samples
#(3) cluster cells
#(4) annotate clusters with cell type labels.
proj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject")
proj
# --- LSI on ATAC (TileMatrix = genome tiled into 500bp bins) ---
proj <- addIterativeLSI(
  ArchRProj    = proj,
  useMatrix    = 'TileMatrix',
  depthCol     = 'nFrags',
  name         = 'LSI_ATAC',
  clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
  saveIterations = FALSE,
  force = TRUE)

# --- LSI on RNA ---
proj <- addIterativeLSI(
  ArchRProj    = proj,
  useMatrix    = 'GeneExpressionMatrix',
  depthCol     = 'Gex_nUMI',
  varFeatures  = 2500,
  firstSelection = 'variable',
  binarize     = FALSE,
  name         = 'LSI_RNA',
  clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
  saveIterations = FALSE,
  force = TRUE)

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined",
                        dimWeights = c(5,5))
#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)


#Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRProject", load = FALSE)
proj
proj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject")

#Harmony Batch Correction
#If your samples were processed on different days or batches, 
#technical variation can cluster cells by sample rather than biology. 

proj <- addHarmony(
  ArchRProj   = proj,
  reducedDims = 'LSI_ATAC',
  name        = 'Harmony_ATAC',
  groupBy     = 'Sample',    # Corrects for sample-level batch effects
  theta       = c(2),        # Diversity penalty - increase if samples over-mix
  lambda      = c(1),        # Ridge regression - increase for stronger correction
  do_pca      = TRUE,
  scaleDims   = TRUE,
  corCutOff   = 0.75,
  force       = TRUE)

# --- Build UMAP ---

proj <- addClusters(
  input       = proj,
  reducedDims = "Harmony_ATAC",
  name        = "Clusters_broad",
  resolution  = 0.1,
  force       = TRUE)
table(proj$Clusters_broad)   # check cluster sizes before proceeding

proj <- addUMAP(
  ArchRProj   = proj,
  reducedDims = 'Harmony_ATAC',
  name        = 'UMAP_Harmony',
  minDist     = 0.8,
  dims        = c(1:30),    # Try 1:20 if UMAP looks noisy
  force       = TRUE)

# Make pseudobulk replicates, then Call Peaks
proj <- addGroupCoverages(ArchRProj = proj, 
                          groupBy = "Clusters_broad",
                          minCells = 80,
                          maxCells = 500,
                          minReplicates = 3,
                          maxReplicates = 5,
                          sampleRatio = 0.8,
                          force = TRUE,
                          returnGroups = F)

proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Clusters_broad", 
  pathToMacs2 = "/Users/yuzhihuang/miniforge3/envs/macs2_env/bin/macs2")

#add a new matrix to Multiproj containing insertion counts within our new merged peak set
ps <- proj@peakSet
ps$name <- paste0(seqnames(ps),"_peak",ps$idx)
proj@peakSet <- ps
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE, threads = 16)

# Verify and save
length(proj@peakSet)

saveRDS(
  proj@peakSet,
  "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject/peakSet.rds"
)
saveArchRProject(
  ArchRProj       = proj,
  outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject",
  load            = FALSE
)

proj <- addImputeWeights(proj)

# ── SAVE 
saveArchRProject(ArchRProj = proj, outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject", load = FALSE)

proj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject")
# --- Plot the UMAP coloured by cluster ---

table(proj$Clusters_broad)   # check cluster sizes before proceeding
#C1   C2   C3   C4   C5   C6 
#1592 1675 4367  281 2276  454 
plotEmbedding(proj,
              name      = 'Clusters_broad',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')

# --- Also plot coloured by sample to check batch correction ---
plotEmbedding(proj,
              name      = 'Sample',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1)



# --- Check marker gene expression on UMAP to guide annotation ---
install.packages("hexbin")
# Replace with markers appropriate for your biological system
marker_genes <- c("PAX8", "CDH6", "SOX9", "COL4A1", "ITGB1", "LAMB1")
plot_list <- plotEmbedding(
  ArchRProj  = proj,
  colorBy    = 'GeneExpressionMatrix',
  name       = marker_genes,
  embedding  = 'UMAP_Harmony',
  imputeWeights = getImputeWeights(proj))

combined <- plot_grid(
  plotlist = lapply(plot_list, function(p) {
    p + theme(legend.position = "right",
              plot.title      = element_text(size = 10, hjust = 0.5),
              axis.text       = element_blank(),
              axis.ticks      = element_blank())
  }),
  ncol   = 3,    # 3 columns → 2 rows of 3
  labels = "AUTO")

combined