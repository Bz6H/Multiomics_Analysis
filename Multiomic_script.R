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
library(enrichplot
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(pheatmap)

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
proj <- addUMAP(proj, 
                reducedDims = "LSI_Combined", 
                name = "UMAP_Combined", 
                minDist = 0.8, 
                force = TRUE)


#Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRProject", load = FALSE)
proj
proj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject")

names(proj@reducedDims)
names(proj@embeddings)
#Harmony Batch Correction
#If your samples were processed on different days or batches, 
#technical variation can cluster cells by sample rather than biology. 

proj <- addHarmony(
  ArchRProj   = proj,
  reducedDims = 'LSI_Combined',  #editable, RNA/ATAC/Combined
  name        = 'Harmony_Combined',
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
  reducedDims = "Harmony_Combined",       #Editable ATAC/RNA/Combined
  name        = "Clusters_test_Combined",
  resolution  = 0.4,
  force       = TRUE)
table(proj$Clusters_test_Combined) 

# check cluster sizes before proceeding
#broad ATAC: C1   C2   C3   C4   C5   C6 
#           1592 1675 4367  281 2276  454 

proj <- addUMAP(
  ArchRProj   = proj,
  reducedDims = 'Harmony_Combined',
  name        = 'UMAP_Combined',
  minDist     = 0.8,
  dims        = c(1:30),    # Try 1:20 if UMAP looks noisy
  force       = TRUE)

#Plot by cluster to check, change name to plot different cluster resolution
plotEmbedding(proj,
              name      = 'Clusters_test_Combined',  #Produced from addClusters
              embedding = 'UMAP_Combined',            #Produced from addUMAP
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')

# --- Also plot coloured by sample to check batch correction ---
plotEmbedding(proj,
              name      = 'Sample',
              embedding = 'UMAP_RNA',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 0.5)


#Subset to remove some intermediate clusters, make it looks cleaner(test)
keep_clusters <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")

keep_cells <- getCellNames(
  proj[proj$Clusters_test %in% keep_clusters, ]
)

length(keep_cells)


setwd("/Users/yuzhihuang/MultiomeAnalysis")

proj_sub <- subsetArchRProject(
  ArchRProj       = proj,
  cells           = keep_cells,
  outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRSubset_cleaner",
  dropCells       = TRUE,
  logFile         = NULL,
  threads         = getArchRThreads(),
  force           = TRUE
)

plotEmbedding(proj_sub,
              name      = 'Clusters_test',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')
plotEmbedding(proj_sub,
              name      = 'Clusters_broad',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')
# --- Also plot coloured by sample to check batch correction ---
plotEmbedding(proj_sub,
              name      = 'Sample',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 0.5)
proj_sub <- addImputeWeights(proj_sub, reducedDims = "Harmony_ATAC")

#############################################
#Block below used to find significant accessible peaks
#used to identify coactivators and find differential accessible regions between
#cell types
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


# ── SAVE 
saveArchRProject(ArchRProj = proj, outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject", load = FALSE)

proj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject") 

saveArchRProject(ArchRProj = proj_sub, outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRSubset_cleaner", load = FALSE)

##################################################################
######################GENE EXPRESSION ON UMAP#####################
##################################################################
# --- Check marker gene expression on UMAP to guide annotation ---
install.packages("hexbin")
# Replace with markers appropriate for your biological system
marker_genes <- c('CDH6','WT1','JAG1','HNF1B','CDH1','PAX8','LAMC1')

plot_list <- plotEmbedding(
  ArchRProj  = proj,
  colorBy    = 'GeneExpressionMatrix',
  name       = marker_genes,
  embedding  = 'UMAP_ATAC',    #Editable ATAC/RNA/Combined
  plotAs = "points",
  size = 0.5,
  rastr = F,
  continuousSet = "greenBlue",
  imputeWeights = getImputeWeights(proj))

combined <- plot_grid(
  plotlist = lapply(plot_list, function(p) {
    p + theme(
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.title       = element_text(size = 4),
      legend.text        = element_text(size = 4),
      legend.key.size    = unit(0.3, "cm"),
      plot.title         = element_text(size = 4, hjust = 0.5),
      axis.text          = element_blank(),
      axis.ticks         = element_blank(),
      axis.title         = element_text(size = 4),
      plot.margin        = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
  }),
  ncol        = 5,
  labels      = "AUTO",
  label_size  = 6,
  rel_widths  = c(1, 1, 1)
)

combined

#############################################################################
####### Annotating the multiomics clusters by first getting DEG of each cluster###
###### Then compare to the nephrogenesis atlas cluster DEG to see which cluster 
###### The multiomics cluster is most similar to

library(readxl)
library(dplyr)
library(pheatmap)

markersRNA <- getMarkerFeatures(
  ArchRProj  = proj_sub,
  useMatrix  = "GeneExpressionMatrix",
  groupBy    = "Clusters_test",
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersRNA,
                         cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

multiome_weighted <- lapply(markerList, function(x) {
  x <- as.data.frame(x)
  x <- x[!is.na(x$FDR) & x$FDR <= 0.05 & x$Log2FC >= 0.5, ]
  list(
    genes   = x$name,
    weights = -log10(x$FDR + 1e-300)
  )
})
multiome_weighted <- multiome_weighted[sapply(multiome_weighted, function(x) length(x$genes) > 0)]

cat("Multiome clusters with DEGs:\n")
print(sapply(multiome_weighted, function(x) length(x$genes)))

atlas <- read_excel("/Users/yuzhihuang/MultiomeAnalysis/nephrogenesis_atlas_DEG.xlsx")

atlas_weighted <- atlas %>%
  dplyr::filter(!is.na(p_val_adj), p_val_adj < 0.05, avg_logFC > 0.5) %>%
  dplyr::mutate(
    weight  = -log10(p_val_adj + 1e-300),
    cluster = as.character(cluster)
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    data = list(data.frame(gene = gene, weight = weight))
  ) %>%
  tibble::deframe()

atlas_weighted <- atlas_weighted[order(as.numeric(names(atlas_weighted)))]

cat("\nAtlas clusters with DEGs:\n")
print(sapply(atlas_weighted, nrow))

weighted_jaccard <- function(genes_a, weights_a, genes_b, weights_b) {
  if (length(genes_a) == 0 | length(genes_b) == 0) return(0)
  if (length(intersect(genes_a, genes_b)) == 0) return(0)
  
  all_genes <- union(genes_a, genes_b)
  
  wa_full <- rep(0, length(all_genes))
  wb_full <- rep(0, length(all_genes))
  
  idx_a <- match(all_genes, genes_a)
  idx_b <- match(all_genes, genes_b)
  
  wa_full[!is.na(idx_a)] <- weights_a[idx_a[!is.na(idx_a)]]
  wb_full[!is.na(idx_b)] <- weights_b[idx_b[!is.na(idx_b)]]
  
  sum(pmin(wa_full, wb_full)) / sum(pmax(wa_full, wb_full))
}

jaccard_mat <- matrix(
  NA,
  nrow = length(multiome_weighted),
  ncol = length(atlas_weighted),
  dimnames = list(names(multiome_weighted), names(atlas_weighted))
)

for (m in names(multiome_weighted)) {
  for (a in names(atlas_weighted)) {
    jaccard_mat[m, a] <- weighted_jaccard(
      genes_a   = multiome_weighted[[m]]$genes,
      weights_a = multiome_weighted[[m]]$weights,
      genes_b   = atlas_weighted[[a]]$gene,
      weights_b = atlas_weighted[[a]]$weight
    )
  }
}

cat("\nJaccard range:", range(jaccard_mat, na.rm = TRUE), "\n")
print(round(jaccard_mat, 3))

best_matches <- data.frame(
  Multiome_cluster = rownames(jaccard_mat),
  Best_atlas_match = colnames(jaccard_mat)[apply(jaccard_mat, 1, which.max)],
  Jaccard_score    = apply(jaccard_mat, 1, max),
  Shared_genes     = NA_integer_
)

for (i in seq_len(nrow(best_matches))) {
  m <- best_matches$Multiome_cluster[i]
  a <- best_matches$Best_atlas_match[i]
  best_matches$Shared_genes[i] <- length(intersect(
    multiome_weighted[[m]]$genes,
    atlas_weighted[[a]]$gene
  ))
}

best_matches <- best_matches[order(-best_matches$Jaccard_score), ]
print(best_matches)

write.csv(best_matches,
          "/Users/yuzhihuang/MultiomeAnalysis/output200526/multiome_atlas_best_matches.csv",
          row.names = FALSE)

write.csv(jaccard_mat,
          "/Users/yuzhihuang/MultiomeAnalysis/output200526/jaccard_full_matrix.csv")

pdf("/Users/yuzhihuang/MultiomeAnalysis/output200526/jaccard_heatmap.pdf",
    width = 12, height = 8)
pheatmap(
  jaccard_mat,
  cluster_rows    = TRUE,
  cluster_cols    = F,
  display_numbers = TRUE,
  number_format   = "%.3f",
  angle_col       = 0,
  fontsize_number = 7,
  color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
  main  = "Weighted Jaccard (p-adj): multiome vs nephrogenesis atlas"
)
dev.off()

##############################################################
#########QUANTIFY GENE EXPRESSION PER CLUSTER###############
plotGroups(
  ArchRProj  = proj,
  groupBy    = "Clusters_broad",
  colorBy    = "GeneExpressionMatrix",
  name       = marker_genes,
  plotAs     = "violin",
  alpha      = 0.4,
  addBoxPlot = TRUE,
  imputeWeights = getImputeWeights(proj)   # drop this arg if you want raw values
)

se <- getGroupSE(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix",
  groupBy   = "Clusters_broad",
  divideN   = TRUE        # mean per cell, not sum
)
assay(se)[which(rowData(se)$name == "CDH4"),]
##############################################################################
#######Pseudotime Trajectory of epithelial clusters based on ATAC data########
###using this to plot the changes in PAX8 target expression over maturity#####
##############################################################################
proj <- addClusters(
  input       = proj,
  reducedDims = "Harmony_ATAC",
  name        = "Clusters_finer",
  resolution  = 1,
  force       = TRUE)

table(proj$Clusters_fine) #Check cluster size
#C1   C2   C3   C4   C5   C6   C7   C8   C9 
#455  281 1686 1947 1586 1681  714 1815  480 
#JP's code used 0.3 resolution to plot trajectory data, but by increasing resolution to 0.5
#(C2,3,4,7)

#"AddCluster" statistically categorise the data, addUMAP visualises it
#After using plotEmbedding to plot a graph using addUMAP output (UMAP_Harmony in this case),
#Each cell is labelled with their assigned cluster (i.e. cluters_fine)

proj <- addUMAP(
  ArchRProj   = proj,
  reducedDims = 'Harmony_ATAC',
  name        = 'UMAP_Harmony',
  minDist     = 0.8,
  dims        = c(1:30),
  force       = TRUE)

#LSI generated top 30 "principal components" the dims setting here specifies which
#components to use to visualise the UMAP. We are using all 30 in this case, JP used
#top 28 (i.e. dims = c(1:28), idk why)

plotEmbedding(proj,
              name      = 'Clusters_finer',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 0.5,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')

#I was able to see 4 subclusters in D14 sample instead of 3, now check for marker genes

marker_genes <- c()
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

# --- Assign cell type labels based on your marker gene analysis ---
# Modify the cluster-to-celltype mapping to match your data
proj$celltype <- NA
proj$celltype[proj$Clusters_fine == 'C1'] <- 'D12_Stromal'
proj$celltype[proj$Clusters_fine == 'C2'] <- 'D14_Stromal'
proj$celltype[proj$Clusters_fine == 'C3'] <- 'D14_Podocyte' #"MAFB","WT1","FOXC2","HEY1","NPHS1"
proj$celltype[proj$Clusters_fine == 'C4'] <- 'D14_Proximal' #"JAG1","FXYD2", "LRP2"
proj$celltype[proj$Clusters_fine == 'C5'] <- 'D10_Induced'
proj$celltype[proj$Clusters_fine == 'C6'] <- 'D12_Induced'
proj$celltype[proj$Clusters_fine == 'C7'] <- 'D14_Distal'   #"POU3F3","CLDN3"
proj$celltype[proj$Clusters_fine == 'C8'] <- 'D10_Stromal?'
proj$celltype[proj$Clusters_fine == 'C9'] <- 'D12_?'
# Add more lines for each cluster

# --- Plot final annotated UMAP ---
plotEmbedding(proj,
              name      = 'celltype',
              embedding = 'UMAP_Harmony',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')


#JP then assigned the trajectory backbone using the epithelial cluster he identified
#earlier, D10epi, D12 epi, D14epi.
Trajectory <- c("C5", "C6", "C4", "C7") 
#In our case, D10 induced --> D12 induced --> D14 Proximal -- > D14 Distal
Trajectory_pod <-c('C5','C6','C3')

proj <- addTrajectory(
  ArchRProj = proj, 
  name = "Trajectory_pod", 
  groupBy = "Clusters_fine",
  trajectory = Trajectory_pod, 
  embedding = "UMAP_Harmony",
  force = TRUE,
  seed = 1,
  useAll = F,
  spar = 1.5,
  preFilterQuantile = 1,
  postFilterQuantile = 1,
  saveDF = paste0("/Users/yuzhihuang/MultiomeAnalysis/Trajectory_pod.rds")
)

####Plot Trajectory
trajectory_cells <- getCellNames(
  proj[complete.cases(proj@cellColData$Trajectory), ]
)

proj@cellColData$Trajectory_plot <- proj@cellColData$Trajectory
proj@cellColData$Trajectory_plot[is.na(proj@cellColData$Trajectory_plot)] <- 0

plotEmbedding(
  ArchRProj     = proj,
  colorBy       = "cellColData",
  name          = "Trajectory_plot",
  embedding     = "UMAP_Harmony",
  plotAs        = "points",
  size          = 1,
  rastr         = FALSE,
  continuousSet = "fireworks2",
  labelMeans    = FALSE,
  highlightCells = trajectory_cells
) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_blank(),
    axis.text.y      = element_blank(),
    axis.title.x     = element_text(hjust = 0, vjust = 2.5),
    axis.title.y     = element_text(hjust = 0, vjust = -1.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    line             = element_blank(),
    plot.title       = element_text(hjust = 0, vjust = -1.5)
  ) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2")

#Now with the trajectory I want to measure gene expression changes along it#
library(dplyr)

ann_colors.2 <- list(
  Cluster = c("C5" = "#2A1E5C",
              "C6" = "#9D2C24",
              "C4" = "#E25822",
              "C7" = "#F2C140")
)
#Assign colours to the clusters along the pseudotime track
label <- data.frame(
  cellName   = proj$cellNames,
  Pseudotime = proj$Trajectory,
  Cluster    = proj$Clusters_fine
)
#Create dataframe with one row per cell storing cell ID, trajectory order (0-100)
#and the cluster the cell belongs to.

label <- label[complete.cases(label), ]
#Drop cells with NA for trajectory, like stromal, podocytes, etc.

label$bins <- as.numeric(cut(label$Pseudotime, 100, labels = FALSE))
#cuts the trajectory value to 100 bins with 1 in width, each cell gets categorise
#into a bin

label <- dplyr::arrange(label, Pseudotime)
#Sort all cells by the pseudotime bin they get assigned to

label <- label %>%
  dplyr::group_by(bins) %>%
  dplyr::count(Cluster)
#Group cells that are on the same trajecotry point together, count at each point,
#number of cells from each cluster

label <- as.data.frame(label)
#Convert into data frame, with bin name, and number of cells per cluster

label <- label[order(label[, "bins"], -label[, "n"]), ]
# Sort by pseudotime coordinate than by cell number of each cluster at the point

label <- label[!duplicated(label$bins), ]
#Isolate the dominant cluster per pseudotime point
label$n    <- NULL
label$bins <- NULL
rownames(label) <- 1:nrow(label)

trajGEM <- getTrajectory(
  ArchRProj = proj, 
  name = "Trajectory", 
  useMatrix = "GeneExpressionMatrix", 
  log2Norm = TRUE
) 
# Create a matrix of all gene's average expression along the pseudotime bin in the trajectory
grep(":PAX8$", rownames(trajGEM), value = TRUE)
#From the master matrix, extract PAX8 expression along this trajectory
#Plot PAX8 expression along the pseudotrajectory

###########Quality check of marker#################
########by plotting a line graph along pseudo######
markers <- c("PAX8",'LAMC1','LAMC2','COL1A1','COL6A1','COL6A2','EZR','ITGA9','ITGB5')

marker_rows <- unlist(lapply(markers, function(g)
  grep(paste0(":", g, "$"), rownames(trajGEM), value = TRUE)))

marker_mat <- assay(trajGEM)[marker_rows, , drop = FALSE]
rownames(marker_mat) <- sub(".*:", "", rownames(marker_mat))

#marker_mat useful, 
#Generates a matrix with genes as row, and pseudotime axis as column, can be used for correlation.

marker_df <- data.frame(
  bin        = rep(1:100, times = nrow(marker_mat)),
  gene       = rep(rownames(marker_mat), each = 100),
  expression = as.vector(t(marker_mat))
)

ggplot(marker_df, aes(x = bin, y = expression, colour = gene)) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Pseudotime bin", y = "log2 expression",
       title = "Marker gene expression along trajectory") +
  theme_classic(base_size = 12)


#############################################
#######Plot expression overlapping UMAP

traj_exp <- plotTrajectory(
  ArchRProj = proj,
  trajectory = "Trajectory",
  colorBy   = "GeneExpressionMatrix",
  name      = "LAMC1",
  embedding = "UMAP_Harmony",
  continuousSet = "fireworks2"
)

traj_exp[[1]]
traj_exp[[2]]

# ── SAVE 
saveArchRProject(ArchRProj = proj, outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject", load = FALSE)

########################################
########################################


#####Get top expressed genes per cluster
markersRNA <- getMarkerFeatures(
  ArchRProj  = proj_sub,
  useMatrix  = "GeneExpressionMatrix",
  groupBy    = "Clusters_test",
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# ── Get top DE genes per cluster ──────────────────────────────────────────────
markerList <- getMarkers(markersRNA, 
                         cutOff = "FDR <= 0.05 & Log2FC >= 1")

view(markerList$C7)

# Print top 20 per cluster
lapply(markerList, function(x) head(x[order(x$Log2FC, decreasing = TRUE), ], 20))


###Gene ontology of each cluter based on differentially expressed genes
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ── Convert marker list to entrez IDs and run GO per cluster ──────────────────
go_results <- lapply(names(markerList), function(cl) {
  
  # Get significant markers for this cluster
  df <- as.data.frame(markerList[[cl]])
  genes <- df$name[df$FDR <= 0.05 & df$Log2FC >= 1]
  
  cat("\nCluster", cl, "— testing", length(genes), "genes\n")
  if (length(genes) < 5) { 
    cat("Too few genes, skipping\n")
    return(NULL) 
  }
  
  # Convert gene symbols → Entrez IDs
  entrez <- bitr(genes, 
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
  
  # GO Biological Process enrichment
  ego <- enrichGO(gene          = entrez$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",           # BP / MF / CC
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)            # converts back to symbols
  return(ego)
})
names(go_results) <- names(markerList)

# ── Print top 15 GO terms per cluster ─────────────────────────────────────────
lapply(names(go_results), function(cl) {
  ego <- go_results[[cl]]
  if (is.null(ego) || nrow(ego) == 0) {
    cat("Cluster", cl, "— no enriched terms\n")
    return(NULL)
  }
  cat("\n========== Cluster", cl, "— top GO:BP terms ==========\n")
  print(as.data.frame(ego)[1:min(15, nrow(ego)), 
                           c("Description", "GeneRatio", "pvalue", 
                             "p.adjust", "geneID")])
})







#######################################################################
#Plotting UMAP by timepoint samples
#######################################################################
# ── CONFIG: change only this block for each timepoint ─────────────────────────
day        <- "3_Day_12"          # options: "3_Day_10", "3_Day_12", "3_Day_14"
short_day  <- "D12"               # used in naming: "D10", "D12", "D14"
resolution <- 0.2                 # adjust per timepoint if needed

# ── Auto-generated names (don't touch) ────────────────────────────────────────
lsi_atac  <- paste0("LSI_ATAC_",  short_day)
lsi_rna   <- paste0("LSI_RNA_",   short_day)
clusters  <- paste0("Clusters_",  short_day)
umap_name <- paste0("UMAP_",      short_day)
out_dir   <- paste0("ArchRProject_", short_day)

# ── 1. Subset ─────────────────────────────────────────────────────────────────
proj_day <- proj[proj$Sample == day, ]
cat("Cells in", day, ":", nCells(proj_day), "\n")

# ── 2. LSI ATAC ───────────────────────────────────────────────────────────────
proj_day <- addIterativeLSI(
  ArchRProj      = proj_day,
  useMatrix      = "TileMatrix",
  depthCol       = "nFrags",
  name           = lsi_atac,
  clusterParams  = list(resolution = 0.2, sampleCells = 5000, n.start = 10),
  saveIterations = FALSE,
  force          = TRUE)

# ── 3. LSI RNA ────────────────────────────────────────────────────────────────
proj_day <- addIterativeLSI(
  ArchRProj      = proj_day,
  useMatrix      = "GeneExpressionMatrix",
  depthCol       = "Gex_nUMI",
  varFeatures    = 2500,
  firstSelection = "variable",
  binarize       = FALSE,
  name           = lsi_rna,
  clusterParams  = list(resolution = 0.2, sampleCells = 5000, n.start = 10),
  saveIterations = FALSE,
  force          = TRUE)

# ── 4. Cluster ────────────────────────────────────────────────────────────────
proj_day <- addClusters(
  input       = proj_day,
  reducedDims = lsi_atac,
  name        = clusters,
  resolution  = resolution,
  force       = TRUE)

table(getCellColData(proj_day)[[clusters]])

# ── 5. UMAP ───────────────────────────────────────────────────────────────────
proj_day <- addUMAP(
  ArchRProj   = proj_day,
  reducedDims = lsi_atac,
  name        = umap_name,
  minDist     = 0.5,
  dims        = 1:30,
  force       = TRUE)

# ── 6. Imputation weights ─────────────────────────────────────────────────────
proj_day <- addImputeWeights(proj_day, reducedDims = lsi_atac)

# ── 7. Save ───────────────────────────────────────────────────────────────────
saveArchRProject(ArchRProj       = proj_day,
                 outputDirectory = out_dir,
                 load            = FALSE)

# ── 8. UMAP plot ──────────────────────────────────────────────────────────────
plotEmbedding(proj_day,
              embedding      = umap_name,
              colorBy        = "cellColData",
              name           = clusters,
              plotAs         = "points",
              size           = 1,
              labelAsFactors = FALSE,
              labelMeans     = TRUE,
              discreteSet    = "stallion")

# ── Plot key kidney marker genes on UMAP ──────────────────────────────────────
genes_of_interest <- c("PAX8", "PDGFRA","WT1", "CGAS", "TMEM173",
                       "UMOD", "PKD1", "PKD2") #TMEM173 is old name for

plot_list <- plotEmbedding(
  ArchRProj     = proj_day,
  colorBy       = "GeneExpressionMatrix",
  name          = genes_of_interest,
  embedding     = umap_name,
  imputeWeights = getImputeWeights(proj_day),
  plotAs        = "points",
  size          = 0.6,
  returnPlot    = TRUE)   # returns list instead of printing

# ── Side-by-side combined plot ────────────────────────────────────────────────
library(ggplot2)
library(cowplot)

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
################################################################################