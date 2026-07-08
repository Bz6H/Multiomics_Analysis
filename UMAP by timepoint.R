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

set.seed(1)                          # Makes results reproducible
addArchRThreads(threads = 4)          # Adjust to your Mac's CPU count
addArchRGenome('hg38')               # Human genome

library(future)
plan('multisession', workers = 4)    # Parallelise where possible

#Check current working directory and files
getwd()
list.files(recursive = TRUE)

proj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject")

#Plotting UMAP by timepoint samples
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
genes_of_interest <- c(
  "PAX8",
  "SOX9",
  "HNF1B",
  "WT1",
  "DCDC2",
  "HOOK1",
  "HOOK2",
  "HOOK3",
  "MAP7",
  "DCX",
  "DCLK1",
  "DCLK2",
  "DCLK3",
  "DCDC1"
) #TMEM173 is old name for

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

# ── Marker genes per cluster (RNA) ────────────────────────────────────────────
markersRNA <- getMarkerFeatures(
  ArchRProj  = proj_day,
  useMatrix  = "GeneExpressionMatrix",
  groupBy    = clusters,          # "Clusters_D10"
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# ── Get top DE genes per cluster ──────────────────────────────────────────────
markerList <- getMarkers(markersRNA, 
                         cutOff = "FDR <= 0.05 & Log2FC >= 1")

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