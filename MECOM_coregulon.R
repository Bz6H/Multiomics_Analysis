## ===== Peak-to-Gene Links ===== ##

library(ArchR)
library(dplyr)
library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1)
addArchRThreads(threads = 4)
addArchRGenome("hg38")

## ── 1. Load project ───────────────────────────────────────────────────────────
proj <- loadArchRProject("/Users/yuzhihuang/MultiomeAnalysis/ArchRProject")
proj@peakSet <- readRDS("/Users/yuzhihuang/MultiomeAnalysis/ArchRProject/peakSet.rds")

proj@peakSet <- proj@peakSet
is.null(proj@peakSet)  # should now be FALSE
length(proj@peakSet)   # should match parent
## ── 2. Compute peak-to-gene links ─────────────────────────────────────────────
proj <- addPeak2GeneLinks(
  ArchRProj       = proj,
  reducedDims     = "Harmony_ATAC",
  useMatrix       = "GeneExpressionMatrix",
  cellsToUse      = NULL,
  addEmpiricalPval = TRUE,
  addPermutedPval  = TRUE,
  nperm            = 10,
  maxDist          = 500000
)

## ── 3. Build annotated peak-to-gene dataframe ─────────────────────────────────
p2geneDF             <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName    <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName    <- (metadata(p2geneDF)$peakSet %>%
                           {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]

annot                    <- readRDS(metadata(p2geneDF)$seATAC)
p2geneDF$peakType        <- annot@rowRanges$peakType[p2geneDF$idxATAC]
p2geneDF$nearestGene     <- annot@rowRanges$nearestGene[p2geneDF$idxATAC]
p2geneDF$GroupReplicate  <- annot@rowRanges$GroupReplicate[p2geneDF$idxATAC]

p2geneDF.peaks           <- as.data.frame(metadata(p2geneDF)[[1]])
p2geneDF.peaks$idxATAC   <- rownames(p2geneDF.peaks)
p2geneDF.merged          <- merge(p2geneDF, p2geneDF.peaks, by = "idxATAC")

p2g.df.obs               <- as.data.frame(p2geneDF.merged)
p2g.df.obs               <- p2g.df.obs[complete.cases(p2g.df.obs), ]
p2g.df.obs$Unique_Peak_ID <- rownames(p2g.df.obs)
p2g.df.obs$strand        <- "+"

names(p2g.df.obs)[names(p2g.df.obs) == "seqnames"] <- "chromosome"
names(p2g.df.obs)[names(p2g.df.obs) == "start"]    <- "starting position"
names(p2g.df.obs)[names(p2g.df.obs) == "end"]      <- "ending position"

p2g.df.obs               <- p2g.df.obs[, c(19, 14, 15, 16, 18, 1, 2, 3, 4,
                                           5, 6, 7, 8, 9, 10, 11, 12, 13, 17)]
p2g.df.obs$GroupReplicate <- str_split_fixed(p2g.df.obs$GroupReplicate, "\\.", n = Inf)[, 1]

## ── 4. Filter to significant links ────────────────────────────────────────────
p2g.df.obs.sub <- p2g.df.obs %>%
  dplyr::filter(FDR < 1e-4, Correlation > 0.4)

cat("Total significant peak-gene links:", nrow(p2g.df.obs.sub), "\n")


###########################
#HOMER TF-SPECIFIC MOTIF SEARCH DONE IN BASH (motif_finding.sh)
###########################

library(ArchR)
library(dplyr)
library(tidyverse)

projPath <- "/Volumes/Elements/Leo_CUT_RUN"

p2g.df.obs.sub <- read.delim(
  file = "/Users/yuzhihuang/MultiomeAnalysis/p2g.df.obs.sub.txt")
cat("Total significant peak-gene links:", nrow(p2g.df.obs.sub), "\n")

setwd("/Users/yuzhihuang/MECOM_predicted_sites")

homer_out <- list.files(pattern = "*.txt") %>%
  lapply(read.delim) %>%
  bind_rows()

cat("Total HOMER motif hits:", nrow(homer_out), "\n")
cat("Unique peaks with MECOM motif:", length(unique(homer_out$PositionID)), "\n")

links <- merge(
  p2g.df.obs.sub,
  homer_out,
  by.x = "Unique_Peak_ID",
  by.y = "PositionID"
)

cat("Peaks after merge:", nrow(links), "\n")

links <- links[, c("Unique_Peak_ID", "geneName", "Correlation",
                   "Motif.Name", "peakType", "peakName")]

links$Motif.Name <- toupper(
  str_split_fixed(links$Motif.Name, "\\(", n = Inf)[, 1]
)

links <- links %>%
  dplyr::filter(!duplicated(Unique_Peak_ID)) %>%
  dplyr::arrange(geneName)

cat("MECOM-linked peaks:", nrow(links), "\n")
cat("Unique predicted MECOM target genes:", length(unique(links$geneName)), "\n")
table(links$peakType)

write.csv(
  links,
  file      = "/Users/yuzhihuang/MECOM_predicted_sites/Predicted_MECOM_regulon.csv",
  row.names = FALSE
)

cat("Saved to: /Users/yuzhihuang/MECOM_predicted_sites/Predicted_MECOM_regulon.csv\n")


peaks_multi <- read.csv("/Users/yuzhihuang/MECOM_predicted_sites/Predicted_MECOM_regulon.csv")
peaks_multi <- peaks_multi[, c("peakName", "geneName", "peakType")]
peaks_multi <- peaks_multi %>%
  separate(peakName, into = c("chr", "pos"), sep = ":") %>%
  separate(pos, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.integer(start),
         end   = as.integer(end))

nrow(peaks_multi)

MECOM_genes <- unique(peaks_multi$geneName)
MECOM_genes <- MECOM_genes[!is.na(MECOM_genes)]

PAX8_genes <- read.table(
  paste0(projPath, "/Validated/All_PAX8_val_genes.txt"),
  header = FALSE
)$V1

PAX8_MECOM_coreg <- intersect(PAX8_genes, MECOM_genes)
length(PAX8_MECOM_coreg)

dir.create(paste0(projPath, "/Validated/MECOM"), recursive = TRUE, showWarnings = FALSE)

write.table(
  PAX8_MECOM_coreg,
  file = paste0(projPath, "/Validated/MECOM/PAX8_MECOM_coregulated_genes.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)


library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

dir.create(paste0(projPath, "/Output/MECOM_PAX8_coregulon"), recursive = TRUE, showWarnings = FALSE)

Timepoints <- c("All")

for (day in Timepoints) {
  
  cat("Processing", day, "...\n")
  
  gene_file <- paste0(projPath, "/Validated/MECOM/PAX8_MECOM_coregulated_genes.txt")
  genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  cat("  Loaded", length(genes), "genes\n")
  
  go_results_BP <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_results_MF <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_results_CC <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_simple_BP <- simplify(go_results_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_simple_MF <- simplify(go_results_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_simple_CC <- simplify(go_results_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  pdf(paste0(projPath, "/Output/MECOM_PAX8_coregulon/", day, "_MECOM+PAX8_GOBP_Dotplot.pdf"),
      width = 16, height = 12)
  print(dotplot(go_simple_BP, showCategory = 30) +
          ggtitle(paste(day, "MECOM/PAX8 Target Genes: Biological Processes")))
  dev.off()
  
  pdf(paste0(projPath, "/Output/MECOM_PAX8_coregulon/", day, "_MECOM+PAX8_GOMF_Dotplot.pdf"),
      width = 16, height = 12)
  print(dotplot(go_simple_MF, showCategory = 30) +
          ggtitle(paste(day, "MECOM/PAX8 Target Genes: Molecular Functions")))
  dev.off()
  
  pdf(paste0(projPath, "/Output/MECOM_PAX8_coregulon/", day, "_MECOM+PAX8_GOCC_Dotplot.pdf"),
      width = 16, height = 8)
  print(dotplot(go_simple_CC, showCategory = 20) +
          ggtitle(paste(day, "MECOM/PAX8 Target Genes: Cellular Components")))
  dev.off()
  
  write.csv(as.data.frame(go_simple_BP),
            paste0(projPath, "/Output/MECOM_PAX8_coregulon/", day, "_MECOM+PAX8_GOBP_Results.csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_simple_MF),
            paste0(projPath, "/Output/MECOM_PAX8_coregulon/", day, "_MECOM+PAX8_GOMF_Results.csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_simple_CC),
            paste0(projPath, "/Output/MECOM_PAX8_coregulon/", day, "_MECOM+PAX8_GOCC_Results.csv"),
            row.names = FALSE)
  
  cat("  Completed", day, "\n\n")
}

cat("All timepoints processed!\n")


gene_file <- paste0(projPath, "/Validated/MECOM/PAX8_MECOM_coregulated_genes.txt")
genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)$V1

go_res <- enrichGO(
  gene         = genes,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pvalueCutoff = 0.01,
  readable     = TRUE
)
go_res_simple <- simplify(go_res, cutoff = 0.7)
go_res_simple <- pairwise_termsim(go_res_simple)

emapplot(
  go_res_simple,
  showCategory   = 10,
  size_category  = 1.5,
  node_label_size = 4,
  nCluster       = 3
)

cnetplot(
  go_res_simple,
  layout         = "circle",
  color_category = "#a50f15",
  size_category  = 0.8,
  color          = "category",
  node_label     = "share",
  color_item     = "#fee5d9",
  size_item      = 0.5,
  showCategory   = 10,
  hilight        = "none",
  hilight_alpha  = 0.3,
  curvature      = 0.1,
  size_edge      = 0.8
)
