## ===== Peak-to-Gene Links for HNF1B ===== ##

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
#HOMER MOTIF ANALYSIS DONE IN BASH
###########################

# Annotate linked peaks with gene names by importing this table into R. (the below code runs in R)

p2g.df.obs.sub <- read.delim(
  file = "/Users/yuzhihuang/MultiomeAnalysis/p2g.df.obs.sub.txt")
cat("Total significant peak-gene links:", nrow(p2g.df.obs.sub), "\n")

##Combine all Homer output peak files into 1 list (called myfiles) - useful if analysing more than one motif
setwd("/Users/yuzhihuang/TEAD1_predicted_sites")

homer_out <- list.files(pattern = "*.txt") %>%
  lapply(read.delim) %>%
  bind_rows()

cat("Total HOMER motif hits:", nrow(homer_out), "\n")
cat("Unique peaks with interested motif:", length(unique(homer_out$PositionID)), "\n")


## ── 3. Merge HOMER hits with p2g links by PositionID ─────────────────────────
## PositionID in HOMER output = row number of the input p2g file
## Unique_Peak_ID in p2g = row number assigned during dataframe construction
## These must match if HOMER was run on the SAME p2g.df.obs.sub.txt
links <- merge(
  p2g.df.obs.sub,
  homer_out,
  by.x = "Unique_Peak_ID",
  by.y = "PositionID"
)

cat("Peaks after merge:", nrow(links), "\n")

## ── 4. Clean and format ───────────────────────────────────────────────────────
links <- links[, c("Unique_Peak_ID", "geneName", "Correlation",
                   "Motif.Name", "peakType", "peakName")]

links$Motif.Name <- toupper(
  str_split_fixed(links$Motif.Name, "\\(", n = Inf)[, 1]
)

links <- links %>%
  dplyr::filter(!duplicated(Unique_Peak_ID)) %>%
  dplyr::arrange(geneName)

## ── 5. Summary ────────────────────────────────────────────────────────────────
cat("HNF1B-linked peaks:", nrow(links), "\n") #16548
cat("Unique predicted HNF1B target genes:", length(unique(links$geneName)), "\n") #6992
table(links$peakType)

## ── 6. Save ───────────────────────────────────────────────────────────────────
write.csv(
  links,
  file      = "/Users/yuzhihuang/TEAD1_predicted_sites/Predicted_TEAD1_regulon.csv",
  row.names = FALSE
)

cat("Saved to: /Users/yuzhihuang/HNF1B_predicted_sites/Predicted_HNF1B_regulon.csv\n")

#Two types of overlaps are done here
#1: Overlap between Multiomics predicted sites of HNF1B & PAX8 CUT&RUN peaks
#2: Overlap between Multiomics predicted sites of HNF1B & PAX8 high confidence list
#I decided to use second overlap for greater stringency, tho first overlap produced interesting results, HNF1B GO showed tubular specific themes.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot
library(data.table)
projPath = "/Volumes/Elements/Leo_CUT_RUN"

##########################################################################
##########Find overlap between JP predicted sites and physical sites######
##########################################################################
sampleList = c("D7_H3K4me", "D7_IgG", "D7_PAX8", 
               "D10_H3K4me", "D10_IgG", "D10_PAX8",
               "D12_H3K4me", "D12_IgG", "D12_PAX8",
               "D14_H3K4me", "D14_IgG", "D14_PAX8")
peaks_SEACR <- list()
Timepoints <- c("D10","D12","D14")
for (day in Timepoints) {
  peaks_SEACR[[day]] <- read.table(
    paste0(projPath, "/peakCalling/SEACR/", day, "_PAX8_seacr_control.relaxed.bed"),
    sep = "\t",
    header = F,
    stringsAsFactors = FALSE
  )
}

# Keep only the first 3 columns and name them
for (day in Timepoints) {
  peaks_SEACR[[day]] <- peaks_SEACR[[day]][, 1:3]
  colnames(peaks_SEACR[[day]]) <- c("chr", "start", "end")
}
# View the result
head(peaks_SEACR[[day]])
tail(peaks_SEACR[[day]])

selected_days <- c("D10", "D12", "D14")
peaks_combined <- rbindlist(
  lapply(selected_days, function(day) {
    dt <- as.data.table(peaks_SEACR[[day]])
    dt <- dt[, .(chr, start, end)]
    dt[, Timepoint := day]   # optional: keep origin
    return(dt)
  })
)

peaks_multi <- read.csv(paste0(projPath, "/Data/predicted_HNF1B_regulon.csv"))
colnames(peaks_multi)
# Keep only the specific columns you need
peaks_multi <- peaks_multi[, c("peakName","geneName", "peakType")]
peaks_multi <- peaks_multi %>%
  separate(peakName, into = c("chr", "pos"), sep = ":") %>%
  separate(pos, into = c("start", "end"), sep = "-")

peaks_multi <- peaks_multi %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  )

nrow(peaks_multi)
#16548 peaks in regulated sites with HNF1B motif
length(unique(peaks_multi$geneName))
#6992 unique genes predicted to be regulated by HNF1B



length(peaks_combined$start)
#Total 16045 physical binding peaks in stringent setting (1st seq)
#Total of 60910 physical binding peaks across 3 timepoints in relaxed (1st seq)
#Total 10819 physical binding peaks in stringent setting (combined 1+2 seq)
#Total 37923 physical binding peaks in stringent setting (combined 1+2 seq)
setDT(peaks_multi)
peaks_multi <- peaks_multi[, .(chr, start, end, geneName, peakType)]
setkey(peaks_multi, chr, start, end)

setDT(peaks_combined)
setkey(peaks_combined, chr, start, end)

# Find overlaps: Which JP peaks are validated by SEACR?
peaks_overlap <- foverlaps(
  peaks_multi,                      # Query: JP predicted peaks
  peaks_combined,               #   Subject: SEACR physical peaks
  nomatch = 0L                   # CRITICAL: Only keep overlaps (0L not NULL!)
)

# Keep only the JP peak info (remove SEACR columns)
peaks_overlap <- peaks_overlap[, .(
  chr,
  start = i.start,
  end = i.end,
  geneName,
  peakType,
  Timepoint
)]

# Remove duplicates
peaks_overlap <- unique(peaks_overlap, by = c("chr", "start", "end","geneName"))

nrow(peaks_combined)
nrow(peaks_multi)

nrow(peaks_overlap) 
#Number of relaxed PHYSICAL SITES overlap with HNF1B 2957

length(unique(peaks_overlap$geneName)) 
#Number of unique genes that have PAX8 physical binding
#& overlap with HNF1B predicted regulon is 1985


peaks_overlap[, .N, by = Timepoint]
#number of peaks per timepoint

#Relaxed 1st seq:
#D12  1295  looks like its affected by sequencing depth
#D10   311  higher sequencing depth - clearer peaks
#D14   377  but we can still run GO analysis on validated peaks

#stringent 1st seq:
#D12   782
#D14   142
#D10   87

#stringent 1+2 combine seq
#D12   850
#D10   106
#D14   10

#relaxed combine seq
#1:       D12  1485
#2:       D10   391
#3:       D14    51

peaks_overlap[, .N, by = peakType][, .(peakType, Percentage = (N/sum(N))*100)]
# % of peaks by peaktype
#Promoter   25.06304
#  Distal   19.46546
#Intronic   41.70449
#  Exonic   13.76702 #Most are intronic???

#1+2 combined seq
#Promoter   32.38905
#Intronic   43.90935
#Distal   11.42587
#Exonic   12.27573


#########Get a txt list for validated peaks#########
for (day in selected_days) {
  
  cat("Processing", day, "\n")
  
  ## Subset by timepoint
  genes <- peaks_overlap[Timepoint == day, unique(geneName)]
  
  ## Remove NA (important)
  genes <- genes[!is.na(genes)]
  
  ## Save to file
  write.table(
    genes,
    file = paste0(projPath, "/validated/HNF1B/", day, "_HNF1B_PAX8_val_genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


###All timepoints target txt
genes <- unique(peaks_overlap$geneName)
genes <- genes[!is.na(genes)]
## Save to file
write.table(
  genes,
  file = paste0(projPath, "/validated/HNF1B/All_HNF1B_PAX8_val_genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
)



#Sanity check
pax8_targets <- read.table(paste0(projPath, "/Validated/All_PAX8_val_genes.txt"))$V1
overlap     <- read.table(paste0(projPath, "/Validated/HNF1B/All_HNF1B_PAX8_val_genes.txt"))$V1

length(pax8_targets)                    # total PAX8 targets
length(overlap)                          # PAX8 ∩ HNF1B
length(intersect(pax8_targets, overlap)) # should equal length(overlap) if subset
length(setdiff(overlap, pax8_targets))   # should be 0 if overlap is true subset

# Fraction of PAX8 targets that are also HNF1B-linked
length(overlap) / length(pax8_targets) * 100


##################
#Code below only overlaps PAX8 high confidence regulon with HNF1B predicted regulon
#Produce a list that is PAX8 regulated, and potentially co-regulated by HNF1B
##################
# HNF1B predicted regulon genes (from peaks_multi already loaded)
HNF1B_genes <- unique(peaks_multi$geneName)
HNF1B_genes <- HNF1B_genes[!is.na(HNF1B_genes)]

# PAX8 high confidence regulon genes
PAX8_genes <- read.table(
  paste0(projPath, "/Validated/All_PAX8_val_genes.txt"),
  header = FALSE
)$V1

# Overlap
PAX8_HNF1B_coreg <- intersect(PAX8_genes, HNF1B_genes)
length(PAX8_HNF1B_coreg)
#1117 genes overlap

#Hypergeometric test: whats the probability of 1117 overlap purely by chance
N     <- 12392      # All peak to gene link unique genes
K     <- 1433       # PAX8 target genes
n     <- 6992       # HNF1B predicted genes
k     <- 1117       # overlap

phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Save
write.table(
  PAX8_HNF1B_coreg,
  file = paste0(projPath, "/Validated/HNF1B/PAX8_HNF1B_coregulated_genes.txt"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)


###################################################
############GENE ONTOLOGY ANALYSIS#################
###################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Set project path
projPath <- "/Volumes/Elements/Leo_CUT_RUN"

# Define timepoints
Timepoints <- c("All")

# Loop through each timepoint
for (day in Timepoints) {
  
  cat("Processing", day, "...\n")
  
  # Load gene list for this timepoint
  gene_file <- paste0(projPath, "/Validated/HNF1B/", "PAX8_HNF1B_coregulated_genes.txt")
  genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)$V1
  
  cat("  Loaded", length(genes), "genes\n")
  
  # Run GO Enrichment
  go_results_BP <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = 'SYMBOL',
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_results_MF <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = 'SYMBOL',
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_results_CC <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = 'SYMBOL',
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # Simplify results
  go_simple_BP <- simplify(go_results_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_simple_MF <- simplify(go_results_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_simple_CC <- simplify(go_results_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # Visualize - Biological Process
  pdf(paste0(projPath, "/Output/HNF1B_PAX8_coregulon/", day, "_HNF1B+PAX8_GOBP_Dotplot.pdf"), width = 16, height = 12)
  print(dotplot(go_simple_BP, showCategory = 30) + 
          ggtitle(paste(day, "HNF1B/PAX8 Target Genes: Biological Processes")))
  dev.off()
  
  # Visualize - Molecular Function
  pdf(paste0(projPath, "/Output/HNF1B_PAX8_coregulon/", day, "_HNF1B+PAX8_GOMF_Dotplot.pdf"), width = 16, height = 12)
  print(dotplot(go_simple_MF, showCategory = 30) + 
          ggtitle(paste(day, "HNF1B/PAX8 Target Genes: Molecular Functions")))
  dev.off()
  
  # Visualize - Cellular Component
  pdf(paste0(projPath, "/Output/HNF1B_PAX8_coregulon/", day, "_HNF1B+PAX8_GOCC_Dotplot.pdf"), width = 16, height = 8)
  print(dotplot(go_simple_CC, showCategory = 20) + 
          ggtitle(paste(day, "HNF1B/PAX8 Target Genes: Cellular Components")))
  dev.off()
  
  # Save results to CSV
  write.csv(as.data.frame(go_simple_BP), 
            paste0(projPath, "/Output/HNF1B_PAX8_coregulon/", day, "_HNF1B+PAX8_GOBP_Results.csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_simple_MF), 
            paste0(projPath, "/Output/HNF1B_PAX8_coregulon/", day, "_HNF1B+PAX8_GOMF_Results.csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_simple_CC), 
            paste0(projPath, "/Output/HNF1B_PAX8_coregulon/", day, "_HNF1B+AX8_GOCC_Results.csv"),
            row.names = FALSE)
  
  cat("  Completed", day, "\n\n")
}

cat("All timepoints processed!\n")


#EMAPPLOT AND CNETPLOT
gene_file <- paste0(projPath, "/Validated/HNF1B/PAX8_HNF1B_coregulated_genes.txt")
genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)$V1
#import using above (already have a list) or below (using existing variable in environment)


go_res <- enrichGO(
  gene = genes,              # your gene list
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.01,
  readable = TRUE
)

go_res_simple <- simplify(go_res, cutoff = 0.8)

library(enrichplot)
go_res_simple <- pairwise_termsim(go_res_simple)
emapplot(
  go_res_simple, 
  showCategory = 20,
  size_category = 1.5,
  node_label_size = 4,
  nCluster = 3,
)

cnetplot(
  go_res_simple, 
  layout = "circle",
  color_category = "#a50f15",
  size_category = 0.8,
  color = "category",
  node_label = "share",
  color_item = "#fee5d9",
  size_item = 0.5,
  showCategory = 10,
  hilight = "none",
  hilight_alpha = 0.3,
  curvature = 0.1,
  size_edge = 0.8
)


####GO analysis figure#####
# 1. Convert to data frames and add a label for each ontology
df_bp <- as.data.frame(go_simple_BP)
if(nrow(df_bp) > 0) df_bp$Ontology <- "BP"

df_mf <- as.data.frame(go_simple_MF)
if(nrow(df_mf) > 0) df_mf$Ontology <- "MF"

df_cc <- as.data.frame(go_simple_CC)
if(nrow(df_cc) > 0) df_cc$Ontology <- "CC"

# 2. Bind them together into one big table
go_full_validated <- rbind(df_bp, df_mf, df_cc)

# 3. Write to a single CSV file
write.csv(go_full_validated, 
          file = paste0(projPath, "/validated/PAX8_val_GO_Combined.csv"), 
          row.names = FALSE)


rm(df_bp)
rm(df_mf)
rm(df_cc)

