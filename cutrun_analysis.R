## ==== CUT&RUN R Analysis ==== ##
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
sampleList = c("D7_H3K4me", "D7_IgG", "D7_PAX8", 
               "D10_H3K4me", "D10_IgG", "D10_PAX8",
               "D12_H3K4me", "D12_IgG", "D12_PAX8",
               "D14_H3K4me", "D14_IgG", "D14_PAX8")
# Collect the alignment results from the bowtie2 alignment
alignResult = NULL
for(hist in sampleList) {
  alignRes = read.table(
    paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2.txt"),
    header = FALSE, 
    fill = TRUE
  )
  
  alignRate = substr(
    alignRes$V1[6], 
    1,
    nchar(as.character(alignRes$V1[6])) - 1
  )
  
  histInfo = strsplit(hist, "_")[[1]]
  
  alignResult = data.frame(
    Timepoint = histInfo[1],
    Target = histInfo[2],
    SequencingDepth = as.numeric(as.character(alignRes$V1[1])),
    MappedFragNum_hg38 = as.numeric(as.character(alignRes$V1[4])) + 
      as.numeric(as.character(alignRes$V1[5])),
    AlignmentRate_hg38 = as.numeric(alignRate)
  ) %>% 
    rbind(alignResult, .)
}
rm(alignRes)
rm(hist)
rm(histInfo)
rm(alignRate)


# Collect spike-in alignment results for summary
spikeAlign = NULL

for(hist in sampleList) {
  spikeRes = read.table(
    paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"),
    header = FALSE, 
    fill = TRUE
  )
  
  alignRate = substr(
    spikeRes$V1[6], 
    1,
    nchar(as.character(spikeRes$V1[6])) - 1
  )
  
  histInfo = strsplit(hist, "_")[[1]]
  
  spikeAlign = data.frame(
    Timepoint = histInfo[1],
    Target = histInfo[2],
    SequencingDepth = spikeRes$V1[1] %>% 
      as.character() %>% 
      as.numeric(),
    MappedFragNum_spikeIn = spikeRes$V1[4] %>% 
      as.character() %>% 
      as.numeric() + 
      spikeRes$V1[5] %>% 
      as.character() %>% 
      as.numeric(),
    AlignmentRate_spikeIn = alignRate %>% 
      as.numeric()
  ) %>% 
    rbind(spikeAlign, .)
}
rm(spikeRes)
spikeAlign$Timepoint = factor(spikeAlign$Timepoint)

spikeAlign %>% 
  mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))

# Join the alignment summary
alignSummary = left_join(alignResult, spikeAlign, 
                         by = c("Timepoint", "Target", "SequencingDepth")) %>%
  mutate(
    AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"),
    AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%")
  )

  
  
rm(alignResult, spikeAlign, hist, histInfo, alignRate)

# Generate plots to visualize sequencing depth and alignment rate to hg38
# and drosophila spike-in

fig1A = alignSummary %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14"))) %>%
  ggplot(aes(x = Timepoint, y = SequencingDepth / 1000000,
             fill = Target,
             label = round(SequencingDepth / 1000000, 1))) +  
  geom_col(position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth (Millions)") +
  xlab("") +
  ggtitle("A. Sequencing Depth")

fig1B = alignSummary %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14"))) %>%
  ggplot(aes(x = Timepoint, y = MappedFragNum_hg38 / 1000000, fill = Target)) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments (Millions)") +
  xlab("") +
  ggtitle("B. Mapped Fragments")

fig1C = alignSummary %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14")),
    AlignmentRate_hg38_num = as.numeric(gsub("%", "", AlignmentRate_hg38))
  ) %>%
  ggplot(aes(x = Timepoint, y = AlignmentRate_hg38_num, fill = Target)) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylim(0, 100) +
  ylab("Hg38 Alignment Rate (%)") +
  xlab("") +
  ggtitle("C. Alignment Rate")

fig1D = alignSummary %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14")),
    AlignmentRate_spikeIn_num = as.numeric(gsub("%", "", AlignmentRate_spikeIn))
  ) %>%
  ggplot(aes(x = Timepoint, y = AlignmentRate_spikeIn_num, fill = Target)) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylim(0, 10) +
  theme_bw(base_size = 18) +
  ylab("Spike-in Alignment Rate (%)") +
  xlab("Timepoint") +
  ggtitle("D. Spike-in Alignment Rate")

# Arrange plots in a grid
ggarrange(
  fig1A, fig1B, fig1C, fig1D,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
) %>%
  annotate_figure(top = text_grob("Alignment Summary", size = 18, face = "bold"))

######Sample Duplication########
##=== R command to measure degree of fragment duplication ===##
## Summarize the duplication information from the picard summary outputs
dupResult = c()
dupResult = data.frame()

for(hist in sampleList) {
  dupRes = read.table(
    paste0(projPath, "/alignment/bam/picard_summary/", hist, "_picard.dupMark.txt"),
    header = TRUE, 
    fill = TRUE, 
    skip = 6
  )
  
  histInfo = strsplit(hist, "_")[[1]]
  
  newRow = data.frame(
    Timepoint = histInfo[1],
    Target = histInfo[2],
    MappedFragNum_hg38 = as.numeric(as.character(dupRes$READ_PAIRS_EXAMINED[1])),
    DuplicationRate = as.numeric(as.character(dupRes$PERCENT_DUPLICATION[1])) * 100,
    EstimatedLibrarySize = as.numeric(as.character(dupRes$ESTIMATED_LIBRARY_SIZE[1]))
  ) %>%
    mutate(UniqueFragNum = MappedFragNum_hg38 * (1 - DuplicationRate / 100))
  
  dupResult = rbind(dupResult, newRow)
}

dupResult$Timepoint = factor(dupResult$Timepoint, levels = c("D7", "D10", "D12", "D14"))

alignDupSummary = left_join(
  alignSummary, 
  dupResult, 
  by = c("Target", "MappedFragNum_hg38")
) %>%
  mutate(DuplicationRate = paste0(DuplicationRate, "%"))

alignDupSummary$Timepoint.y = NULL
colnames(alignDupSummary)[colnames(alignDupSummary) == "Timepoint.x"] = "Timepoint"

names(alignDupSummary)
alignDupSummary



###### Fragment Length Distribution ########

# Initialize dataframe
fragLen <- data.frame()

for(hist in sampleList){
  
  histInfo <- strsplit(hist, "_")[[1]]
  
  tmp <- read.table(
    paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"),
    header = FALSE
  ) %>%
    mutate(
      fragLen = as.numeric(V1),
      fragCount = as.numeric(V2),
      Weight = fragCount / sum(fragCount),
      Timepoint = histInfo[1],
      Target = histInfo[2],
      sampleInfo = hist
    )
  
  fragLen <- rbind(fragLen, tmp)
}

# Factor ordering
fragLen$sampleInfo <- factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Timepoint <- factor(fragLen$Timepoint, levels = c("D7", "D10", "D12", "D14"))

  
  ## Violin plot (fragment distribution)
  
  fig5A <- fragLen %>%
  ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Target)) +
  geom_violin(bw = 5, scale = "width") +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  theme_bw(base_size = 18) +
  ggpubr::rotate_x_text(angle = 30) +
  ylab("Fragment Length (bp)") +
  xlab("") +
  ggtitle("A. Fragment Size Distribution")

  ## Line plot (fragment profile)
  
  fig5B <- fragLen %>%
  ggplot(aes(x = fragLen, y = fragCount,
             color = Target,
             group = sampleInfo,
             linetype = Timepoint)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                      option = "magma") +
  theme_bw(base_size = 18) +
  xlab("Fragment Length (bp)") +
  ylab("Fragment Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("B. Fragment Length Profile")

  ## Combine plots
  
  ggarrange(fig5A, fig5B, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  
########################################################
## ===== FRiP (Fraction of Reads in Peaks) ===== ##
library(chromVAR)

bamDir <- paste0(projPath, "/alignment/bam")
sampleList = c("D7_H3K4me", "D7_PAX8", 
               "D10_H3K4me", "D10_PAX8",
               "D12_H3K4me", "D12_PAX8",
               "D14_H3K4me", "D14_PAX8")
## Build FRiP counts across all samples
inPeakData <- data.frame()  # initialise as data.frame, not c()

for (hist in sampleList) {
  
  histInfo <- strsplit(hist, "_")[[1]]
  
  ## Read SEACR peaks and convert to GRanges
  peakRes <- read.table(
    paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_control.stringent.bed"),
    header = FALSE,
    fill   = TRUE
  )
  
  peak.gr <- GRanges(
    seqnames = peakRes$V1,
    ranges   = IRanges(start = peakRes$V2, end = peakRes$V3),
    strand   = "*"
  )
  
  ## Count fragments overlapping peaks
  bamFile <- paste0(bamDir, "/", hist, ".filtered.clean.bam")
  
  fragment_counts <- getCounts(
    bamFile,
    peak.gr,
    paired = TRUE,
    by_rg  = FALSE,
    format = "bam"
  )
  
  inPeakData <- rbind(
    inPeakData,
    data.frame(
      Timepoint   = histInfo[1],
      Target      = histInfo[2],
      FragInPeakN = counts(fragment_counts)[, 1] %>% sum()
    )
  )
}

## Join with alignment summary and compute FRiP
fripSummary <- alignSummary %>%
  left_join(inPeakData, by = c("Timepoint", "Target")) %>%
  mutate(
    FRiP = round(FragInPeakN / MappedFragNum_hg38 * 100, 2)
  ) %>%
  dplyr::select(                          # ← force dplyr's select
    Timepoint,
    Target,
    SequencingDepth,
    MappedFragNum_hg38,
    AlignmentRate_hg38,
    FragInPeakN,
    FRiP
  )

fripSummary

## ===== FRiP Figures ===== ##
fripSummary
fig6A <- fripSummary %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14"))) %>%
  filter(Target != "IgG") %>%
  ggplot(aes(x = Timepoint, y = FragInPeakN, fill = Target,
             label = scales::comma(FragInPeakN))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Fragments in Peaks") + xlab("") +
  ggtitle("A. Fragments in Peaks")

fig6B <- rbindlist(
  lapply(names(peaks_SEACR), function(d)
    data.table(Timepoint = d,
               width = peaks_SEACR[[d]]$end - peaks_SEACR[[d]]$start)
  )
) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14"))) %>%
  ggplot(aes(x = Timepoint, y = width, fill = Timepoint)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9,
                     option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log10",
                     breaks = c(200, 500, 1000, 3000, 10000)) +
  theme_bw(base_size = 18) +
  ylab("Peak Width (bp)") + xlab("") +
  ggtitle("B. Peak Width Distribution")

fig6C <- fripSummary %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14"))) %>%
  filter(Target != "IgG") %>%
  ggplot(aes(x = Timepoint, y = FragInPeakN, fill = Target,
             label = scales::comma(FragInPeakN))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55,
                     option = "magma", alpha = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 18) +
  ylab("Fragments in Peaks") + xlab("") +
  ggtitle("C. Fragments in Peaks")

fig6D <- fripSummary %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D7", "D10", "D12", "D14"))) %>%
  filter(Target != "IgG") %>%
  ggplot(aes(x = Timepoint, y = FRiP, fill = Target,
             label = paste0(FRiP, "%"))) +
  geom_bar(position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.3, size = 3.5) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55,
                     option = "magma", alpha = 0.8) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks") + xlab("") +
  ggtitle("D. FRiP Score")

ggarrange(
  fig6A, fig6B, fig6C, fig6D,
  ncol          = 2,
  nrow          = 2,
  common.legend = TRUE,
  legend        = "bottom"
) %>%
  annotate_figure(top = text_grob("Peak QC Summary", size = 18, face = "bold"))
  
##########################################################################
##########Find overlap between JP predicted sites and physical sites######
##########################################################################
sampleList = c("D7_H3K4me", "D7_IgG", "D7_PAX8", 
               "D10_H3K4me", "D10_IgG", "D10_PAX8",
               "D12_H3K4me", "D12_IgG", "D12_PAX8",
               "D14_H3K4me", "D14_IgG", "D14_PAX8")
peaks_SEACR <- list()
Timepoints <- c("D7","D10","D12","D14")
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

peaks_JP <- read.csv(paste0(projPath, "/Data/PAX8_regulon.csv"))
colnames(peaks_JP)
# Keep only the specific columns you need
peaks_JP <- peaks_JP[, c("peakName","geneName", "peakType")]
peaks_JP <- peaks_JP %>%
  separate(peakName, into = c("chr", "pos"), sep = ":") %>%
  separate(pos, into = c("start", "end"), sep = "-")

peaks_JP <- peaks_JP %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  )

nrow(peaks_JP)
#8383 peaks in regulated sites with PAX8 motif
length(unique(peaks_JP$geneName))
#4640 unique genes predicted to be regulated by PAX8



length(peaks_combined$start)
#Total 16045 physical binding peaks in stringent setting (1st seq)
#Total of 60910 physical binding peaks across 3 timepoints in relaxed (1st seq)
#Total 10819 physical binding peaks in stringent setting (combined 1+2 seq)
#Total 37923 physical binding peaks in relaxed setting (combined 1+2 seq)
setDT(peaks_JP)
peaks_JP <- peaks_JP[, .(chr, start, end, geneName, peakType)]
setkey(peaks_JP, chr, start, end)

setDT(peaks_combined)
setkey(peaks_combined, chr, start, end)

# Find overlaps: Which JP peaks are validated by SEACR?
peaks_overlap <- foverlaps(
  peaks_JP,                      # Query: JP predicted peaks
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
nrow(peaks_JP)

nrow(peaks_overlap) 
#Number of stringent overlap 1011 (1st seq)
#Number of relaxed overlap peaks 1983 (1st seq)
#Number of stringent overlap 1059 (1+2 combined seq)
#Number of relaxed overlap 1927 (1+2 combined seq)
length(unique(peaks_overlap$geneName)) 
#Number of genes 819 (1st seq stringent)
#Number of genes 1500 (1st seq relaxed)
#Number of genes 857 (1+2 combined stringent)
#Number of genes 1433 (1+2 combined relaxed)

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
#D12  1485
#D10   391
#D14    51

peaks_overlap[, .N, by = peakType][, .(peakType, Percentage = (N/sum(N))*100)]
# % of peaks by peaktype
#peakType on nearest gene not necessarily on target gene
#Promoter   25.06304
#  Distal   19.46546
#Intronic   41.70449
#  Exonic   13.76702 #Most are intronic??? because they might be false positive

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
    file = paste0(projPath, "/validated/", day, "_PAX8_val_genes.txt"),
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
  file = paste0(projPath, "/validated/All_PAX8_val_relaxed_genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
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
  gene_file <- paste0(projPath, "/Validated/", day, "_PAX8_val_genes.txt")
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
  pdf(paste0(projPath, "/Output/", day, "_PAX8_GOBP_Dotplot.pdf"), width = 16, height = 12)
  print(dotplot(go_simple_BP, showCategory = 30) + 
          ggtitle(paste(day, "PAX8 Target Genes: Biological Processes")))
  dev.off()
  
  # Visualize - Molecular Function
  pdf(paste0(projPath, "/Output/", day, "_PAX8_GOMF_Dotplot.pdf"), width = 16, height = 12)
  print(dotplot(go_simple_MF, showCategory = 30) + 
          ggtitle(paste(day, "PAX8 Target Genes: Molecular Functions")))
  dev.off()
  
  # Visualize - Cellular Component
  pdf(paste0(projPath, "/Output/", day, "_PAX8_GOCC_Dotplot.pdf"), width = 16, height = 8)
  print(dotplot(go_simple_CC, showCategory = 20) + 
          ggtitle(paste(day, "PAX8 STarget Genes: Cellular Components")))
  dev.off()
  
  # Save results to CSV
  write.csv(as.data.frame(go_simple_BP), 
            paste0(projPath, "/Output/", day, "_PAX8_GOBP_Results.csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_simple_MF), 
            paste0(projPath, "/Output/", day, "_PAX8_GOMF_Results.csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_simple_CC), 
            paste0(projPath, "/Output/", day, "_PAX8_GOCC_Results.csv"),
            row.names = FALSE)
  
  cat("  Completed", day, "\n\n")
}

cat("All timepoints processed!\n")

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
          file = paste0(projPath, "/Output/PAX8_val_GO_Combined.csv"), 
          row.names = FALSE)


rm(df_bp)
rm(df_mf)
rm(df_cc)

##########Generate Visualisation graphs for GO analysis across time##############
library(ggplot2)
library(dplyr)
library(readr)

# 1. Load your files (Update the paths to where your files are saved)
d7  <- read_csv(paste0(projPath,"/Output/D7_PAX8_GOCC_Results.csv")) %>% mutate(Day = "D7")
d10 <- read_csv(paste0(projPath,"/Output/D10_PAX8_GOCC_Results.csv")) %>% mutate(Day = "D10")
d12 <- read_csv(paste0(projPath,"/Output/D12_PAX8_GOCC_Results.csv")) %>% mutate(Day = "D12")
d14 <- read_csv(paste0(projPath,"/Output/D14_PAX8_GOCC_Results.csv")) %>% mutate(Day = "D14")

# 2. Combine them
all_data <- bind_rows(d7, d10, d12, d14)

# 3. Process GeneRatio and Significance
# We convert "278/7180" into a decimal number (0.038)
all_data <- all_data %>%
  mutate(
    RatioNumeric = sapply(GeneRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "/")))
      return(parts[1] / parts[2])
    }),
    LogP = -log10(p.adjust),
    Day = factor(Day, levels = c("D7", "D10", "D12", "D14")) # Keeps them in order
  )

# 4. Filter for specific terms you want to show (e.g., top 5 per day)
# For this example, let's just take a few interesting ones
target_terms <- c("focal adhesion",
                  "cell leading edge",
                  "cell cortex",
                  "adherence junction",
                  "ubiquitin ligase complex",
                  "cell-cell contact zone",
                  "lamellipodium",
                  "actin filament bundle",
                  "distal axon",
                  "synaptic membrame",
                  "neuron to neuron synapse",
                  "neuronal cell body",
                  "postsynaptic density")

plot_df <- all_data %>% filter(Description %in% target_terms)

#############PLOTTING across timepoints##############
plot_df$Description <- factor(plot_df$Description, levels = rev(target_terms))

ggplot(plot_df, aes(x = Day, y = Description)) +
  # Use GeneRatio for size and LogP for color
  geom_point(aes(size = RatioNumeric, color = LogP)) +
  
  # Set the color gradient to Red
  scale_color_gradient(low = "#fee5d9", high = "#a50f15", name = "-log10(p.adj)") +
  
  # Set the size legend name
  scale_size_continuous(name = "Gene Ratio") +
  
  # Aesthetic tweaks
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90")
  ) +
  labs(title = "PAX8 GOCC Enrichment",
       subtitle = "Size = Gene Ratio | Color = Significance",
       x = "Timepoint",
       y = "Cellular components")


######SNETPLOT AND EMAPPLOT#################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#Temporary code to genearte a input variable##
genes <- unique(peaks_overlap$geneName)
genes <- genes[!is.na(genes)]
####


go_res <- enrichGO(
  gene = genes,              # your gene list
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.01,
  readable = TRUE
)
go_res_simple <- simplify(go_res, cutoff = 0.7)

library(enrichplot)
go_res_simple <- pairwise_termsim(go_res_simple)
emapplot(
  go_res_simple, 
  showCategory = 15,
  size_category = 1.5,
  node_label_size = 5,
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
?emapplot
?cnetplot




sessionInfo()
