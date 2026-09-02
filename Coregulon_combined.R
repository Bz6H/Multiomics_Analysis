# =============================================================================
#  Coregulon_combined.R
#  Merge of HNF1B_coregulon.R + WT1_coregulon.R + SOX17_coregulon.R +
#  MECOM_coregulon.R into a single TF-parameterised script.
#
#  How it was combined:
#    * TFs to process are set once at the top (`TFs`). Everything TF-specific
#      runs inside `for (TF in TFs) { ... }`.
#    * Part 1 (peak-to-gene links) is identical in all four scripts, so it runs
#      once, before the loop.
#    * TF-specific object names (HNF1B_genes / WT1_genes / PAX8_SOX17_coreg ...)
#      are replaced by generic names (TF_genes, PAX8_TF_coreg) inside the loop.
#    * All paths are derived from `TF`; the "<TF>_predicted_sites" convention
#      from WT1/SOX17/MECOM is used (fixes HNF1B's stray TEAD1 paths).
#    * `motif_input` is added as requested and left empty.
#    * The CUT&RUN physical-overlap block existed only in HNF1B_coregulon.R;
#      kept as Part 3, parameterised by TF, toggled with `run_cutrun_validation`.
#    * emapplot/cnetplot trailing comma removed (was a parse error); those
#      plots are now written to PDF like the dotplots so the loop keeps them.
#    * The hypergeometric test uses counts computed from the data instead of
#      HNF1B's hard-coded N/K/n/k.
# =============================================================================

library(ArchR)
library(dplyr)
library(tidyverse)      # tidyr::separate, stringr::str_split_fixed used below
library(JASPAR2020)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)     # Part 3: foverlaps / rbindlist / setDT
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
## HNF1B_coregulon.R also load()-ed these in its CUT&RUN section; uncomment
## if a downstream edit needs them:
# library(viridis); library(GenomicRanges); library(chromVAR)
# library(DESeq2); library(ggpubr); library(corrplot)

set.seed(1)
addArchRThreads(threads = 4)
addArchRGenome("hg38")

# =============================================================================
#  USER SETTINGS
# =============================================================================
TFs <- c("HNF1B", "WT1", "SOX17", "MECOM")   # transcription factor(s) to process

## Motif input for the HOMER motif scan (run separately in bash - see
## Motif_finding.sh). Left empty on purpose: set this to your motif file before
## running HOMER. For several TFs use a named vector, e.g.
##   motif_input <- c(HNF1B = "", WT1 = "", SOX17 = "", MECOM = "")
motif_input <- ""

project_dir    <- "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject"
peakset_rds    <- "/Users/yuzhihuang/MultiomeAnalysis/ArchRProject/peakSet.rds"
p2g_txt        <- "/Users/yuzhihuang/MultiomeAnalysis/p2g.df.obs.sub.txt"
sites_base     <- "/Users/yuzhihuang/MultiomeAnalysis"  # parent of <TF>_predicted_sites/ (kept inside the project)
projPath       <- "/Volumes/Elements/Leo_CUT_RUN"     # CUT&RUN validation data
pax8_val_genes <- paste0(projPath, "/Validated/All_PAX8_val_genes.txt")

run_peak2gene         <- TRUE                 # Part 1 (slow); FALSE to reuse p2g_txt
run_cutrun_validation <- TRUE                 # Part 3 (was HNF1B-only)
seacr_days            <- c("D10", "D12", "D14")

# =============================================================================
#  PART 1 - Peak-to-gene links   (TF-independent, run once)
# =============================================================================
if (run_peak2gene) {

  ## ── 1. Load project ─────────────────────────────────────────────────────────
  proj <- loadArchRProject(project_dir)
  proj@peakSet <- readRDS(peakset_rds)
  is.null(proj@peakSet)   # should be FALSE
  length(proj@peakSet)    # should match parent

  ## ── 2. Compute peak-to-gene links ──────────────────────────────────────────
  proj <- addPeak2GeneLinks(
    ArchRProj        = proj,
    reducedDims      = "Harmony_ATAC",
    useMatrix        = "GeneExpressionMatrix",
    cellsToUse       = NULL,
    addEmpiricalPval = TRUE,
    addPermutedPval  = TRUE,
    nperm            = 10,
    maxDist          = 500000
  )

  ## ── 3. Build annotated peak-to-gene dataframe ──────────────────────────────
  p2geneDF          <- metadata(proj@peakSet)$Peak2GeneLinks
  p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>%
                          {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]

  annot                   <- readRDS(metadata(p2geneDF)$seATAC)
  p2geneDF$peakType       <- annot@rowRanges$peakType[p2geneDF$idxATAC]
  p2geneDF$nearestGene    <- annot@rowRanges$nearestGene[p2geneDF$idxATAC]
  p2geneDF$GroupReplicate <- annot@rowRanges$GroupReplicate[p2geneDF$idxATAC]

  p2geneDF.peaks         <- as.data.frame(metadata(p2geneDF)[[1]])
  p2geneDF.peaks$idxATAC <- rownames(p2geneDF.peaks)
  p2geneDF.merged        <- merge(p2geneDF, p2geneDF.peaks, by = "idxATAC")

  p2g.df.obs                <- as.data.frame(p2geneDF.merged)
  p2g.df.obs                <- p2g.df.obs[complete.cases(p2g.df.obs), ]
  p2g.df.obs$Unique_Peak_ID <- rownames(p2g.df.obs)
  p2g.df.obs$strand         <- "+"

  names(p2g.df.obs)[names(p2g.df.obs) == "seqnames"] <- "chromosome"
  names(p2g.df.obs)[names(p2g.df.obs) == "start"]    <- "starting position"
  names(p2g.df.obs)[names(p2g.df.obs) == "end"]      <- "ending position"

  p2g.df.obs <- p2g.df.obs[, c(19, 14, 15, 16, 18, 1, 2, 3, 4,
                               5, 6, 7, 8, 9, 10, 11, 12, 13, 17)]
  p2g.df.obs$GroupReplicate <- str_split_fixed(p2g.df.obs$GroupReplicate, "\\.", n = Inf)[, 1]

  ## ── 4. Filter to significant links ────────────────────────────────────────
  p2g.df.obs.sub <- p2g.df.obs %>%
    dplyr::filter(FDR < 1e-4, Correlation > 0.4)
  cat("Total significant peak-gene links:", nrow(p2g.df.obs.sub), "\n")

  ## HOMER (bash) reads this table; PositionID in HOMER output = row number here.
  ## NOTE (Motif_finding.sh): before findMotifsGenome.pl -find, this file is
  ## re-saved with MS-DOS line endings and its first column removed
  ## (changeNewLine.pl).
  write.table(p2g.df.obs.sub, p2g_txt, sep = "\t", quote = FALSE, row.names = FALSE)
}

# =============================================================================
#  >>> RUN HOMER IN BASH NOW  (Motif_finding.sh)
#      For each TF in `TFs`, run findMotifsGenome.pl -find <motif_input for TF>
#      on p2g_txt and drop the output .txt into
#      file.path(sites_base, paste0(TF, "_predicted_sites")).
# =============================================================================

# =============================================================================
#  PART 2 - per-TF predicted regulon + PAX8 co-regulon + GO
# =============================================================================
p2g.df.obs.sub <- read.delim(p2g_txt)
cat("Total significant peak-gene links:", nrow(p2g.df.obs.sub), "\n")

for (TF in TFs) {

  cat("\n=====================  ", TF, "  =====================\n")

  sites_dir   <- file.path(sites_base, paste0(TF, "_predicted_sites"))
  regulon_csv <- file.path(sites_dir, paste0("Predicted_", TF, "_regulon.csv"))

  ## ── merge HOMER motif hits with p2g links by PositionID ───────────────────
  ## PositionID (HOMER) = row number of p2g_txt = Unique_Peak_ID.
  setwd(sites_dir)
  homer_out <- list.files(sites_dir, pattern = "*.txt", full.names = TRUE) %>%
    lapply(read.delim) %>%
    bind_rows()
  cat("Total HOMER motif hits:", nrow(homer_out), "\n")
  cat("Unique peaks with", TF, "motif:", length(unique(homer_out$PositionID)), "\n")

  links <- merge(
    p2g.df.obs.sub, homer_out,
    by.x = "Unique_Peak_ID", by.y = "PositionID"
  )
  cat("Peaks after merge:", nrow(links), "\n")

  ## ── clean and format ─────────────────────────────────────────────────────
  links <- links[, c("Unique_Peak_ID", "geneName", "Correlation",
                     "Motif.Name", "peakType", "peakName")]
  links$Motif.Name <- toupper(str_split_fixed(links$Motif.Name, "\\(", n = Inf)[, 1])
  links <- links %>%
    dplyr::filter(!duplicated(Unique_Peak_ID)) %>%
    dplyr::arrange(geneName)

  ## ── summary + save predicted regulon ─────────────────────────────────────
  cat(TF, "-linked peaks:", nrow(links), "\n")
  cat("Unique predicted", TF, "target genes:", length(unique(links$geneName)), "\n")
  print(table(links$peakType))

  dir.create(sites_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(links, file = regulon_csv, row.names = FALSE)
  cat("Saved to:", regulon_csv, "\n")

  ## ── overlap predicted regulon with PAX8 high-confidence regulon ──────────
  peaks_multi <- read.csv(regulon_csv)
  peaks_multi <- peaks_multi[, c("peakName", "geneName", "peakType")]
  peaks_multi <- peaks_multi %>%
    separate(peakName, into = c("chr", "pos"),   sep = ":") %>%
    separate(pos,      into = c("start", "end"), sep = "-") %>%
    mutate(start = as.integer(start),
           end   = as.integer(end))
  cat("Rows in peaks_multi:", nrow(peaks_multi), "\n")

  TF_genes <- unique(peaks_multi$geneName)
  TF_genes <- TF_genes[!is.na(TF_genes)]

  PAX8_genes <- read.table(pax8_val_genes, header = FALSE)$V1

  PAX8_TF_coreg <- intersect(PAX8_genes, TF_genes)
  cat("PAX8 &", TF, "co-regulated genes:", length(PAX8_TF_coreg), "\n")

  ## hypergeometric test for the overlap (counts derived from the data)
  N <- length(unique(p2g.df.obs.sub$geneName))   # background: all p2g-linked genes
  K <- length(unique(PAX8_genes))                # PAX8 target genes
  n <- length(TF_genes)                          # predicted TF regulon genes
  k <- length(PAX8_TF_coreg)                     # overlap
  cat("Hypergeometric p =", phyper(k - 1, K, N - K, n, lower.tail = FALSE), "\n")

  coreg_dir <- paste0(projPath, "/Validated/", TF)
  dir.create(coreg_dir, recursive = TRUE, showWarnings = FALSE)
  coreg_txt <- paste0(coreg_dir, "/PAX8_", TF, "_coregulated_genes.txt")
  write.table(PAX8_TF_coreg, file = coreg_txt,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## ── GENE ONTOLOGY on the co-regulated gene list ─────────────────────────
  out_dir <- paste0(projPath, "/Output/", TF, "_PAX8_coregulon")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  Timepoints <- c("All")
  for (day in Timepoints) {

    cat("Processing", day, "...\n")
    gene_file <- coreg_txt
    genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)$V1
    cat("  Loaded", length(genes), "genes\n")

    go_results_BP <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                              ont = "BP", pAdjustMethod = "BH",
                              pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
    go_results_MF <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                              ont = "MF", pAdjustMethod = "BH",
                              pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
    go_results_CC <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                              ont = "CC", pAdjustMethod = "BH",
                              pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

    go_simple_BP <- simplify(go_results_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
    go_simple_MF <- simplify(go_results_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
    go_simple_CC <- simplify(go_results_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)

    pdf(paste0(out_dir, "/", day, "_", TF, "+PAX8_GOBP_Dotplot.pdf"), width = 16, height = 12)
    print(dotplot(go_simple_BP, showCategory = 30) +
            ggtitle(paste(day, paste0(TF, "/PAX8 Target Genes: Biological Processes"))))
    dev.off()

    pdf(paste0(out_dir, "/", day, "_", TF, "+PAX8_GOMF_Dotplot.pdf"), width = 16, height = 12)
    print(dotplot(go_simple_MF, showCategory = 30) +
            ggtitle(paste(day, paste0(TF, "/PAX8 Target Genes: Molecular Functions"))))
    dev.off()

    pdf(paste0(out_dir, "/", day, "_", TF, "+PAX8_GOCC_Dotplot.pdf"), width = 16, height = 8)
    print(dotplot(go_simple_CC, showCategory = 20) +
            ggtitle(paste(day, paste0(TF, "/PAX8 Target Genes: Cellular Components"))))
    dev.off()

    write.csv(as.data.frame(go_simple_BP),
              paste0(out_dir, "/", day, "_", TF, "+PAX8_GOBP_Results.csv"), row.names = FALSE)
    write.csv(as.data.frame(go_simple_MF),
              paste0(out_dir, "/", day, "_", TF, "+PAX8_GOMF_Results.csv"), row.names = FALSE)
    write.csv(as.data.frame(go_simple_CC),
              paste0(out_dir, "/", day, "_", TF, "+PAX8_GOCC_Results.csv"), row.names = FALSE)

    cat("  Completed", day, "\n\n")
  }
  cat("All timepoints processed!\n")

  ## ── EMAPPLOT / CNETPLOT on the co-regulated gene list ──────────────────
  gene_file <- coreg_txt
  genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)$V1

  go_res        <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                            ont = "BP", pvalueCutoff = 0.01, readable = TRUE)
  go_res_simple <- simplify(go_res, cutoff = 0.7)
  go_res_simple <- pairwise_termsim(go_res_simple)

  pdf(paste0(out_dir, "/All_", TF, "+PAX8_emapplot.pdf"), width = 12, height = 10)
  print(emapplot(go_res_simple,
                 showCategory     = 10,
                 size_category    = 1.5,
                 node_label_size  = 4,
                 nCluster         = 3))
  dev.off()

  pdf(paste0(out_dir, "/All_", TF, "+PAX8_cnetplot.pdf"), width = 12, height = 10)
  print(cnetplot(go_res_simple,
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
                 size_edge      = 0.8))
  dev.off()
}

# =============================================================================
#  PART 3 - CUT&RUN physical-site overlap   (was HNF1B_coregulon.R only)
#  Overlap predicted TF regulon peaks with PAX8 SEACR peaks, per timepoint.
#  The paper used the sequence-overlap result from Part 2 for stringency; this
#  is the alternative "physical overlap" approach.
# =============================================================================
if (run_cutrun_validation) {

  ## PAX8 SEACR peaks, all timepoints combined
  peaks_SEACR <- list()
  for (day in seacr_days) {
    peaks_SEACR[[day]] <- read.table(
      paste0(projPath, "/peakCalling/SEACR/", day, "_PAX8_seacr_control.relaxed.bed"),
      sep = "\t", header = FALSE, stringsAsFactors = FALSE
    )[, 1:3]
    colnames(peaks_SEACR[[day]]) <- c("chr", "start", "end")
  }

  peaks_combined <- rbindlist(lapply(seacr_days, function(day) {
    dt <- as.data.table(peaks_SEACR[[day]])[, .(chr, start, end)]
    dt[, Timepoint := day]
    dt
  }))
  setkey(peaks_combined, chr, start, end)

  for (TF in TFs) {

    cat("\n---------------  CUT&RUN overlap:", TF, "  ---------------\n")

    sites_dir   <- file.path(sites_base, paste0(TF, "_predicted_sites"))
    regulon_csv <- file.path(sites_dir, paste0("Predicted_", TF, "_regulon.csv"))

    peaks_multi <- read.csv(regulon_csv)[, c("peakName", "geneName", "peakType")]
    peaks_multi <- peaks_multi %>%
      separate(peakName, into = c("chr", "pos"),   sep = ":") %>%
      separate(pos,      into = c("start", "end"), sep = "-") %>%
      mutate(start = as.integer(start),
             end   = as.integer(end))

    nrow(peaks_multi)
    length(unique(peaks_multi$geneName))

    setDT(peaks_multi)
    peaks_multi <- peaks_multi[, .(chr, start, end, geneName, peakType)]
    setkey(peaks_multi, chr, start, end)

    ## Which predicted peaks are backed by a physical PAX8 peak?
    peaks_overlap <- foverlaps(
      peaks_multi,        # Query: predicted peaks
      peaks_combined,     # Subject: SEACR physical peaks
      nomatch = 0L        # keep overlaps only (0L, not NULL)
    )
    peaks_overlap <- peaks_overlap[, .(
      chr,
      start = i.start,
      end   = i.end,
      geneName,
      peakType,
      Timepoint
    )]
    peaks_overlap <- unique(peaks_overlap, by = c("chr", "start", "end", "geneName"))

    cat("Physical sites overlapping", TF, "predicted regulon:", nrow(peaks_overlap), "\n")
    cat("Unique genes (PAX8 physical binding &", TF, "regulon):",
        length(unique(peaks_overlap$geneName)), "\n")
    print(peaks_overlap[, .N, by = Timepoint])                                  # peaks per timepoint
    print(peaks_overlap[, .N, by = peakType][, .(peakType, Percentage = (N / sum(N)) * 100)])

    ## Per-timepoint and combined validated-gene lists
    val_dir <- paste0(projPath, "/validated/", TF)
    dir.create(val_dir, recursive = TRUE, showWarnings = FALSE)

    for (day in seacr_days) {
      cat("Processing", day, "\n")
      genes <- peaks_overlap[Timepoint == day, unique(geneName)]
      genes <- genes[!is.na(genes)]
      write.table(
        genes,
        file = paste0(val_dir, "/", day, "_", TF, "_PAX8_val_genes.txt"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )
    }

    genes <- unique(peaks_overlap$geneName)
    genes <- genes[!is.na(genes)]
    write.table(
      genes,
      file = paste0(val_dir, "/All_", TF, "_PAX8_val_genes.txt"),
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    ## Sanity check vs the PAX8 high-confidence gene list
    pax8_targets  <- read.table(pax8_val_genes)$V1
    overlap_genes <- read.table(paste0(val_dir, "/All_", TF, "_PAX8_val_genes.txt"))$V1
    cat("  PAX8 targets:", length(pax8_targets),
        "| physical overlap:", length(overlap_genes),
        "| not in PAX8 list:", length(setdiff(overlap_genes, pax8_targets)), "\n")
  }
}

sessionInfo()
