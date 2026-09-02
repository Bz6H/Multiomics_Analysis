# Multiome & CUT&RUN Analysis Pipeline

Analysis of iPSC-derived kidney organoid multiome (scATAC-seq + scRNA-seq) and CUT&RUN data across Day 7/10/12/14 timepoints.

## Pipeline Order

Run scripts in this order:

### 1. CUT&RUN (bash + R)

| Step | Script | Description |
|------|--------|-------------|
| 1a | `cutrun_all.sh` | Full CUT&RUN pipeline: FastQC, Trim Galore, Bowtie2 alignment (hg38 + dm6 spike-in), SAM/BAM conversion, Picard dedup, BED/bedGraph conversion, spike-in normalisation, ENCODE blacklist filtering, SEACR peak calling (stringent + relaxed), BigWig generation, deepTools heatmaps, GSEA gene signal extraction |
| 1b | `cutrun_analysis.R` | R companion to cutrun_all.sh: alignment QC plots, spike-in summary, duplication rates, fragment size distributions, FRiP analysis, PAX8 peak overlap with predicted regulon sites, Gene Ontology (BP/MF/CC), emapplot/cnetplot |

### 2. Multiome (R)

| Step | Script | Description |
|------|--------|-------------|
| 2a | `Multiomic_script.R` | Core processing: Arrow file creation from 10x fragment files, ArchR project setup, QC filtering, LSI dimensionality reduction (ATAC + RNA), Harmony batch correction, Louvain clustering, UMAP embedding, MACS2 peak calling, marker gene identification, atlas DEG matching, pseudotime, per-timepoint UMAPs |
| 2b | `JP_analysis_clauderefined.R` | Downstream analysis on `ArchRSubset_big3`: combined dims, Harmony, clustering, marker gene expression, JASPAR + Vierstra motif annotation, chromVAR deviations, cell type labelling, pseudotime trajectory assignment (dist/pod/med), gene expression along pseudotime, PAX8 target correlation, trajectory heatmaps, Monocle3 |
| 2c | `Motif_finding.sh` | HOMER motif scanning on peak-to-gene links. Run between Part 1 and Part 2 of Coregulon_combined.R |
| 2d | `Coregulon_combined.R` | Co-regulon analysis for HNF1B/WT1/SOX17/MECOM: peak-to-gene links (Part 1), HOMER motif merge + predicted regulon + PAX8 co-regulon + hypergeometric test + GO enrichment (Part 2, per-TF loop), CUT&RUN physical overlap validation (Part 3) |

### Utility

| Script | Description |
|--------|-------------|
| `jaspar_homer_conversion.sh` | Converts JASPAR `.jaspar` format to HOMER `.motif` format (used by `Motif_finding.sh`) |

## Data Requirements

### CUT&RUN (`cutrun_all.sh` / `cutrun_analysis.R`)
- Paired-end FASTQ files for 12 samples: D7/D10/D12/D14 x H3K4me/IgG/PAX8
- Bowtie2 indexes: hg38 (human) and dm6 (Drosophila spike-in)
- SEACR peak caller (v1.3)
- ENCODE hg38 blacklist (downloaded by the script)

### Multiome (`Multiomic_script.R`)
- 10x Multiome data per timepoint (D10/D12/D14):
  - `*_atac_fragments.tsv.gz` (ATAC fragment file)
  - `*_filtered_feature_bc_matrix.h5` (RNA counts)

### Coregulon (`Coregulon_combined.R`)
- ArchR project with peak-to-gene links (`ArchRSubset_big3`)
- `p2g.df.obs.sub.txt` (peak-to-gene links table, exported from Part 1)
- HOMER installed and on PATH
- TF motif files in HOMER format (via `jaspar_homer_conversion.sh`)

## Dependencies

### R packages
- **ArchR** v1.0.3, presto, Harmony
- **Bioconductor**: JASPAR2020, TFBSTools, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38, chromVAR, DESeq2, GenomicRanges, clusterProfiler, org.Hs.eg.db, enrichplot
- **CRAN**: tidyverse, dplyr, ggplot2, cowplot, viridis, pheatmap, openxlsx, readxl, Hmisc, data.table, ggpubr, corrplot, scales
- **GitHub**: monocle3, BPCells

### Command-line tools
- Bowtie2, samtools, bedtools, Picard
- HOMER (findMotifsGenome.pl, changeNewLine.pl)
- SEACR v1.3
- deepTools (bamCoverage, computeMatrix, plotHeatmap)
- Trim Galore, FastQC

## Genome

All analyses use **hg38** (GRCh38).

## Project Structure

```
MultiomeAnalysis/
  3_Day_10/               # Batch 3 raw 10x data (D10)
  3_Day_12/               # Batch 3 raw 10x data (D12)
  3_Day_14/               # Batch 3 raw 10x data (D14)
  ArchRProject/           # Full ArchR project (all cells)
  ArchRSubset_big3/       # Main analysis subset
  JasparMotif/            # Manual JASPAR motif files (WT1, PAX2, SALL1, MAFB, TWIST2)
  p2g.df.obs.sub.txt     # Peak-to-gene links (HOMER input)
  Nephrogenesis_atlas_DEG.xlsx  # Reference atlas DEGs for cluster annotation
  Trajectory_dist.rds     # Saved trajectory data
  Trajectory_pod.rds
  Trajectory_med.rds
```
