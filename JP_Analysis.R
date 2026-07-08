library(ArchR) #Version 1.0.3
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
plan('multisession', workers = 4)    # Parallelise where possible

#Check current working directory and files
getwd()

JPproj@peakSet


######SAVING, LOADING, CHECKING AVAILABLE EMBEDDINGS AND CLUSTERINGS##############
JPproj <- loadArchRProject(path = "/Users/yuzhihuang/MultiomeAnalysis/ArchRSubset_big3") 
saveArchRProject(ArchRProj = JPproj, outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRSubset_big3", load = FALSE)
names(JPproj@reducedDims)
#"LSI_ATAC"         "LSI_RNA"          "LSI_Combined"     "Harmony_ATAC"     "Harmony_Combined"
names(JPproj@embeddings)
#"UMAP_ATAC" "UMAP_RNA" <- these two dont know if Harmony'ed or not 
#"UMAP_Combined"<-(Made from LSI_combined - Harmony - AddUMAP)      "UMAP_Harmony_ATAC" <- JP paper used
names(JPproj@cellColData)
#Leo made Clusters_broad and Clusters_test_Combined (res0.4)

####################################################################################
####################################################################################
####################################################################################
JPproj <- addCombinedDims(
                        JPproj, 
                        reducedDims = c("LSI_ATAC", "LSI_RNA"), 
                        name =  "LSI_Combined",
                        dimWeights = c(5,5),
                        )
#Quality control of LSI###
rd <- getReducedDims(JPproj,"LSI_Combined")

dim(rd)
# cells × dimensions

any(duplicated(colnames(rd)))
# should be FALSE
rm(rd)
##########################

JPproj <- addHarmony(
  ArchRProj   = JPproj,
  reducedDims = 'LSI_RNA',  #editable, RNA/ATAC/Combined
  name        = 'Harmony_RNA',
  groupBy     = 'Sample',    # Corrects for sample-level batch effects
  theta       = c(2),        # Diversity penalty - increase if samples over-mix
  lambda      = c(1),        # Ridge regression - increase for stronger correction
  do_pca      = TRUE,
  scaleDims   = TRUE,
  corCutOff   = 0.75,
  force       = TRUE)

JPproj <- addClusters(
  input       = JPproj,
  reducedDims = "Harmony_Combined",       #Editable ATAC/RNA/Combined
  name        = "Clusters_test_Combined_prelabel",
  resolution  = 0.4,
  force       = TRUE)
table(JPproj$Clusters_test_Combined_prelabel) 

JPproj <- addUMAP(
  ArchRProj   = JPproj,
  reducedDims = 'Harmony_RNA',
  name        = 'UMAP_RNA',
  minDist     = 0.8,
  dims        = c(1:30),    # Try 1:20 if UMAP looks noisy
  force       = TRUE)

#Plot by cluster to check, change name to plot different cluster resolution
plotEmbedding(JPproj,
              name      = 'Clusters_test_Combined',  #Produced from addClusters
              embedding = 'UMAP_Harmony_ATAC',            #Produced from addUMAP
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 0.5,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')

# --- Also plot coloured by sample to check batch correction ---
plotEmbedding(JPproj,
              name      = 'Sample',
              embedding = 'UMAP_Combined',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 0.5)

#Impute weight
JPproj <- addImputeWeights(JPproj, reducedDims = 'Harmony_Combined')

####################################################################################
#######################Test marker gene expression dotplot##########################
library(ggplot2)
library(dplyr)
library(tidyr)

clust_col <- "Clusters_test_Combined"   # <-- change to the higher-res clustering column

# Confirm the codes you're about to filter on actually exist in that column
table(getCellColData(JPproj, select = clust_col, drop = TRUE))   # expect C1 … C10

clusters_test <- c("C1", 'C2',"C3",'C4','C5','C6','C7','C8','C9','C10')
genes <- c("MAFB","PTPRO","PODXL","WT1","CDH6",   # podocyte / precursor
           "HNF4A","LRP2","SLC3A1","GATA3","SPP1") # proximal vs distal
#  NPHS1, WT1 = podocyte (+) | HNF4A = proximal (+) | PAX8 = shared epithelial | TWIST1 = stromal (–) control

# ── Expression matrices ─────────────────────────────────────────────────────
seRNA    <- getMatrixFromProject(JPproj, useMatrix = "GeneExpressionMatrix")
expr_raw <- assay(seRNA)
expr_imp <- imputeMatrix(assay(seRNA), getImputeWeights(JPproj))  # needs addImputeWeights() to have been run
rownames(expr_raw) <- rowData(seRNA)$name
rownames(expr_imp) <- rowData(seRNA)$name   # names on BOTH matrices, not just one

genes       <- intersect(genes, rownames(expr_raw))
cluster_vec <- getCellColData(JPproj, select = clust_col, drop = TRUE)

# ── Long-format helper, shared by both matrices ─────────────────────────────
to_long <- function(mat) {
  t(as.matrix(mat[genes, , drop = FALSE])) %>%
    as.data.frame() %>%
    dplyr::mutate(cluster = cluster_vec) %>%
    dplyr::filter(cluster %in% clusters_test) %>%
    tidyr::pivot_longer(-cluster, names_to = "gene", values_to = "value")
}

pct_df <- to_long(expr_raw) %>%                       # % expressed  <- RAW
  dplyr::group_by(cluster, gene) %>%
  dplyr::summarise(pct = mean(value > 0) * 100, .groups = "drop")

avg_df <- to_long(expr_imp) %>%                       # magnitude    <- IMPUTED
  dplyr::group_by(cluster, gene) %>%
  dplyr::summarise(avg = mean(value), .groups = "drop") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(avg_scaled = avg / max(avg)) %>%      # each gene scaled to its own max across the 2 clusters
  dplyr::ungroup()

df <- dplyr::left_join(pct_df, avg_df, by = c("cluster", "gene"))
df$gene    <- factor(df$gene,    levels = genes)
df$cluster <- factor(df$cluster, levels = clusters_test)

# ── Plot ────────────────────────────────────────────────────────────────────
ggplot(df, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct, colour = avg_scaled)) +
  scale_colour_gradient(low = "white", high = "blue", name = "Scaled\nexpression") +
  scale_size_continuous(range = c(1, 8), name = "% expressed") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y     = element_text(size = 10),
    legend.position = "right",
    plot.title      = element_text(hjust = 0.5)
  ) +
  labs(x = NULL, y = NULL,
       title = "D14 podocyte (C1) vs proximal (C3) markers")

print(df[, c("cluster", "gene", "avg", "pct")])
view(df)
####################################################################################
####################################################################################
####################################################################################
# --- Check marker gene expression on UMAP to guide annotation ---
install.packages("hexbin")
# Replace with markers appropriate for your biological system
marker_genes <- c('LEF1','PAX8','HNF1B')

plot_list <- plotEmbedding(
  ArchRProj  = JPproj,
  colorBy    = 'GeneExpressionMatrix',
  name       = marker_genes,
  embedding  = 'UMAP_Harmony_ATAC',    #Editable ATAC/RNA/Combined
  plotAs = "points",
  size = 0.5,
  rastr = F,
  continuousSet = "greenBlue",
  imputeWeights = getImputeWeights(JPproj))

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
####################################################################
####################ATAC annotation#################################
####################################################################
#Manual loading of motifs for Pax2, Sall1, Mafb, Wt1 and Twist2 
# Load manual JASPAR files as PFM (raw counts)
WT1_pfm    <- readJASPARMatrix("~/MultiomeAnalysis/JasparMotif/MA1627.1.jaspar", matrixClass = "PFM")
PAX2_pfm   <- readJASPARMatrix("~/MultiomeAnalysis/JasparMotif/MA0067.2.jaspar", matrixClass = "PFM")
SALL1_pfm  <- readJASPARMatrix("~/MultiomeAnalysis/JasparMotif/UN1278.1.jaspar", matrixClass = "PFM")
MAFB_pfm   <- readJASPARMatrix("~/MultiomeAnalysis/JasparMotif/MA0117.2.jaspar", matrixClass = "PFM")
TWIST2_pfm <- readJASPARMatrix("~/MultiomeAnalysis/JasparMotif/MA0633.2.jaspar", matrixClass = "PFM")

# Convert PFM to PWM (log-odds scores)
WT1_pwm    <- toPWM(WT1_pfm,    pseudocounts = 0.8)
PAX2_pwm   <- toPWM(PAX2_pfm,   pseudocounts = 0.8)
SALL1_pwm  <- toPWM(SALL1_pfm,  pseudocounts = 0.8)
MAFB_pwm   <- toPWM(MAFB_pfm,   pseudocounts = 0.8)
TWIST2_pwm <- toPWM(TWIST2_pfm, pseudocounts = 0.8)

# Fix names
names(WT1_pwm)    <- "WT1"
names(PAX2_pwm)   <- "PAX2"
names(SALL1_pwm)  <- "SALL1"
names(MAFB_pwm)   <- "MAFB"
names(TWIST2_pwm) <- "TWIST2"


# Build JASPAR2020 set
pwm_set1 <- getMatrixSet(x = JASPAR2020,
                         opts = list(all_versions = FALSE, species = 9606,
                                     collection = "CORE", matrixtype = "PWM"))
pwm_set2 <- getMatrixSet(x = JASPAR2020,
                         opts = list(all_versions = FALSE, species = 9606,
                                     collection = "UNVALIDATED", matrixtype = "PWM"))
pwm_set <- c(pwm_set1, pwm_set2)
names(pwm_set) <- name(pwm_set)

# Merge manual motifs (only ONCE)
pwm_set <- c(pwm_set, WT1_pwm, PAX2_pwm, SALL1_pwm, MAFB_pwm, TWIST2_pwm)
names(pwm_set) <- make.unique(names(pwm_set))

# Verify all 5 manual motifs added
c("WT1", "PAX2", "SALL1", "MAFB", "TWIST2") %in% names(pwm_set)
length(pwm_set)



JPproj <- addMotifAnnotations(ArchRProj = JPproj,
                            motifPWMs = pwm_set,
                            name      = "JASPAR",
                            force     = TRUE)

JPproj <- addBgdPeaks(JPproj, force = TRUE)


JPproj <- addDeviationsMatrix(ArchRProj      = JPproj,
                            peakAnnotation = "Motif",
                            matrixName     = "JASPAR",
                            force          = TRUE)
#This step takes >20 minutes

JPproj <- addMotifAnnotations(ArchRProj  = JPproj,
                            motifSet   = "vierstra",
                            name       = "Vierstra",
                            collection = "archetype",
                            cutOff     = 5e-05,
                            width      = 7,
                            version    = 2,
                            force      = TRUE)


plotEmbedding(
  ArchRProj = JPproj,
  colorBy = "JASPAR",          # the deviations matrix
  name = paste0("deviations:", "SALL1"),              # specific TF motif name
  embedding = "UMAP_Harmony_ATAC",
  plotAs = "points",
  size = 0.5,
  rastr = F,
  continuousSet = "solarExtra",
  quantCut = c(0.05,0.95),
  imputeWeights = getImputeWeights(JPproj)
) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust = 0, vjust = 2.5),
        axis.title.y = element_text(hjust = 0, vjust = -1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        line = element_blank(),
        plot.title = element_text(hjust = 0, vjust = 0)) + 
  labs(x = 'UMAP 1', y = 'UMAP 2')




####################################################################################
###################Pseudotime gene expression pattern###############################
####################################################################################

# --- Assign cell type labels based on your marker gene analysis ---
# Modify the cluster-to-celltype mapping to match your data
cluster_map <- c(
  C1  = 'D14_Prox_Pod',
  C2  = 'D12_Prox',
  C3  = 'D14_Prox_Pod',
  C4  = 'D14_Distal',
  C5  = 'D10_Induced',
  C6  = 'D12_Distal',
  C7  = 'D14_Medial',
  C8  = 'D14_Stromal',
  C9  = 'D10_NPC',
  C10 = 'D12_Stromal'
)
new_labels <- cluster_map[as.character(JPproj$Clusters_test_Combined)]

JPproj <- addCellColData(
  ArchRProj = JPproj,
  data      = new_labels,
  name      = 'Clusters_test_Combined',
  cells     = getCellNames(JPproj),
  force     = TRUE
)
# --- Plot final annotated UMAP ---
plotEmbedding(JPproj,
              name      = 'Clusters_test_Combined',
              embedding = 'UMAP_Combined',
              colorBy   = 'cellColData',
              plotAs    = 'points',
              size      = 1,
              labelAsFactors = FALSE,
              labelMeans     = FALSE,
              discreteSet    = 'stallion')


#JP then assigned the trajectory backbone using the epithelial cluster he identified
#earlier, D10epi, D12 epi, D14epi.
Trajectory_dist <- c('D10_Induced','D12_Distal','D14_Distal') 
#In our case, D10 induced --> D12 induced --> D14 Proximal -- > D14 Distal
Trajectory_pod <-c('D10_Induced','D12_Prox','D14_Prox_Pod')
Trajectory_med <-c('D10_Induced','D12_Distal', 'D14_Medial')

JPproj <- addTrajectory( #3 edits
  ArchRProj = JPproj, 
  name = "Trajectory_dist",              #Edit to change trajectory
  groupBy = "Clusters_test_Combined",    #Edit to change, which clustering contain labels to identify items in trajectory
  trajectory = Trajectory_dist,          #Ordered vector of cluster names (specified earlier)
  embedding = 'UMAP_Harmony_ATAC',       #Which UMAP is used for spline fitting
  force = TRUE,
  seed = 1,
  useAll = F,
  spar = 2,
  preFilterQuantile = 1,
  postFilterQuantile = 1,
  saveDF = paste0("/Users/yuzhihuang/MultiomeAnalysis/Trajectory_dist.rds") #Editable
)

####Plot Trajectory
trajectory_cells <- getCellNames(
  JPproj[complete.cases(JPproj@cellColData$Trajectory_dist), ] #edit
)

JPproj@cellColData$Trajectory_plot <- JPproj@cellColData$Trajectory_dist #edit
JPproj@cellColData$Trajectory_plot[is.na(JPproj@cellColData$Trajectory_plot)] <- 0

traj_plot <- plotEmbedding(
  ArchRProj     = JPproj,
  colorBy       = "cellColData",
  name          = "Trajectory_plot",
  embedding     = "UMAP_Harmony_ATAC",
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

traj_plot
#Now with the trajectory I want to measure gene expression changes along it#
library(dplyr)

ann_colors.2 <- list(
  Cluster = c("D10_Induced" = "#9D2C24",
              "D12_Distal" = "#E25822",
              "D14_Distal" = "#F2C140")
)
#Assign colours to the clusters along the pseudotime track
label <- data.frame(
  cellName   = JPproj$cellNames,
  Pseudotime = JPproj$Trajectory_dist, #edit
  Cluster    = JPproj$Clusters_test_Combined
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
  ArchRProj = JPproj, 
  name = "Trajectory_dist",   #EDIT
  useMatrix = "GeneExpressionMatrix", 
  log2Norm = TRUE
) 
# Create a matrix of all gene's average expression along the pseudotime bin in the trajectory

####################################################################################
###########Quality check of marker#################
########by plotting a line graph along pseudo######
markers <- c('HNF1B','PAX8','LEF1')

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
  scale_colour_viridis_d(option = "B") + 
  labs(x = "Pseudotime bin", y = "log2 expression",
       title = "Marker gene expression along trajectory") +
  theme_classic(base_size = 12)

ggplot(marker_df, aes(x = bin, y = expression, colour = gene)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(HNF1B = "yellow", PAX8 = "red", LEF1 = "blue")) +
  labs(x = "Pseudotime bin", y = "log2 expression",
       title = "Marker gene expression along trajectory") +
  theme_classic(base_size = 12)


####################################################################################
#######Plot expression overlapping UMAP
####################################################################################

traj_exp <- plotTrajectory(
  ArchRProj = JPproj,
  trajectory = "Trajectory_dist",
  colorBy   = "GeneExpressionMatrix",
  name      = 'CDC42EP1',
  embedding = "UMAP_Harmony_ATAC",
  continuousSet = "fireworks2"
)


traj_exp[[1]]
traj_exp[[2]]
plotEmbedding(JPproj, name = "Clusters_test_Combined", embedding = "UMAP_Combined",
              colorBy = "cellColData", plotAs = "points", size = 1)
# ── SAVE 
saveArchRProject(ArchRProj = JPproj, outputDirectory = "/Users/yuzhihuang/MultiomeAnalysis/ArchRSubset_annotated", load = FALSE)
####################################################################################
####################################################################################
####################################################################################
library(dplyr)
library(openxlsx)

# --- Load validated gene list ---
val_genes <- readLines("/Volumes/Elements/Leo_CUT_RUN/Validated/All_PAX8_val_genes_stringent.txt")
val_genes <- unique(trimws(val_genes))
val_genes <- val_genes[val_genes != ""]

# --- Extract gene symbols from trajGEM rownames (format: "chr:GENE") ---
##If not ran before, run this
trajGEM <- getTrajectory(
  ArchRProj = JPproj, 
  name = "Trajectory_pod", 
  useMatrix = "GeneExpressionMatrix", 
  log2Norm = TRUE
) 
##
traj_mat   <- assay(trajGEM)                       # genes x 100 bins
gene_syms  <- sub(".*:", "", rownames(traj_mat))

# --- PAX8 reference vector --- (editable)
FOXC2_row <- which(gene_syms == "FOXC2")
stopifnot(length(FOXC2_row) == 1)
FOXC2_vec <- as.numeric(traj_mat[FOXC2_row, ])

# --- Subset matrix to validated genes present in trajGEM ---
keep      <- gene_syms %in% val_genes
sub_mat   <- traj_mat[keep, , drop = FALSE]
sub_syms  <- gene_syms[keep]

missing   <- setdiff(val_genes, sub_syms)
cat("Validated genes found in trajGEM:", length(sub_syms),
    "/ missing:", length(missing), "\n")

# --- Filter out flat/zero genes (correlation undefined) ---
row_var   <- apply(sub_mat, 1, var)
ok        <- row_var > 0
sub_mat   <- sub_mat[ok, , drop = FALSE]
sub_syms  <- sub_syms[ok]

# --- Correlation with PAX8 ---
pear <- apply(sub_mat, 1, function(x) cor(x, FOXC2_vec, method = "pearson"))
spear<- apply(sub_mat, 1, function(x) cor(x, FOXC2_vec, method = "spearman"))

# --- Simple classification ---
# co-peak  : strong positive  (Pearson >= 0.6)
# inverse  : strong negative  (Pearson <= -0.6)
# weak     : |r| < 0.3
# complex  : everything else (moderate / non-monotonic — Spearman vs Pearson disagree)
classify <- function(p, s) {
  if (is.na(p) || is.na(s)) return(NA_character_)
  if (p >=  0.6) return("co-peak")
  if (p <= -0.6) return("inverse")
  if (abs(p) < 0.3) return("weak")
  if (abs(p - s) > 0.25) return("complex")
  "moderate"
}
cls <- mapply(classify, pear, spear)

# --- Peak bin (where each gene is maximal) — useful for ordering ---   (editable)
peak_bin <- apply(sub_mat, 1, which.max)
min_bin <- apply(sub_mat, 1, which.min)
FOXC2_peak <- which.max(FOXC2_vec)

res <- data.frame(
  gene          = sub_syms,
  pearson_r     = round(pear, 4),
  spearman_r    = round(spear, 4),
  abs_pearson   = round(abs(pear), 4),
  peak_bin      = peak_bin,
  min_bin       = min_bin,
  bin_offset_vs_FOXC2 = peak_bin - FOXC2_peak,
  class         = cls,
  row.names     = NULL,
  stringsAsFactors = FALSE
)
res <- res[order(-res$abs_pearson), ]

# --- Write to Excel: ranked sheet + class-split sheets + missing list ---
install.packages("openxlsx")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "ranked_all")
writeData(wb, "ranked_all", res)

for (k in c("co-peak", "inverse", "complex", "moderate", "weak")) {
  sub <- res[!is.na(res$class) & res$class == k, ]
  addWorksheet(wb, k)
  writeData(wb, k, sub)
}

addWorksheet(wb, "missing_from_trajGEM")
writeData(wb, "missing_from_trajGEM", data.frame(gene = missing))

addWorksheet(wb, "FOXC2_trajectory") #editable
writeData(wb, "FOXC2_trajectory",
          data.frame(bin = 1:length(FOXC2_vec), FOXC2_log2 = FOXC2_vec))

out_path <- "/Users/yuzhihuang/MultiomeAnalysis/output170626/FOXC2_pseudotime_correlation_Trajectory_pod_stringent.xlsx" #filtered for targets of PAX8
saveWorkbook(wb, out_path, overwrite = TRUE)
cat("Written:", out_path, "\n")
####################################################################################
##############################Heatmaps##############################################
####################################################################################
trajMM  <- getTrajectory(JPproj, name = "Trajectory_pod", useMatrix = "JASPAR",                log2Norm = FALSE)  # your motif deviations
trajGEM <- getTrajectory(JPproj, name = "Trajectory_pod", useMatrix = "GeneExpressionMatrix",  log2Norm = TRUE)

plotTrajectoryHeatmap(trajGEM, pal = paletteContinuous("horizonExtra"))   # genes along prox-pod pseudotime
plotTrajectoryHeatmap(trajMM,  pal = paletteContinuous("solarExtra"))     # TF motifs

corr <- correlateTrajectories(trajMM, trajGEM)   # where does PAX8 motif track PAX8 expression?
corr[[1]]
####################################################################################
##########################Monocle 3.0 in ArchR#####################################
####################################################################################
# In R - point to homebrew hdf5
Sys.setenv(HDF5_DIR = "/opt/homebrew/opt/hdf5")
devtools::install_github("bnprks/BPCells/r")
devtools::install_github('cole-trapnell-lab/monocle3')
pseudoDist <- getMonocleTrajectories(
  ArchRProj = JPproj, 
  groupBy = "Clusters_test_Combined", 
  clusterParams = list(k = 100),
  useGroups = c("D10_Induced",
                "D12_Distal",
                "D14_Distal"),
  principalGroup = "D10_Induced", 
  embedding = "UMAP_Harmony_ATAC"
)

head(JPproj$pseudoDist[!is.na(JPproj$pseudoDist)])


