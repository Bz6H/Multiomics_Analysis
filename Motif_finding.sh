# In terminal:
/bin/bash

# Add Homer to PATH
export PATH=/Users/yuzhihuang/Homer/bin:$PATH

# ── Settings ────────────────────────────────────────────────────────────────
# All HOMER output stays INSIDE the project folder - nothing is written to the
# home directory or other outer folders.
PROJ=/Users/yuzhihuang/MultiomeAnalysis
TF=TEAD1                              # transcription factor to scan
MOTIF=                                # leave empty; path to the TF .motif file (set below)
P2G=$PROJ/p2g.df.obs.sub.txt

SITES_DIR=$PROJ/${TF}_predicted_sites # per-TF output dir, inside $PROJ

# Peak-to-gene-links file p2g.df.obs.sub needs to be in correct format (MS-DOS .txt).
# FIRST - open the .csv, delete the first column, save as MS-DOS .txt
# THEN  - fix line endings:
changeNewLine.pl "$P2G"

# Create output folder for HOMER results (inside the project)
mkdir -p "$SITES_DIR"

# Provide the motif file, written into $SITES_DIR. Either copy a local HOMER known motif:
#   cp /Users/yuzhihuang/homer/data/knownTFs/motifs/${TF}.motif "$SITES_DIR/${TF}.motif"
# or download from JASPAR and convert to HOMER format with jaspar_homer_conversion.sh:
#   curl -L -o "$SITES_DIR/${TF}.jaspar" "https://jaspar.elixir.no/api/v1/matrix/MA0078.2.jaspar"
#   chmod +x "$PROJ/jaspar_homer_conversion.sh"
#   "$PROJ/jaspar_homer_conversion.sh" "$SITES_DIR/${TF}.jaspar" "$SITES_DIR/${TF}.motif"
MOTIF=${MOTIF:-$SITES_DIR/${TF}.motif}
head "$MOTIF"
# Verify it looks correct

# Now search peak-to-gene-links for instances of the TF motif.
# HOMER writes its preparsed/ tmp files into the -o directory, so run from $SITES_DIR.
cd "$SITES_DIR"
findMotifsGenome.pl \
    "$P2G" \
    hg38 \
    "$SITES_DIR" \
    -find "$MOTIF" \
    -size given \
    -mask \
    > "$SITES_DIR/${TF}.txt"

echo "HOMER complete. Output: $SITES_DIR/${TF}.txt"

# This generates a table of peaks containing over-represented TF motifs, annotated
# by a corresponding PositionID.

# MOVE ANALYSIS TO R TO JOIN WITH GENE NAMES (Coregulon_combined.R, Part 2)
