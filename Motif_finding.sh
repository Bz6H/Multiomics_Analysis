# In terminal:
/bin/bash

# Add Homer to PATH
export PATH=/Users/yuzhihuang/Homer/bin:$PATH

# Peak-to-gene-links file p2g.df.obs.sub needs to be in correct format (MS-DOS .txt.). 
# FIRST - OPEN .csv file, delete first column and save as MS-DOS .txt
# THEN - Run this in terminal to change format:

changeNewLine.pl /Users/yuzhihuang/MultiomeAnalysis/p2g.df.obs.sub.txt

#Create output folder for Homer results
mkdir -p /Users/yuzhihuang/SOX17_predicted_sites
    
# Copy the local HOMER transcription factor motif file
cp /Users/yuzhihuang/homer/data/knownTFs/motifs/mecom.motif \
   /Users/yuzhihuang/SOX17_predicted_sites/SOX17.motif

head /Users/yuzhihuang/SOX17_predicted_sites/SOX17.motif
#Or download from JASPAR and convert to HOMER format using jaspar_homer_conversion.sh script
cd /Users/yuzhihuang/SOX17_predicted_sites

curl -L -o SOX17_MA0078.jaspar "https://jaspar.elixir.no/api/v1/matrix/MA0078.2.jaspar"

head SOX17_MA0078.jaspar

chmod +x /Users/yuzhihuang/MultiomeAnalysis/jaspar_homer_conversion.sh
/Users/yuzhihuang/MultiomeAnalysis/jaspar_homer_conversion.sh \
  /Users/yuzhihuang/SOX17_predicted_sites/SOX17_MA0078.jaspar \
  /Users/yuzhihuang/SOX17_predicted_sites/SOX17.motif
head SOX17.motif
# Verify it looks correct

# Now search Peak-to-gene-links for instances of TF motifs (make sure TF.motif file is in output_folder)
cd /Users/yuzhihuang/SOX17_predicted_sites
findMotifsGenome.pl \
    /Users/yuzhihuang/MultiomeAnalysis/p2g.df.obs.sub.txt \
    hg38 \
    /Users/yuzhihuang/SOX17_predicted_sites \
    -find /Users/yuzhihuang/SOX17_predicted_sites/SOX17.motif \
    -size given \
    -mask \
    > /Users/yuzhihuang/SOX17_predicted_sites/SOX17.txt
 
echo "HOMER complete. Output: /Users/yuzhihuang/SOX17_predicted_sites/SOX17.txt"

# This generates a table of peaks containing overrepresentations of transcription factor motifs annonated by a corresponding PositionID

#MOVE ANALYSIS TO R TO JOIN WITH GENE NAMES