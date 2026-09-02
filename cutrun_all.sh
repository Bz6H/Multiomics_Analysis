#Downloading data from CRUK FTP server
# ftp ftp1.cruk.cam.ac.uk
# wget -r -np -nH --cut-dirs=1 -R "index.html*" \
# "ftp://<USERNAME>:<PASSWORD>@ftp1.cruk.cam.ac.uk/" \
# -P $projPath/data/

#Merging fastq files from 2 sequencing runs:
projPath="/Volumes/Elements/Leo_CUT_RUN"
samples=("D7_H3K4me" "D7_IgG" "D7_PAX8" "D10_H3K4me" "D10_IgG" "D10_PAX8" "D12_H3K4me" "D12_IgG" "D12_PAX8" "D14_H3K4me" "D14_IgG" "D14_PAX8")

for histName in "${samples[@]}"; do
    cat $projPath/data/2_${histName}_R1.fastq.gz > $projPath/data/${histName}_R1.fastq.gz
    echo "Merged R1 for ${histName}"
    cat $projPath/data/2_${histName}_R2.fastq.gz > $projPath/data/${histName}_R2.fastq.gz
    echo "Merged R2 for ${histName}"
done


!/bin/bash
# CUT&RUN Analysis Pipeline
# Date: 18 March 2026

# Set project paths
projPath="/Volumes/Elements/Leo_CUT_RUN"
cd "$projPath"
ls
cores=8

samples=("D14_H3K4me")

# Check if it worked by printing the array
echo ${samples[@]}
cores=8

# Run FastQC for all samples, need to download FASTQC and put it in a tool folder.
for histName in "${samples[@]}"; do
    mkdir -p ${projPath}/fastqFileQC/${histName}
    $projPath/tools/FastQC/fastqc -o ${projPath}/fastqFileQC/${histName} -f fastq /Volumes/Elements/Leo_CUT_RUN/data/${histName}_R1.fastq.gz
    $projPath/tools/FastQC/fastqc -o ${projPath}/fastqFileQC/${histName} -f fastq /Volumes/Elements/Leo_CUT_RUN/data/${histName}_R2.fastq.gz
done

#Adapter contamination found in all samples, so we will trim the reads before alignment.
#Do some trimming then realign.
# One-time setup: install trim_galore
# conda config --add channels conda-forge
# conda config --set channel_priority strict
# conda install -c bioconda trim-galore




cores=8
samples=("D7_H3K4me" "D7_IgG" "D7_PAX8" "D10_H3K4me" "D10_IgG" "D10_PAX8" "D12_H3K4me" "D12_IgG" "D12_PAX8" "D14_H3K4me" "D14_IgG" "D14_PAX8")
##perform alignment with trimmed fastq files
ref="$projPath/bowtie2Index/hg38"
mkdir -p $projPath/alignment/sam/bowtie2_summary
for histName in "${samples[@]}";
do
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
-I 10 -X 700 \
-p ${cores} \
-x ${ref} \
-1 ${projPath}/trimmed/${histName}_R1_val_1.fq.gz \
-2 ${projPath}/trimmed/${histName}_R2_val_2.fq.gz \
-S ${projPath}/alignment/sam/${histName}_bowtie2.sam \
&> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
done


#########Spike in alignment#####
#Create path to drosophila genome
spikeInRef="$projPath/bowtie2Index/dm6"

# Bowtie2 alignment to spike-in genome
for histName in "${samples[@]}";
do
echo "Processing Spike-in for: ${histName}..."
bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail \
--no-mixed --no-discordant --phred33 \
-I 10 -X 700 \
-p ${cores} \
-x ${spikeInRef} \
-1 ${projPath}/trimmed/${histName}_R1_val_1.fq.gz \
-2 ${projPath}/trimmed/${histName}_R2_val_2.fq.gz \
-S ${projPath}/alignment/sam/${histName}_bowtie2_spikeIn.sam \
&> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt
done


##Use samtools to look at mapped size distributions
mkdir -p $projPath/alignment/sam/fragmentLen
## Extract the 9th column from the alignment sam file which is the fragment length
for histName in "${samples[@]}"; do
    samtools view -F 0x04 \
        $projPath/alignment/sam/${histName}_bowtie2.sam | \
        awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
        sort | uniq -c | \
        awk -v OFS="\t" '{print $2, $1/2}' \
        > $projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt
done

#Convert to bam file for downstream processing
mkdir -p $projPath/alignment/bam
for histName in "${samples[@]}"; do
    samtools view -b -F 0x04 \
        $projPath/alignment/sam/${histName}_bowtie2.sam \
        > $projPath/alignment/bam/${histName}_bowtie2.bam
        echo "${histName} converted to BAM."
done
#Use picard to look at duplication rates, but not removing (yet) to preserve biological information
which picard
picard MarkDuplicates --version
mkdir -p $projPath/alignment/bam/picard_summary

for histName in "D14_H3K4me"; do
echo "Processing Picard MarkDuplicates for: ${histName}..."
    ## Sort by coordinate
    picard SortSam \
        I=$projPath/alignment/bam/${histName}_bowtie2.bam \
        O=$projPath/alignment/bam/${histName}_bowtie2.sorted.bam \
        SORT_ORDER=coordinate

    ## Add read group tags (required for MarkDuplicates)
    picard AddOrReplaceReadGroups \
        I=$projPath/alignment/bam/${histName}_bowtie2.sorted.bam \
        O=$projPath/alignment/bam/${histName}_bowtie2.sorted.RG.bam \
        RGID=${histName} \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=${histName}

    ## Mark duplicates
    picard MarkDuplicates \
        I=$projPath/alignment/bam/${histName}_bowtie2.sorted.RG.bam \
        O=$projPath/alignment/bam/${histName}_bowtie2.sorted.dupMarked.bam \
        METRICS_FILE=$projPath/alignment/bam/picard_summary/${histName}_picard.dupMark.txt

done

for histName in "${samples[@]}"; do
    # Index sorted BAM for IGV
    samtools index \
        $projPath/alignment/bam/${histName}_bowtie2.sorted.bam
    echo "${histName} indexed."
done

##############################################
##FILE CONVERSION: Convert BAM to BED##
##############################################
# Parameters
minQualityScore=30
cores=6
maxFragmentSize=1000

#=========================================
#BAM to BED Conversion Pipeline, with filtering for unmapped reads, duplicate and quality score higher than 30, fragment size less than 1000bp
#=========================================
mkdir -p $projPath/alignment/bed
for histName in "${samples[@]}"; do
    echo ""
    echo "Processing $histName..."
    echo "-----------------------------------------"
    
    inputBam="$projPath/alignment/bam/${histName}_bowtie2.sorted.dupMarked.bam"
    cleanedBam="$projPath/alignment/bam/${histName}_bowtie2_cleaned.bam"
    bedFile="$projPath/alignment/bed/${histName}_bowtie2.bed"
    fragmentFile="$projPath/alignment/bed/${histName}_bowtie2.fragments.bed"
    
    ## 1. Filter: Remove unmapped reads, duplicates, and low MAPQ reads
    echo "  Filtering reads (removing unmapped, duplicates, MAPQ < $minQualityScore)..."
    echo "    Total reads before filtering:"
    samtools view -c "$inputBam"
    
    samtools view -@ $cores -b -F 0x04 -F 0x400 -q $minQualityScore \
        "$inputBam" \
        -o "$cleanedBam"
    
    echo "    Total reads after filtering:"
    samtools view -c "$cleanedBam"
    
    ## 2. Convert BAM to BEDPE
    echo "  Converting BAM to BEDPE..."
    samtools sort -n -@ $cores "$cleanedBam" | \
        bedtools bamtobed -i stdin -bedpe \
        > "$bedFile"
    
    ## 3. Create Fragment File
    # Keep only: same chromosome pairs, fragment size < maxFragmentSize
    echo "  Creating fragment file (max size: $maxFragmentSize bp)..."
    awk -v max=$maxFragmentSize '$1==$4 && $6-$2 < max {print $1"\t"$2"\t"$6}' \
        "$bedFile" | \
        sort -k1,1 -k2,2n \
        > "$fragmentFile"
    
    fragmentCount=$(wc -l < "$fragmentFile")
    echo "    Fragments retained: $fragmentCount"
    echo "  ✓ Finished $histName"
done

##############This block generates fragment bam files for Bigwig##########
blacklist=$projPath/bowtie2Index/hg38-blacklist.v2.bed

for histName in "${samples[@]}"; do

    input_bam=$projPath/alignment/bam/${histName}_bowtie2_cleaned.bam
    output_bam=$projPath/alignment/bam/${histName}.filtered.clean.bam
    samtools view -h $input_bam | \
    awk '($1 ~ /^@/) || ($9 > -1000 && $9 < 1000)' | \
    samtools view -bS - | \
    bedtools intersect -v -abam stdin -b $blacklist \
    > $output_bam

    echo "Finished $histName"

done

###################################
#######Spike-in calibration########
#Scaling factor is defined as S, and use constant C = 10,000 to avoid small factors.
#S = C / (fragments_mapped_to_drosophila_genome)
#Normalized coverage = (primary_genome_coverage) * S

mkdir -p $projPath/alignment/bedgraph
chromSize=$projPath/bowtie2Index/hg38.chrom.sizes

for histName in "${samples[@]}"; do

    ## Get spike-in sequencing depth from spike-in SAM file
    seqDepth=$(samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l)
    echo "Spike-in sequencing depth for $histName is: $seqDepth"

    ## Normalize fragments using spike-in scale factor
    if [[ "$seqDepth" -gt "1" ]]; then
        mkdir -p $projPath/alignment/bedgraph

        scale_factor=$(echo "10000 / $seqDepth" | bc -l)
        echo "Spike-in scaling factor for $histName is: $scale_factor"

        bedtools genomecov \
            -bg \
            -scale $scale_factor \
            -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed \
            -g $chromSize \
            > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
    fi

done

#Get spike in factors for all samples

for histName in "${samples[@]}"; do
    seqDepth=$(samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l)
    if [[ "$seqDepth" -gt "1" ]]; then
        scale_factor=$(echo "10000 / $seqDepth" | bc -l)
        echo "${histName}: ${scale_factor}"
    else
        echo "${histName}: No spike-in reads mapped, cannot calculate scale factor."
    fi
done
#Spike-in sequencing depth for D7_H3K4me is:   146606
#Spike-in scaling factor for D7_H3K4me is: .06821003233155532515
#Spike-in sequencing depth for D7_IgG is:   770248
#Spike-in scaling factor for D7_IgG is: .01298283150361961342
#Spike-in sequencing depth for D7_PAX8 is:   150832
#Spike-in scaling factor for D7_PAX8 is: .06629892860931367349
#Spike-in sequencing depth for D10_H3K4me is:  1136500
#Spike-in scaling factor for D10_H3K4me is: .00879894412670479542
#Spike-in sequencing depth for D10_IgG is:   389360
#Spike-in scaling factor for D10_IgG is: .02568317238545305116
#Spike-in sequencing depth for D10_PAX8 is:   238496
#Spike-in scaling factor for D10_PAX8 is: .04192942439286193479
#Spike-in sequencing depth for D12_H3K4me is:   689200
#Spike-in scaling factor for D12_H3K4me is: .01450957632037144515
#Spike-in sequencing depth for D12_IgG is:   554978
#Spike-in scaling factor for D12_IgG is: .01801873227407212538
#Spike-in sequencing depth for D12_PAX8 is:  1773742
#Spike-in scaling factor for D12_PAX8 is: .00563779850733646719
#Spike-in sequencing depth for D14_H3K4me is:   399542
#Spike-in scaling factor for D14_H3K4me is: 0.02502865781319610954
#Spike-in sequencing depth for D14_IgG is:   518706
#Spike-in scaling factor for D14_IgG is: 0.01927874364283428377
#Spike-in sequencing depth for D14_PAX8 is:   170608
#Spike-in scaling factor for D14_PAX8 is: 0.05861389852761886898

###Filtering technical artifacts using the ENCODE blacklist#####
# Download hg38 blacklist
wget -P "$projPath/bowtie2Index" https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
gunzip $projPath/bowtie2Index/hg38-blacklist.v2.bed.gz

##Now include D12 samples for removing blacklist regions
samples=("D7_H3K4me" "D7_IgG" "D7_PAX8" "D10_H3K4me" "D10_IgG" "D10_PAX8"
    "D12_H3K4me" "D12_IgG" "D12_PAX8" "D14_H3K4me" "D14_IgG" "D14_PAX8")
# Show top 10 highest signal regions (excluding chrM) 
for histName in "${samples[@]}"; do
    echo -e "\nTop 10 peaks in ${histName}:"
    grep -v "chrM" $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph | \
        sort -k4,4nr | head -n 10
done

###before filtering blacklist regions (regions with high accessibility causing a peak in control + PAX8)
###we see chr3    93470000        93471000        chr3    91516200        93749200        High Signal Region
###having the highest signal across all samples, filter blacklist to remove such aritifacts

###############Filter blacklist from all samples###############
mkdir -p $projPath/alignment/bedgraph/filtered

for histName in "${samples[@]}"; do
    echo "Filtering $histName..."
    bedtools subtract \
        -a $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
        -b $projPath/bowtie2Index/hg38-blacklist.v2.bed \
        > $projPath/alignment/bedgraph/filtered/${histName}_bowtie2.fragments.normalized.filtered.bedgraph
done

# Now check D12_PAX8 top peaks after filtering
for histName in "${samples[@]}"; do
    echo -e "\nTop 10 peaks in ${histName} after blacklist filtering:"
    grep -v "chrM" $projPath/alignment/bedgraph/filtered/${histName}_bowtie2.fragments.normalized.filtered.bedgraph | \
        sort -k4,4nr | head -n 10
done

##########################################
############PEAK CALLING##################
##########################################
mkdir -p $projPath/peakCalling/SEACR
seacr="$projPath/SEACR/SEACR_1.3.sh"


########### This block Ensure consistent sorting behavior across different environments#######
export LC_ALL=C 
# 1. Create a "Mac-safe" copy of the SEACR script
# This removes any hidden non-ASCII characters that might be breaking it
tr -cd '\11\12\15\40-\176' < "$projPath/SEACR/SEACR_1.3.sh" > "$projPath/SEACR/SEACR_1.3_fixed.sh"
# 2. Set the encoding fix (Crucial for macOS)
export LC_ALL=C
# 3. Update your variable to use the new fixed script
seacr_fixed="$projPath/SEACR/SEACR_1.3_fixed.sh"
    chmod +x /Volumes/Elements/Leo_CUT_RUN/SEACR/SEACR_1.3_fixed.sh
##################################################################################################
##--- PEAK CALL Method 1: Using IgG as Control ---##
Timepoints=("D14")

# Peak calling for H3K4me using IgG control
for time in "${Timepoints[@]}"; do
    echo "Calling peaks for ${time}_H3K4me..."
    bash $seacr_fixed \
        $projPath/alignment/bedgraph/filtered/${time}_H3K4me_bowtie2.fragments.normalized.filtered.bedgraph \
        $projPath/alignment/bedgraph/filtered/${time}_IgG_bowtie2.fragments.normalized.filtered.bedgraph \
        non \
        stringent \
        $projPath/peakCalling/SEACR/${time}_H3K4me_seacr_control.stringent
done

# Peak calling for PAX8 using IgG control
for time in "${Timepoints[@]}"; do
    echo "Calling peaks for ${time}_PAX8..."
    bash $seacr_fixed \
        $projPath/alignment/bedgraph/filtered/${time}_PAX8_bowtie2.fragments.normalized.filtered.bedgraph \
        $projPath/alignment/bedgraph/filtered/${time}_IgG_bowtie2.fragments.normalized.filtered.bedgraph \
        non \
        stringent \
        $projPath/peakCalling/SEACR/${time}_PAX8_seacr_control.stringent
done

# Count the number of peaks found
for time in "${Timepoints[@]}"; do
    echo "Number of peaks for ${time}_H3K4me:"
    wc -l $projPath/peakCalling/SEACR/${time}_H3K4me_seacr_control.stringent.bed
    echo "Number of peaks for ${time}_PAX8:"
    wc -l $projPath/peakCalling/SEACR/${time}_PAX8_seacr_control.stringent.bed
done
#D7_H3K4me peaks:6648   13389
#D7_PAX8 peaks:15341    9718
#D10_H3K4me peaks:5679  12678
#D10_PAX8 peaks:17484   2853
#D12_H3K4me peaks:8066  9797
#D12_PAX8 peaks:3510    4684
#D14_H3K4me peaks:0     7227
#D14_PAX8 peaks:9264    3282
##############################################################Relaxed peak calling
Timepoints=("D7" "D10" "D12" "D14")

# Peak calling for H3K4me using IgG control
for time in "${Timepoints[@]}"; do
    echo "Calling peaks for ${time}_H3K4me..."
    bash $seacr_fixed \
        $projPath/alignment/bedgraph/filtered/${time}_H3K4me_bowtie2.fragments.normalized.filtered.bedgraph \
        $projPath/alignment/bedgraph/filtered/${time}_IgG_bowtie2.fragments.normalized.filtered.bedgraph \
        non \
        relaxed \
        $projPath/peakCalling/SEACR/${time}_H3K4me_seacr_control
done
# Peak calling for PAX8 using IgG control
for time in "${Timepoints[@]}"; do
    echo "Calling peaks for ${time}_PAX8..."
    bash $seacr_fixed \
        $projPath/alignment/bedgraph/filtered/${time}_PAX8_bowtie2.fragments.normalized.filtered.bedgraph \
        $projPath/alignment/bedgraph/filtered/${time}_IgG_bowtie2.fragments.normalized.filtered.bedgraph \
        non \
        relaxed \
        $projPath/peakCalling/SEACR/${time}_PAX8_seacr_control
done

# Count the number of peaks found
for time in "${Timepoints[@]}"; do
    echo "Number of peaks for ${time}_H3K4me:"
    wc -l $projPath/peakCalling/SEACR/${time}_H3K4me_seacr_control.relaxed.bed
    echo "Number of peaks for ${time}_PAX8:"
    wc -l $projPath/peakCalling/SEACR/${time}_PAX8_seacr_control.relaxed.bed
done
#####Stringent########
#D7_PAX8: 8276  9718
#D10_PAX8: 3545 2853
#D12_PAX8: 4884 4684
#D14_PAX8: 7616 3282


#####Relaxed########
#D7_H3K4me peaks:6648   15589
#D7_PAX8 peaks:15341    9890
#D10_H3K4me peaks:5679  14763
#D10_PAX8 peaks:17484   24815
#D12_H3K4me peaks:8066  11885
#D12_PAX8 peaks:3510    9685
#D14_H3K4me peaks:0     10319
#D14_PAX8 peaks:9264    3423


####################################################
#####GENERATING BIGWIG FILES FOR VISUALIZATION######
####################################################
mkdir -p $projPath/alignment/bigwig
# Create the deeptools environment
# One-time setup: create deeptools environment
# conda create -n deeptools -c bioconda -c conda-forge deeptools python=3.9 -y
# source /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh
conda activate deeptools
deeptools --version
cores=6

# 1. Sort and Index BAM files (MANDATORY for bamCoverage)
for histName in "${samples[@]}"; do
    echo "Sorting and indexing ${histName}..."
    samtools sort -@ 8 \
        -o $projPath/alignment/bam/${histName}.sorted.bam \
        $projPath/alignment/bam/${histName}.filtered.clean.bam
    samtools index $projPath/alignment/bam/${histName}.sorted.bam
done
echo ${samples[@]}
# 2. Generate Normalized BigWig files using spike-in calibration
for histName in "${samples[@]}"; do
    echo "Generating BigWig for ${histName}..."
    
    # Get spike-in depth and calculate scale factor (same as bedGraph)
    seqDepth=$(samtools view -c -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam)
    
    if [[ "$seqDepth" -gt "1" ]]; then
        scale_factor=$(echo "10000 / $seqDepth" | bc -l)
        echo "  Scale factor: ${scale_factor}"
        
        bamCoverage \
            -b $projPath/alignment/bam/${histName}.sorted.bam \
            -o $projPath/alignment/bigwig/${histName}_normalized.bw \
            --binSize 10 \
            --normalizeUsing None \
            --scaleFactor ${scale_factor} \
            -p 8
    fi
done


###GENERATING HEATMAPS FOR PEAK SIGNAL INTENSITY#####
conda activate deeptools
cores=6
#Create a gene bed to compute matrix for heatmap
# 1. Ensure directory exists
mkdir -p $projPath/data/hg38_gene/

# 2. Download RefSeq annotations
echo "Downloading hg38 RefSeq gene annotations..."
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz" | \
    gunzip -c > $projPath/data/hg38_gene/refGene_raw.txt

# 3. Convert to clean BED format in ONE stepx
echo "Converting to BED format..."
awk 'BEGIN{OFS="\t"} {
    # Extract: chr, start, end, gene_name, score, strand
    print $3, $5, $6, $13, "0", $4
}' $projPath/data/hg38_gene/refGene_raw.txt | \
    grep "^chr" | \
    grep -v "chrUn" | \
    grep -v "_random" | \
    awk '$2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3' \
    > $projPath/data/hg38_gene/hg38_genes_clean.bed

# 4. Verify the result
echo -e "\n=== Verification ==="
echo "Total genes: $(wc -l < $projPath/data/hg38_gene/hg38_genes_clean.bed)"
echo -e "\nFirst 5 genes:"
head -n 5 $projPath/data/hg38_gene/hg38_genes_clean.bed
echo -e "\nColumn check (should see numbers in columns 2 and 3):"
head -n 3 $projPath/data/hg38_gene/hg38_genes_clean.bed | column -t

###############Generate heatmaps for each timepoint##############
conda activate deeptools
days=("D14")

for day in "${days[@]}"; do
    echo "Processing $day..."
    
    # UNIQUE FILENAMES PER LOOP
    matrix_out="$projPath/peakCalling/SEACR/${day}_PAX8_H3K4me_genes.mat.gz"
    output_png="$projPath/peakCalling/SEACR/${day}_PAX8_H3K4me_Heatmap.png"
    gene_bed="$projPath/data/hg38_gene/hg38_genes_clean.bed"
    output_bed="$projPath/peakCalling/SEACR/${day}_PAX8_genes_sorted.bed"
    
    # 1. Compute Matrix for specific day
    computeMatrix scale-regions \
        -S "$projPath/alignment/bigwig/${day}_H3K4me_normalized.bw" \
           "$projPath/alignment/bigwig/${day}_PAX8_normalized.bw" \
           "$projPath/alignment/bigwig/${day}_IgG_normalized.bw" \
        -R "$gene_bed" \
        --beforeRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 3000 \
        --skipZeros \
        -o "$matrix_out" \
        -p "$cores"
    
    # 2. Plot Heatmap for specific day and save sorted regions to BED
    plotHeatmap \
        -m "$matrix_out" \
        -out "$output_png" \
        --colorMap 'Blues' 'Oranges' \
        --samplesLabel "H3K4me" "PAX8" "IgG" \
        --plotTitle "${day} Enrichment" \
        --sortUsingSamples 2 \
        --sortRegions descend \
        --sortUsing mean \
        --outFileSortedRegions "$output_bed"\
        --yMin 0 0 \
        --yMax 0.3 0.3
done


# Plot Heatmap for specific matrix
plotHeatmap \
    -m "$projPath/peakCalling/SEACR/D14_PAX8_H3K4me_genes.mat.gz" \
    -out "$projPath/peakCalling/SEACR/D14_PAX8_H3K4me_Heatmap.png" \
    --colorMap 'Blues' 'Oranges' \
    --samplesLabel "H3K4me" "PAX8" "IgG" \
    --plotTitle "${day} Enrichment" \
    --sortUsingSamples 2 \
    --sortRegions descend \
    --sortUsing mean \
    --outFileSortedRegions "$output_bed" \
    --colorMap 'Blues' 'Oranges' 'Greys'\
    --zMin 0 0 0\
    --zMax 0.3 0.3 0.3\
    --yMin 0 0 0\
    --yMax 0.3 0.3 0.3 



##################################################################################
##### PAX8 Peak Heatmaps to check CUT&RUN quality & get a target gene list########
###################################################################################

cores=8
Timepoints=("D7" "D10" "D12" "D14")

for day in "${Timepoints[@]}"; do
    echo "Processing PAX8 $day..."

    baseName="${day}_PAX8"

    peak_bed="$projPath/peakCalling/SEACR/${baseName}_seacr_control.stringent.bed"
    summit_bed="$projPath/peakCalling/SEACR/${baseName}.summitRegion.bed"
    bigwig="$projPath/alignment/bigwig/${baseName}_normalized.bw"
    matrix_out="$projPath/peakCalling/SEACR/${baseName}.mat.gz"
    heatmap_out="$projPath/peakCalling/SEACR/${baseName}_heatmap.png"

    ## 1. Get peak center (summit approximation)
    awk 'BEGIN{OFS="\t"} {
        center = int(($2 + $3) / 2)
        print $1, center, center+1
    }' $peak_bed > $summit_bed

    ## 2. Compute matrix (+/- 3kb)
    computeMatrix reference-point \
        -S $bigwig \
        -R $summit_bed \
        --referencePoint center \
        -b 3000 -a 3000 \
        --skipZeros \
        -o $matrix_out \
        -p $cores

    ## 3. Plot heatmap
    plotHeatmap \
        -m $matrix_out \
        -out $heatmap_out \
        --sortUsing sum \
        --regionsLabel "PAX8 Peaks" \
        --samplesLabel "PAX8 ${day}" \
        --colorMap Reds \
        --zMin 0 \
        --zMax 1 \
        --yMax 0.3 0.3

    echo "Completed $day"
done

echo "All PAX8 heatmaps done!"


#########################################
#####GSEA analysis for PAX8 peaks########
#########################################
#Generate a txt file with gene names and intensity, sorted.
gene_bed="$projPath/data/hg38_gene/hg38_genes_clean.bed"
#####################creating hg38 promoter region bed file for gene mapping#####################
# Sort the gene BED by chromosome (alphabetical) then start position (numeric)
sort -k1,1 -k2,2n "$gene_bed" > "$gene_bed_sorted"
gene_bed_sorted="$projPath/data/hg38_gene/hg38_genes_sorted.bed"
# 1. Define the TSS points (Start of gene)
# Column 2 is Start, Column 3 is End, Column 6 is Strand
awk 'BEGIN{OFS="\t"} {
    if($6=="+") {print $1, $2, $2+1, $4, $5, $6} 
    else {print $1, $3-1, $3, $4, $5, $6}
}' "$projPath/data/hg38_gene/hg38_genes_sorted.bed" > $projPath/data/hg38_gene/hg38_tss_points.bed

# 2. Expand 500 bp in both directions (Total 1kb window)
# Note: You need a chrom.sizes file so it doesn't create negative coordinates
# If you don't have one, you can skip -g but it might throw errors at chromosome edges
bedtools slop -i $projPath/data/hg38_gene/hg38_tss_points.bed -g "$projPath/Data/hg38_gene/hg38.chrom.sizes" -b 500 > $projPath/data/hg38_gene/hg38_promoters_1kb.bed

# 3. Final Sort (Crucial for bedtools map)
sort -k1,1 -k2,2n $projPath/data/hg38_gene/hg38_promoters_1kb.bed > "$projPath/data/hg38_gene/hg38_promoters_sorted.bed"
#########################################

promotor_bed_sorted="$projPath/data/hg38_gene/hg38_promoters_sorted.bed"
# Filter out the ALT contigs from your promoter file
grep -E "^chr([0-9]+|[XYM])[[:space:]]" "$projPath/data/hg38_gene/hg38_promoters_sorted.bed" > "$projPath/data/hg38_gene/hg38_promoters_final.bed"
promotor_bed_final="$projPath/data/hg38_gene/hg38_promoters_final.bed"

for day in D7 D10 D12 D14; do
    echo "Processing $day..."
    
    bedgraph="$projPath/alignment/bedgraph/filtered/${day}_PAX8_bowtie2.fragments.normalized.filtered.bedgraph"
    sorted_bedgraph="/tmp/${day}_PAX8_sorted.bedgraph"
    output="$projPath/peakCalling/SEACR/${day}_PAX8_gene_signal_ranked.txt"
    
    echo "  Filtering and sorting bedGraph..."
    # This keeps only chr1-22, X, Y, and M
    grep -E "^chr([0-9]+|[XYM])[[:space:]]" "$bedgraph" | \
    LC_ALL=C sort -k1,1 -k2,2n > "$sorted_bedgraph"
    
    # Step 2: Map normalized signal to genes
    echo "  Mapping signal to genes..."
    bedtools map \
        -a $promotor_bed_final \
        -b $sorted_bedgraph \
        -g $projPath/data/hg38_gene/hg38.chrom.sizes \
        -c 4 \
        -o mean | \
        awk 'BEGIN{OFS="\t"} {
            gene = $4
            signal = $7
            if (signal != "." && signal > 0) {
                print gene, signal
            }
        }' | \
        sort -k2,2nr \
        > $output
    
    # Clean up
    rm $sorted_bedgraph

    #Sum signals for repetitive genes
    awk '{sum[$1] += $2} END {for (gene in sum) print gene, sum[gene]}' $projPath/peakCalling/SEACR/${day}_PAX8_gene_signal_ranked.txt | sort -k2,2nr > $projPath/peakCalling/SEACR/${day}_PAX8_summed.txt
    
    echo "  Genes with PAX8 signal: $(wc -l < $projPath/peakCalling/SEACR/${day}_PAX8_summed.txt)"
    echo "  Top 5 genes:"
    head -n 5 $projPath/peakCalling/SEACR/${day}_PAX8_summed.txt
    echo ""
done

echo "Gene signal extraction complete"



###########Detecting overlap between predicted and actual PAX8 binding sites in R#######