#!/bin/bash
#SBATCH --job-name=igv_subset
#SBATCH --output=logs/igv_subset_%j.log
#SBATCH --error=logs/igv_subset_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# =============================================================================
# Script: 9_create_igv_subset.sh
# Purpose: Create subset BigWig and BAM files for efficient IGV visualization
#          This allows downloading smaller files while retaining full
#          functionality including interactive sashimi plot generation.
# =============================================================================

set -e  # Exit on error

# Directories
PROJ_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
BED_DIR="${PROJ_DIR}/results/07_igv/bed"
BAM_DIR="${PROJ_DIR}/results/02_aligned"
BIGWIG_DIR="${PROJ_DIR}/results/03_bigwig"
OUTPUT_DIR="${PROJ_DIR}/results/08_igv_subset"
CHROM_SIZES="${OUTPUT_DIR}/mm39.chrom.sizes"

# Create output directories
mkdir -p ${OUTPUT_DIR}/{regions,bigwig,bam}

# Initialize conda and activate bigwig environment (has all required tools)
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/bigwig

# Verify tools are available
echo "Using samtools: $(which samtools)"
echo "Using bedtools: $(which bedtools)"
echo "Using bigWigToBedGraph: $(which bigWigToBedGraph)"
echo "Using bedGraphToBigWig: $(which bedGraphToBigWig)"

echo "=============================================="
echo "Starting IGV Subset Creation"
echo "Date: $(date)"
echo "=============================================="

# -----------------------------------------------------------------------------
# Step 0: Generate chromosome sizes from BAM header
# -----------------------------------------------------------------------------
echo ""
echo "Step 0: Generating chromosome sizes file..."

# Extract from any BAM header
samtools view -H ${BAM_DIR}/1/1_Aligned.sortedByCoord.out.bam | \
    grep "@SQ" | \
    sed 's/@SQ\tSN:\([^\t]*\)\tLN:\([0-9]*\)/\1\t\2/' > ${CHROM_SIZES}

echo "Created ${CHROM_SIZES} with $(wc -l < ${CHROM_SIZES}) chromosomes"

# -----------------------------------------------------------------------------
# Step 1: Merge all significant splicing regions and expand
# -----------------------------------------------------------------------------
echo ""
echo "Step 1: Creating merged regions files..."

# For BigWig: expand by 50kb for genomic context
cat ${BED_DIR}/*_significant_splicing.bed | \
    grep -v "^track" | \
    cut -f1-3 | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin | \
    awk -v OFS="\t" '{
        start = ($2 - 50000 < 0) ? 0 : $2 - 50000;
        print $1, start, $3 + 50000
    }' | \
    bedtools merge -i stdin > ${OUTPUT_DIR}/regions/splicing_regions_expanded.bed

echo "BigWig regions: $(wc -l < ${OUTPUT_DIR}/regions/splicing_regions_expanded.bed) merged regions (±50kb)"

# For BAM: smaller padding (5kb is enough for sashimi plots)
cat ${BED_DIR}/*_significant_splicing.bed | \
    grep -v "^track" | \
    cut -f1-3 | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i stdin | \
    awk -v OFS="\t" '{
        start = ($2 - 5000 < 0) ? 0 : $2 - 5000;
        print $1, start, $3 + 5000
    }' | \
    bedtools merge -i stdin > ${OUTPUT_DIR}/regions/splicing_regions_bam.bed

echo "BAM regions: $(wc -l < ${OUTPUT_DIR}/regions/splicing_regions_bam.bed) merged regions (±5kb)"

# -----------------------------------------------------------------------------
# Step 2: Create subset BigWig files
# -----------------------------------------------------------------------------
echo ""
echo "Step 2: Creating subset BigWig files..."

# Check if bigWigToBedGraph and bedGraphToBigWig are available
if ! command -v bigWigToBedGraph &> /dev/null || ! command -v bedGraphToBigWig &> /dev/null; then
    echo "WARNING: UCSC tools not found in PATH"
    echo "Available in conda: $(find /beegfs/scratch/ric.broccoli/kubacki.michal/conda -name 'bigWigToBedGraph' -type f 2>/dev/null | head -1)"
fi
echo "Using bigWigToBedGraph: $(which bigWigToBedGraph 2>/dev/null || echo 'not found')"
echo "Using bedGraphToBigWig: $(which bedGraphToBigWig 2>/dev/null || echo 'not found')"

# Process each CPM-normalized BigWig file
for bw in ${BIGWIG_DIR}/*_CPM.bw; do
    if [[ -f "$bw" ]]; then
        sample=$(basename $bw _CPM.bw)
        echo "  Processing ${sample}..."

        # Method 1: Using bigWigToBedGraph + bedtools + bedGraphToBigWig
        if command -v bigWigToBedGraph &> /dev/null && command -v bedGraphToBigWig &> /dev/null; then
            bigWigToBedGraph $bw /dev/stdout | \
                bedtools intersect -a stdin -b ${OUTPUT_DIR}/regions/splicing_regions_expanded.bed > \
                ${OUTPUT_DIR}/bigwig/${sample}_subset.bedGraph

            # Sort bedGraph (required for bedGraphToBigWig)
            sort -k1,1 -k2,2n ${OUTPUT_DIR}/bigwig/${sample}_subset.bedGraph > \
                ${OUTPUT_DIR}/bigwig/${sample}_subset.sorted.bedGraph

            bedGraphToBigWig ${OUTPUT_DIR}/bigwig/${sample}_subset.sorted.bedGraph \
                ${CHROM_SIZES} \
                ${OUTPUT_DIR}/bigwig/${sample}_subset.bw

            # Clean up intermediate files
            rm -f ${OUTPUT_DIR}/bigwig/${sample}_subset.bedGraph \
                  ${OUTPUT_DIR}/bigwig/${sample}_subset.sorted.bedGraph
        else
            echo "    WARNING: UCSC tools not available, copying full BigWig"
            cp $bw ${OUTPUT_DIR}/bigwig/${sample}_subset.bw
        fi
    fi
done

echo "BigWig files created:"
ls -lh ${OUTPUT_DIR}/bigwig/*.bw 2>/dev/null | awk '{print "  " $9 ": " $5}' || echo "  No BigWig files found"

# -----------------------------------------------------------------------------
# Step 3: Create subset BAM files for sashimi plots
# -----------------------------------------------------------------------------
echo ""
echo "Step 3: Creating subset BAM files..."

# Sample mapping: sample_number -> condition_replicate
declare -A SAMPLE_MAP=(
    ["1"]="Parental_1"
    ["2"]="Parental_2"
    ["3"]="Parental_3"
    ["4"]="Neg_1"
    ["5"]="Neg_2"
    ["6"]="Neg_3"
    ["7"]="Pos_1"
    ["8"]="Pos_2"
    ["9"]="Pos_3"
    ["13"]="KO_1"
    ["14"]="KO_2"
    ["15"]="KO_3"
)

for sample_num in "${!SAMPLE_MAP[@]}"; do
    sample_name="${SAMPLE_MAP[$sample_num]}"
    bam_path="${BAM_DIR}/${sample_num}/${sample_num}_Aligned.sortedByCoord.out.bam"

    if [[ -f "$bam_path" ]]; then
        echo "  Processing ${sample_name} (sample ${sample_num})..."

        # Extract reads from regions of interest
        samtools view -b -h -L ${OUTPUT_DIR}/regions/splicing_regions_bam.bed \
            -@ 8 \
            $bam_path > ${OUTPUT_DIR}/bam/${sample_name}_subset.bam

        # Index the subset BAM
        samtools index ${OUTPUT_DIR}/bam/${sample_name}_subset.bam
    else
        echo "  WARNING: BAM file not found: $bam_path"
    fi
done

echo ""
echo "BAM files created:"
ls -lh ${OUTPUT_DIR}/bam/*.bam 2>/dev/null | awk '{print "  " $9 ": " $5}' || echo "  No BAM files found"

# -----------------------------------------------------------------------------
# Step 4: Create IGV session file
# -----------------------------------------------------------------------------
echo ""
echo "Step 4: Creating IGV session file..."

cat > ${OUTPUT_DIR}/splicing_subset_session.xml << 'XMLEOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <!-- Coverage tracks (BigWig) - grouped by condition -->
        <!-- Parental -->
        <Resource path="bigwig/Parental_1_subset.bw"/>
        <Resource path="bigwig/Parental_2_subset.bw"/>
        <Resource path="bigwig/Parental_3_subset.bw"/>
        <!-- Neg -->
        <Resource path="bigwig/Neg_1_subset.bw"/>
        <Resource path="bigwig/Neg_2_subset.bw"/>
        <Resource path="bigwig/Neg_3_subset.bw"/>
        <!-- Pos -->
        <Resource path="bigwig/Pos_1_subset.bw"/>
        <Resource path="bigwig/Pos_2_subset.bw"/>
        <Resource path="bigwig/Pos_3_subset.bw"/>
        <!-- KO -->
        <Resource path="bigwig/KO_1_subset.bw"/>
        <Resource path="bigwig/KO_2_subset.bw"/>
        <Resource path="bigwig/KO_3_subset.bw"/>

        <!-- BAM files for sashimi plots - grouped by condition -->
        <!-- Parental -->
        <Resource path="bam/Parental_1_subset.bam"/>
        <Resource path="bam/Parental_2_subset.bam"/>
        <Resource path="bam/Parental_3_subset.bam"/>
        <!-- Neg -->
        <Resource path="bam/Neg_1_subset.bam"/>
        <Resource path="bam/Neg_2_subset.bam"/>
        <Resource path="bam/Neg_3_subset.bam"/>
        <!-- Pos -->
        <Resource path="bam/Pos_1_subset.bam"/>
        <Resource path="bam/Pos_2_subset.bam"/>
        <Resource path="bam/Pos_3_subset.bam"/>
        <!-- KO -->
        <Resource path="bam/KO_1_subset.bam"/>
        <Resource path="bam/KO_2_subset.bam"/>
        <Resource path="bam/KO_3_subset.bam"/>

        <!-- Splicing event BED tracks -->
        <Resource path="bed/KO_vs_Parental_significant_splicing.bed"/>
        <Resource path="bed/Neg_vs_Parental_significant_splicing.bed"/>
        <Resource path="bed/Pos_vs_Parental_significant_splicing.bed"/>
        <Resource path="bed/KO_vs_Neg_significant_splicing.bed"/>
        <Resource path="bed/Pos_vs_Neg_significant_splicing.bed"/>
        <Resource path="bed/top_regions.bed"/>
    </Resources>
</Session>
XMLEOF

echo "IGV session created: ${OUTPUT_DIR}/splicing_subset_session.xml"

# -----------------------------------------------------------------------------
# Step 5: Copy BED files to output directory
# -----------------------------------------------------------------------------
echo ""
echo "Step 5: Copying BED files..."

mkdir -p ${OUTPUT_DIR}/bed
cp ${BED_DIR}/*.bed ${OUTPUT_DIR}/bed/
cp ${PROJ_DIR}/results/07_igv/top_splicing_regions.txt ${OUTPUT_DIR}/

echo "BED files copied: $(ls ${OUTPUT_DIR}/bed/*.bed | wc -l) files"

# -----------------------------------------------------------------------------
# Step 6: Create download package
# -----------------------------------------------------------------------------
echo ""
echo "Step 6: Creating download package..."

cd ${OUTPUT_DIR}

# Create tarball with all subset files (using relative paths)
tar -cvzf igv_subset_package.tar.gz \
    bigwig/*.bw \
    bam/*.bam \
    bam/*.bam.bai \
    bed/*.bed \
    splicing_subset_session.xml \
    top_splicing_regions.txt \
    mm39.chrom.sizes

echo ""
echo "=============================================="
echo "COMPLETE!"
echo "=============================================="
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Files created:"
du -sh ${OUTPUT_DIR}/*
echo ""
echo "Download package:"
ls -lh ${OUTPUT_DIR}/igv_subset_package.tar.gz
echo ""
echo "To download:"
echo "  scp user@server:${OUTPUT_DIR}/igv_subset_package.tar.gz ."
echo ""
echo "After extraction, open splicing_subset_session.xml in IGV (with mm39 genome)."
echo "=============================================="
