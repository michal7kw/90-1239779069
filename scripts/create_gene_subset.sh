#!/bin/bash
# =============================================================================
# Create Gene-Specific IGV Subset
# =============================================================================
# Creates a small subset of BAM/BigWig files for specific genes only.
# Much smaller than the full splicing subset (~50-200 MB vs 10 GB).
#
# Usage:
#   ./create_gene_subset.sh <genes_file> [output_dir]
#
# genes_file format (tab-separated):
#   gene_name  chr  start  end
#   Myc        chr15  61800000  61900000
#   Tp53       chr11  69380000  69450000
#
# Or just gene names (will look up coordinates from BED files):
#   Myc
#   Tp53
#   Foxo1
# =============================================================================

set -e

# Configuration
PROJ_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
BAM_DIR="${PROJ_DIR}/results/02_aligned"
BIGWIG_DIR="${PROJ_DIR}/results/03_bigwig"
BED_DIR="${PROJ_DIR}/results/07_igv/bed"

GENES_FILE="${1:-genes.txt}"
OUTPUT_DIR="${2:-${PROJ_DIR}/results/gene_subset}"

# Initialize conda
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/bigwig

echo "=============================================="
echo "Creating Gene-Specific IGV Subset"
echo "=============================================="
echo "Input genes: ${GENES_FILE}"
echo "Output dir: ${OUTPUT_DIR}"
echo ""

# Create output directories
mkdir -p ${OUTPUT_DIR}/{regions,bigwig,bam,bed}

# Step 1: Create regions file
echo "Step 1: Creating regions file..."

if [[ $(head -1 "${GENES_FILE}" | awk -F'\t' '{print NF}') -ge 4 ]]; then
    # File has coordinates
    echo "  Using provided coordinates..."
    awk -F'\t' -v OFS="\t" '{
        start = ($3 - 10000 < 0) ? 0 : $3 - 10000;
        end = $4 + 10000;
        print $2, start, end, $1
    }' "${GENES_FILE}" > ${OUTPUT_DIR}/regions/genes.bed
else
    # File has only gene names - look up in BED files
    echo "  Looking up gene coordinates from splicing BED files..."

    # Extract gene names
    GENES=$(cat "${GENES_FILE}" | tr '\n' '|' | sed 's/|$//')

    # Search in all significant splicing BED files
    grep -h -E "(${GENES})" ${BED_DIR}/*_significant_splicing.bed 2>/dev/null | \
        grep -v "^track" | \
        awk -F'\t' -v OFS="\t" '{
            split($4, a, "_");
            gene = a[1];
            start = ($2 - 10000 < 0) ? 0 : $2 - 10000;
            end = $3 + 10000;
            print $1, start, end, gene
        }' | \
        sort -k1,1 -k2,2n | \
        bedtools merge -i stdin -c 4 -o distinct > ${OUTPUT_DIR}/regions/genes.bed
fi

N_REGIONS=$(wc -l < ${OUTPUT_DIR}/regions/genes.bed)
echo "  Created ${N_REGIONS} regions"

if [[ ${N_REGIONS} -eq 0 ]]; then
    echo "ERROR: No regions found for specified genes!"
    echo "Make sure gene names match those in the BED files."
    exit 1
fi

# Generate chrom.sizes
samtools view -H ${BAM_DIR}/1/1_Aligned.sortedByCoord.out.bam | \
    grep "@SQ" | \
    sed 's/@SQ\tSN:\([^\t]*\)\tLN:\([0-9]*\)/\1\t\2/' > ${OUTPUT_DIR}/mm39.chrom.sizes

# Step 2: Create subset BigWig files
echo ""
echo "Step 2: Creating subset BigWig files..."

for bw in ${BIGWIG_DIR}/*_CPM.bw; do
    sample=$(basename $bw _CPM.bw)
    echo "  Processing ${sample}..."

    bigWigToBedGraph $bw /dev/stdout | \
        bedtools intersect -a stdin -b ${OUTPUT_DIR}/regions/genes.bed > \
        ${OUTPUT_DIR}/bigwig/${sample}_subset.bedGraph

    sort -k1,1 -k2,2n ${OUTPUT_DIR}/bigwig/${sample}_subset.bedGraph > \
        ${OUTPUT_DIR}/bigwig/${sample}_subset.sorted.bedGraph

    bedGraphToBigWig ${OUTPUT_DIR}/bigwig/${sample}_subset.sorted.bedGraph \
        ${OUTPUT_DIR}/mm39.chrom.sizes \
        ${OUTPUT_DIR}/bigwig/${sample}_subset.bw

    rm -f ${OUTPUT_DIR}/bigwig/${sample}_subset.bedGraph \
          ${OUTPUT_DIR}/bigwig/${sample}_subset.sorted.bedGraph
done

# Step 3: Create subset BAM files
echo ""
echo "Step 3: Creating subset BAM files..."

declare -A SAMPLE_MAP=(
    ["1"]="Parental_1" ["2"]="Parental_2" ["3"]="Parental_3"
    ["4"]="Neg_1" ["5"]="Neg_2" ["6"]="Neg_3"
    ["7"]="Pos_1" ["8"]="Pos_2" ["9"]="Pos_3"
    ["13"]="KO_1" ["14"]="KO_2" ["15"]="KO_3"
)

for sample_num in "${!SAMPLE_MAP[@]}"; do
    sample_name="${SAMPLE_MAP[$sample_num]}"
    bam_path="${BAM_DIR}/${sample_num}/${sample_num}_Aligned.sortedByCoord.out.bam"

    if [[ -f "$bam_path" ]]; then
        echo "  Processing ${sample_name}..."
        samtools view -b -h -L ${OUTPUT_DIR}/regions/genes.bed \
            -@ 4 $bam_path > ${OUTPUT_DIR}/bam/${sample_name}_subset.bam
        samtools index ${OUTPUT_DIR}/bam/${sample_name}_subset.bam
    fi
done

# Step 4: Copy relevant BED annotations
echo ""
echo "Step 4: Extracting relevant BED annotations..."

for bed in ${BED_DIR}/*_significant_splicing.bed; do
    name=$(basename $bed)
    # Extract header and matching regions
    grep "^track" $bed > ${OUTPUT_DIR}/bed/${name}
    bedtools intersect -a $bed -b ${OUTPUT_DIR}/regions/genes.bed -wa >> ${OUTPUT_DIR}/bed/${name}
done

# Step 5: Create IGV session
echo ""
echo "Step 5: Creating IGV session file..."

cat > ${OUTPUT_DIR}/gene_subset_session.xml << 'EOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <Resource path="bigwig/Parental_1_subset.bw"/>
        <Resource path="bigwig/Parental_2_subset.bw"/>
        <Resource path="bigwig/Parental_3_subset.bw"/>
        <Resource path="bigwig/Neg_1_subset.bw"/>
        <Resource path="bigwig/Neg_2_subset.bw"/>
        <Resource path="bigwig/Neg_3_subset.bw"/>
        <Resource path="bigwig/Pos_1_subset.bw"/>
        <Resource path="bigwig/Pos_2_subset.bw"/>
        <Resource path="bigwig/Pos_3_subset.bw"/>
        <Resource path="bigwig/KO_1_subset.bw"/>
        <Resource path="bigwig/KO_2_subset.bw"/>
        <Resource path="bigwig/KO_3_subset.bw"/>
        <Resource path="bam/Parental_1_subset.bam"/>
        <Resource path="bam/Parental_2_subset.bam"/>
        <Resource path="bam/Parental_3_subset.bam"/>
        <Resource path="bam/Neg_1_subset.bam"/>
        <Resource path="bam/Neg_2_subset.bam"/>
        <Resource path="bam/Neg_3_subset.bam"/>
        <Resource path="bam/Pos_1_subset.bam"/>
        <Resource path="bam/Pos_2_subset.bam"/>
        <Resource path="bam/Pos_3_subset.bam"/>
        <Resource path="bam/KO_1_subset.bam"/>
        <Resource path="bam/KO_2_subset.bam"/>
        <Resource path="bam/KO_3_subset.bam"/>
    </Resources>
</Session>
EOF

# Step 6: Create package
echo ""
echo "Step 6: Creating download package..."

cd ${OUTPUT_DIR}
tar -czf gene_subset_package.tar.gz \
    bigwig/*.bw \
    bam/*.bam \
    bam/*.bam.bai \
    bed/*.bed \
    regions/genes.bed \
    gene_subset_session.xml \
    mm39.chrom.sizes

echo ""
echo "=============================================="
echo "COMPLETE!"
echo "=============================================="
echo ""
echo "Output:"
du -sh ${OUTPUT_DIR}/*
echo ""
echo "Download package:"
ls -lh ${OUTPUT_DIR}/gene_subset_package.tar.gz
echo ""
echo "To download:"
echo "  scp user@server:${OUTPUT_DIR}/gene_subset_package.tar.gz ."
