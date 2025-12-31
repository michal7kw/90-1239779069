#!/bin/bash
#SBATCH --job-name=3_bigwig_rerun
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3_bigwig_rerun_%a.err"
#SBATCH --output="./logs/3_bigwig_rerun_%a.out"

# =============================================================================
# BigWig Generation Script - RERUN for missing samples
# Project: 90-1239779069
# Missing samples: 2 (Parental_2), 5 (Neg_2), 7 (Pos_1), 15 (KO_3)
# =============================================================================

set -euo pipefail

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
INPUT_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/results/03_bigwig"

# Only the missing samples
MISSING_SAMPLES=(2 5 7 15)
SAMPLE=${MISSING_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Sample name mapping
declare -A SAMPLE_NAMES
SAMPLE_NAMES[2]="Parental_2"
SAMPLE_NAMES[5]="Neg_2"
SAMPLE_NAMES[7]="Pos_1"
SAMPLE_NAMES[15]="KO_3"

SAMPLE_NAME=${SAMPLE_NAMES[$SAMPLE]}

echo "============================================"
echo "BigWig Generation (RERUN) - Sample ${SAMPLE} (${SAMPLE_NAME})"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment (working deeptools)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/bigwig

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Define input BAM file
BAM_FILE="${INPUT_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

# Verify input file exists
if [[ ! -f "${BAM_FILE}" ]]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "Processing: ${BAM_FILE}"

# Generate CPM-normalized BigWig (main track for cross-sample comparison)
echo "Generating CPM-normalized track..."
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_CPM.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --extendReads \
    --ignoreDuplicates

# Generate strand-specific tracks (forward)
echo "Generating forward strand track..."
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_forward.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --extendReads \
    --ignoreDuplicates \
    --filterRNAstrand forward

# Generate strand-specific tracks (reverse)
echo "Generating reverse strand track..."
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_reverse.bw \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --extendReads \
    --ignoreDuplicates \
    --filterRNAstrand reverse

# Generate RPKM-normalized track (alternative normalization)
echo "Generating RPKM-normalized track..."
bamCoverage \
    --bam ${BAM_FILE} \
    --outFileName ${OUTPUT_DIR}/${SAMPLE_NAME}_RPKM.bw \
    --outFileFormat bigwig \
    --normalizeUsing RPKM \
    --binSize 10 \
    --numberOfProcessors 16 \
    --extendReads \
    --ignoreDuplicates

echo "============================================"
echo "BigWig generation completed for sample ${SAMPLE} (${SAMPLE_NAME})"
echo "Output files:"
echo "  ${OUTPUT_DIR}/${SAMPLE_NAME}_CPM.bw"
echo "  ${OUTPUT_DIR}/${SAMPLE_NAME}_forward.bw"
echo "  ${OUTPUT_DIR}/${SAMPLE_NAME}_reverse.bw"
echo "  ${OUTPUT_DIR}/${SAMPLE_NAME}_RPKM.bw"
echo "End time: $(date)"
echo "============================================"
