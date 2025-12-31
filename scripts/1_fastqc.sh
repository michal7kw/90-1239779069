#!/bin/bash
#SBATCH --job-name=1_fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-11
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/1_fastqc_%a.err"
#SBATCH --output="./logs/1_fastqc_%a.out"

# =============================================================================
# FastQC Quality Control Script
# Project: 90-1239779069
# Samples: Parental (1-3), Neg (4-6), Pos (7-9), KO (13-15)
# =============================================================================

set -euo pipefail

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
INPUT_DIR="${BASE_DIR}/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/01_fastqc"

# Load sample array from config
SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "============================================"
echo "FastQC Analysis - Sample ${SAMPLE}"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/quality

# Create output directory for this sample
mkdir -p ${OUTPUT_DIR}/${SAMPLE}

# Define input files
R1_FILE="${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz"
R2_FILE="${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz"

# Verify input files exist
if [[ ! -f "${R1_FILE}" ]]; then
    echo "ERROR: R1 file not found: ${R1_FILE}"
    exit 1
fi

if [[ ! -f "${R2_FILE}" ]]; then
    echo "ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

echo "Processing files:"
echo "  R1: ${R1_FILE}"
echo "  R2: ${R2_FILE}"

# Run FastQC
fastqc \
    --threads 8 \
    --outdir ${OUTPUT_DIR}/${SAMPLE} \
    --format fastq \
    ${R1_FILE} ${R2_FILE}

echo "============================================"
echo "FastQC completed for sample ${SAMPLE}"
echo "End time: $(date)"
echo "============================================"
