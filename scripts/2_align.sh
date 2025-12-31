#!/bin/bash
#SBATCH --job-name=2_align
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-11
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/2_align_%a.err"
#SBATCH --output="./logs/2_align_%a.out"

# =============================================================================
# STAR Alignment Script
# Project: 90-1239779069
# Samples: Parental (1-3), Neg (4-6), Pos (7-9), KO (13-15)
# =============================================================================

set -euo pipefail

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
INPUT_DIR="${BASE_DIR}/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/02_aligned"
GENOME_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCm39_v2"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Load sample array from config
SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "============================================"
echo "STAR Alignment - Sample ${SAMPLE}"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

# Create output directory for this sample
SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p ${SAMPLE_OUTPUT_DIR}

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
echo "  Genome Index: ${GENOME_INDEX}"
echo "  GTF: ${GTF_FILE}"

# Run STAR alignment (two-pass mode)
STAR \
    --runThreadN 32 \
    --genomeDir ${GENOME_INDEX} \
    --sjdbGTFfile ${GTF_FILE} \
    --readFilesIn ${R1_FILE} ${R2_FILE} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --quantMode GeneCounts \
    --twopassMode Basic

echo "Indexing BAM file..."
samtools index -@ 16 ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam

echo "Generating alignment statistics..."
samtools flagstat ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_flagstat.txt
samtools stats ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_stats.txt

# Create convenient symlinks in the main output directory
ln -sf ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam ${OUTPUT_DIR}/${SAMPLE}_sorted.bam
ln -sf ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam.bai ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai
ln -sf ${SAMPLE_OUTPUT_DIR}/${SAMPLE}_ReadsPerGene.out.tab ${OUTPUT_DIR}/${SAMPLE}_counts.tab

echo "============================================"
echo "STAR alignment completed for sample ${SAMPLE}"
echo "End time: $(date)"
echo "============================================"
