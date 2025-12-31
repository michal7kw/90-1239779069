#!/bin/bash
#SBATCH --job-name=build_star_mouse
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/0_build_star_index.err"
#SBATCH --output="./logs/0_build_star_index.out"

# =============================================================================
# Build STAR Index for Mouse GRCm39
# Uses GENCODE vM33 annotation
# =============================================================================

set -euo pipefail

# Configuration
GENOME_FASTA="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/fasta/genome.fa"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"
OUTPUT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCm39_v2"

echo "============================================"
echo "Building STAR Index for Mouse GRCm39"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

# Check STAR version
echo "STAR version:"
STAR --version

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo ""
echo "Input files:"
echo "  Genome FASTA: ${GENOME_FASTA}"
echo "  GTF: ${GTF_FILE}"
echo "  Output: ${OUTPUT_DIR}"
echo ""

# Build STAR index
STAR \
    --runMode genomeGenerate \
    --runThreadN 32 \
    --genomeDir ${OUTPUT_DIR} \
    --genomeFastaFiles ${GENOME_FASTA} \
    --sjdbGTFfile ${GTF_FILE} \
    --sjdbOverhang 100

echo ""
echo "============================================"
echo "STAR index generation complete!"
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"
echo "============================================"
