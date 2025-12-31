#!/bin/bash
#SBATCH --job-name=5_deseq2
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5_deseq2.err"
#SBATCH --output="./logs/5_deseq2.out"

# =============================================================================
# DESeq2 Differential Expression Analysis - SLURM Wrapper
# Project: 90-1239779069
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"

echo "============================================"
echo "DESeq2 Analysis"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rna_seq_analysis

# Run DESeq2 R script
Rscript ${BASE_DIR}/scripts/5_deseq2.R

echo "============================================"
echo "DESeq2 analysis complete"
echo "End time: $(date)"
echo "============================================"
