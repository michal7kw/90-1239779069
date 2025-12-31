#!/bin/bash
#SBATCH --job-name=4_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/4_multiqc.err"
#SBATCH --output="./logs/4_multiqc.out"

# =============================================================================
# MultiQC Report Aggregation Script
# Project: 90-1239779069
# Combines FastQC and STAR alignment reports into a single HTML report
# =============================================================================

set -euo pipefail

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"

echo "============================================"
echo "MultiQC Report Generation"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/multiqc

# Run MultiQC to aggregate all QC reports
multiqc \
    --force \
    --outdir ${BASE_DIR}/results \
    --filename multiqc_report \
    --title "RNA-seq Analysis - Project 90-1239779069" \
    --comment "Samples: Parental (1-3), Neg (4-6), Pos (7-9), KO (13-15)" \
    ${BASE_DIR}/results/01_fastqc \
    ${BASE_DIR}/results/02_aligned

echo "============================================"
echo "MultiQC report generated"
echo "Output: ${BASE_DIR}/results/multiqc_report.html"
echo "End time: $(date)"
echo "============================================"
