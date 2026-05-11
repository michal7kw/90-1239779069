#!/bin/bash
#SBATCH --job-name=10_sig_events
#SBATCH --output=logs/10_extract_significant_events.out
#SBATCH --error=logs/10_extract_significant_events.out
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq

# =============================================================================
# Script: 10_extract_significant_events.sh
# Purpose: SLURM wrapper for extracting significant splicing events
# Usage: sbatch scripts/10_extract_significant_events.sh
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create logs directory if needed
mkdir -p logs

echo "========================================"
echo "Extracting Significant Splicing Events"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

# Run the R script
Rscript scripts/10_extract_significant_events.R

echo ""
echo "End time: $(date)"
echo "Done!"
