#!/bin/bash
#SBATCH --job-name=11_igv_tracks
#SBATCH --output=logs/11_create_igv_tracks.out
#SBATCH --error=logs/11_create_igv_tracks.out
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq

# =============================================================================
# Script: 11_create_igv_tracks.sh
# Purpose: SLURM wrapper for creating IGV BED tracks
# Usage: sbatch scripts/11_create_igv_tracks.sh
# Output: results/IGV_tracks/*.bed
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create logs directory if needed
mkdir -p logs

echo "========================================"
echo "Creating IGV BED Tracks"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

# Run the R script
Rscript scripts/11_create_igv_tracks.R

echo ""
echo "End time: $(date)"
echo "Done!"
