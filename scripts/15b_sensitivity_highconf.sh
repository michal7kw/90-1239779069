#!/bin/bash
#SBATCH --job-name=15b_highconf
#SBATCH --output=logs/15b_sensitivity_highconf.out
#SBATCH --error=logs/15b_sensitivity_highconf.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

# =============================================================================
# Script: 15b_sensitivity_highconf.sh
# Purpose: SLURM wrapper for high-confidence (gene+length) subset analysis
# Usage: sbatch scripts/15b_sensitivity_highconf.sh
#
# Reads existing step 15 matched data and produces high-confidence subset
# CSVs and enhanced comparison plots.
#
# Output: results/15_paper_microexon_sensitivity/high_confidence/
#       + 3 new PDFs in results/15_paper_microexon_sensitivity/
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create directories
mkdir -p logs
mkdir -p results/15_paper_microexon_sensitivity/high_confidence

echo "========================================"
echo "High-Confidence Subset Analysis"
echo "Gene+Length Matches Only"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

# Run the R script
Rscript scripts/15b_sensitivity_highconf.R

echo ""
echo "========================================"
echo "Done: $(date)"
echo "========================================"
