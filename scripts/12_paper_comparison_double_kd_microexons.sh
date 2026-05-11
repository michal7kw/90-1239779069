#!/bin/bash
#SBATCH --job-name=12_cmp_dkd_me
#SBATCH --output=logs/12_paper_comparison_double_kd_microexons.out
#SBATCH --error=logs/12_paper_comparison_double_kd_microexons.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

# =============================================================================
# Script: 12_paper_comparison_double_kd_microexons.sh
# Purpose: SLURM wrapper for double KD paper comparison — microexons only
# Usage: sbatch scripts/12_paper_comparison_double_kd_microexons.sh
#
# Microexon-filtered variant: only considers exons ≤30bp from the VastDB
# double knockdown dataset (~700 events vs ~27K in the full analysis).
#
# Output:
#   results/12_paper_comparison_double_kd_microexons/{comparison}/
#   results/12_paper_comparison_double_kd_microexons/cross_comparison_statistics.csv
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create logs directory if needed
mkdir -p logs

echo "========================================"
echo "Paper Comparison - Double KD Microexons"
echo "Gonatopoulos-Pournatzis et al. 2018"
echo "Exons ≤30bp only"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

# Run the R script
Rscript scripts/12_paper_comparison_double_kd_microexons.R

echo ""
echo "End time: $(date)"
echo "Done!"
