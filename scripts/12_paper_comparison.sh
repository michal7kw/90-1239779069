#!/bin/bash
#SBATCH --job-name=12_paper_cmp
#SBATCH --output=logs/12_paper_comparison.out
#SBATCH --error=logs/12_paper_comparison.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

# =============================================================================
# Script: 12_paper_comparison.sh
# Purpose: SLURM wrapper for paper comparison analysis
# Usage: sbatch scripts/12_paper_comparison.sh
#
# This script compares SRRM3 splicing data with Gonatopoulos-Pournatzis 2018
# paper for all three vs-Parental comparisons:
#   - Neg_vs_Parental (partial SRRM3 knockdown)
#   - Pos_vs_Parental (SRRM3 overexpression)
#   - KO_vs_Parental  (SRRM3 knockout)
#
# Output:
#   results/12_paper_comparison/{comparison}/
#     - scatter_srrm3_vs_*.pdf, scatter_srrm4_vs_*.pdf
#     - scatter_combined.pdf, dpsi_heatmap.pdf, summary_bar.pdf
#     - paper_comparison_full.csv, matched_events.csv
#     - unmatched_events.csv, statistics_summary.csv
#     - multi_match_events.txt, multi_match_events.csv
#   results/12_paper_comparison/cross_comparison_statistics.csv
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create logs directory if needed
mkdir -p logs

echo "========================================"
echo "Paper Comparison Analysis"
echo "Gonatopoulos-Pournatzis et al. 2018"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

# Run the R script
Rscript scripts/12_paper_comparison.R

echo ""
echo "End time: $(date)"
echo "Done!"
