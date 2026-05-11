#!/bin/bash
#SBATCH --job-name=12_paper_cmp_dkd
#SBATCH --output=logs/12_paper_comparison_double_kd.out
#SBATCH --error=logs/12_paper_comparison_double_kd.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

# =============================================================================
# Script: 12_paper_comparison_double_kd.sh
# Purpose: SLURM wrapper for double knockdown paper comparison analysis
# Usage: sbatch scripts/12_paper_comparison_double_kd.sh
#
# This script compares SRRM3 RNA-seq splicing data with Gonatopoulos-Pournatzis
# 2018 double knockdown (Srrm3+Srrm4) data from VastDB (54K events) for all
# three vs-Parental comparisons:
#   - Neg_vs_Parental (partial SRRM3 knockdown)
#   - Pos_vs_Parental (SRRM3 overexpression)
#   - KO_vs_Parental  (SRRM3 knockout)
#
# Output:
#   results/12_paper_comparison_double_kd/{comparison}/
#     - scatter_srrm3_srrm4_vs_*.pdf, scatter_srrm3_vs_*.pdf, scatter_srrm4_vs_*.pdf
#     - scatter_combined.pdf, dpsi_heatmap.pdf, summary_bar.pdf
#     - paper_comparison_full.csv, matched_events.csv
#     - unmatched_events.csv, statistics_summary.csv
#     - multi_match_events.txt, multi_match_events.csv
#   results/12_paper_comparison_double_kd/cross_comparison_statistics.csv
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create logs directory if needed
mkdir -p logs

echo "========================================"
echo "Paper Comparison Analysis - Double KD"
echo "Gonatopoulos-Pournatzis et al. 2018"
echo "Srrm3+Srrm4 Double Knockdown Data"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

# Run the R script
Rscript scripts/12_paper_comparison_double_kd.R

echo ""
echo "End time: $(date)"
echo "Done!"
