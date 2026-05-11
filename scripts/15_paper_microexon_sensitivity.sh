#!/bin/bash
#SBATCH --job-name=15_sensitivity
#SBATCH --output=logs/15_paper_microexon_sensitivity.out
#SBATCH --error=logs/15_paper_microexon_sensitivity.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=workq

# =============================================================================
# Script: 15_paper_microexon_sensitivity.sh
# Purpose: SLURM wrapper for paper microexon sensitivity cross-reference
# Usage: sbatch scripts/15_paper_microexon_sensitivity.sh
#
# Cross-references 234 microexons classified by sensitivity tier (HS, LS, CS,
# CR, NonResponding) from a 2025 NSMB paper with our SRRM3 significant
# splicing events for Neg/KO/Pos vs Parental comparisons.
#
# Output: results/15_paper_microexon_sensitivity/
# =============================================================================

set -e

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create directories
mkdir -p logs
mkdir -p results/15_paper_microexon_sensitivity

echo "========================================"
echo "Paper Microexon Sensitivity Analysis"
echo "NSMB 2025 - Sensitivity Tier Cross-Ref"
echo "========================================"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo ""

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

# Run the R script
Rscript scripts/15_paper_microexon_sensitivity.R

echo ""
echo "========================================"
echo "Done: $(date)"
echo "========================================"
