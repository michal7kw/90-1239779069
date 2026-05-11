#!/bin/bash
#SBATCH --job-name=14f_pos_sig_deep
#SBATCH --account=kubacki.michal
#SBATCH --output=logs/14f_neg_excl_pos_significant_%j.out
#SBATCH --error=logs/14f_neg_excl_pos_significant_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=workq

# Deep characterization of significant Pos events in Neg-exclusive microexon genes
# Requires: 14b + 14c to have been run

set -euo pipefail

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

echo "=== 14f: Pos Significant Deep Characterization ==="
echo "Start time: $(date)"
echo ""

# Validate required input files from 14b/14c
INPUT_DIR="results/14_todo_analysis/task2_neg_exclusive_microexon"

REQUIRED_FILES=(
    "${INPUT_DIR}/neg_exclusive_genes_pos_events.csv"
    "${INPUT_DIR}/neg_all_events_for_genes.csv"
    "${INPUT_DIR}/neg_exclusive_genes_pos_summary.csv"
    "${INPUT_DIR}/neg_exclusive_microexon_events.csv"
)

for f in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required input file not found: $f"
        echo "Please run 14b and 14c first."
        exit 1
    fi
done

echo "All input files found. Running analysis..."
echo ""

Rscript scripts/14f_neg_excl_pos_significant.R

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "=== 14f completed successfully ==="
    echo "Output: ${INPUT_DIR}/pos_significant_deep/"
    ls -la "${INPUT_DIR}/pos_significant_deep/"
else
    echo "=== 14f FAILED with exit code $EXIT_CODE ==="
fi

echo "End time: $(date)"
exit $EXIT_CODE
