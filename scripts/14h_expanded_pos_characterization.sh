#!/bin/bash
#SBATCH --job-name=14h_expanded_char
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14h_expanded_pos_characterization_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14h_expanded_pos_characterization_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=workq

echo "============================================"
echo "14h: Expanded Pos Characterization (3 Groups)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

# Validate 14g output exists
INPUT_BASE="results/14_todo_analysis/task5_expanded_groups"
for group in group1_neg_microexon group2_neg_micro_small group3_neg_all_sizes; do
    if [ ! -f "${INPUT_BASE}/${group}/pos_all_events_for_genes.csv" ]; then
        echo "ERROR: Missing ${INPUT_BASE}/${group}/pos_all_events_for_genes.csv"
        echo "Please run 14g first."
        exit 1
    fi
done

echo "All input files found. Running analysis..."
echo ""

Rscript scripts/14h_expanded_pos_characterization.R

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "=== 14h completed successfully ==="
    for group in group1_neg_microexon group2_neg_micro_small group3_neg_all_sizes; do
        echo ""
        echo "Output: ${INPUT_BASE}/${group}/pos_significant_deep/"
        ls -la "${INPUT_BASE}/${group}/pos_significant_deep/"
    done
else
    echo "=== 14h FAILED with exit code $EXIT_CODE ==="
fi

echo "End time: $(date)"
exit $EXIT_CODE
