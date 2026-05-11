#!/bin/bash
#SBATCH --job-name=14i_se_negonly
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14i_se_negonly_expanded_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14i_se_negonly_expanded_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

echo "============================================"
echo "14i_se: Neg-Only Expanded Groups (SE only)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

Rscript scripts/14i_se_negonly_expanded.R

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "=== 14i_se completed successfully ==="
    for group in group2b_negonly_micro_small group3b_negonly_all_sizes; do
        echo ""
        echo "Output: results/14_todo_analysis/task5_expanded_groups_se/${group}/"
        ls -la "results/14_todo_analysis/task5_expanded_groups_se/${group}/"
        echo "Deep analysis:"
        ls -la "results/14_todo_analysis/task5_expanded_groups_se/${group}/pos_significant_deep/"
    done
else
    echo "=== 14i_se FAILED with exit code $EXIT_CODE ==="
fi

echo "End time: $(date)"
exit $EXIT_CODE
