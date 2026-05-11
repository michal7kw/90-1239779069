#!/bin/bash
#SBATCH --job-name=14g_se_expand
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14g_se_expanded_gene_selection_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14g_se_expanded_gene_selection_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

echo "============================================"
echo "14g_se: Expanded Gene Selection (SE only)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

Rscript scripts/14g_se_expanded_gene_selection.R

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "Results in: results/14_todo_analysis/task5_expanded_groups_se/"
    ls -la results/14_todo_analysis/task5_expanded_groups_se/*/
else
    echo "=== 14g_se FAILED with exit code $EXIT_CODE ==="
fi

echo "End time: $(date)"
exit $EXIT_CODE
