#!/bin/bash
#SBATCH --job-name=14d_srrm
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14d_srrm_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14d_srrm_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

echo "============================================"
echo "14d: Srrm3/Srrm4 Expression"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

Rscript scripts/14d_srrm_expression.R

if [ $? -eq 0 ]; then
    echo ""
    echo "Results in: results/14_todo_analysis/task3_srrm_expression/"
    echo "End time: $(date)"
else
    echo "ERROR: 14d failed!"
    exit 1
fi
