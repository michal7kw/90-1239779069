#!/bin/bash
#SBATCH --job-name=14a_venn
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14a_venn_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14a_venn_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

echo "============================================"
echo "14a: Venn Diagrams — Splicing Event Overlap"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

Rscript scripts/14a_venn_diagrams.R

if [ $? -eq 0 ]; then
    echo ""
    echo "Results in: results/14_todo_analysis/task1_venn_diagrams/"
    echo "End time: $(date)"
else
    echo "ERROR: 14a failed!"
    exit 1
fi
