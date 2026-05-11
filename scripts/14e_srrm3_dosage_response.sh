#!/bin/bash
#SBATCH --job-name=14e_dosage
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14e_dosage_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14e_dosage_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

echo "============================================"
echo "14e: SRRM3 Dosage-Microexon PSI Correlation"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

Rscript scripts/14e_srrm3_dosage_response.R

if [ $? -eq 0 ]; then
    echo ""
    echo "Results in: results/14_todo_analysis/task4_dosage_response/"
    echo "End time: $(date)"
else
    echo "ERROR: 14e failed!"
    exit 1
fi
