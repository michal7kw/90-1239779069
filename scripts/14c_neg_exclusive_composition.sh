#!/bin/bash
#SBATCH --job-name=14c_composition
#SBATCH --account=kubacki.michal
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14c_composition_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/14c_composition_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=workq

echo "============================================"
echo "14c: Neg-Exclusive Microexon Gene Composition"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Requires: 14b output"
echo "Start time: $(date)"
echo "============================================"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate rna_seq_analysis

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069
mkdir -p logs

# Check that 14b output exists
if [ ! -f results/14_todo_analysis/task2_neg_exclusive_microexon/neg_exclusive_gene_ids.txt ]; then
    echo "ERROR: 14b output not found. Run 14b_neg_exclusive_genes.sh first."
    exit 1
fi

Rscript scripts/14c_neg_exclusive_composition.R

if [ $? -eq 0 ]; then
    echo ""
    echo "Results in: results/14_todo_analysis/task2_neg_exclusive_microexon/"
    echo "End time: $(date)"
else
    echo "ERROR: 14c failed!"
    exit 1
fi
