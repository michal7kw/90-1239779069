#!/bin/bash
# Master script to submit downstream analysis jobs in the correct order

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"

# Ensure we are in the base directory so ./logs paths in SBATCH directives work
cd ${BASE_DIR}
mkdir -p logs

echo "============================================"
echo "Submitting Analysis Pipeline"
echo "============================================"

# 1. Splicing Downstream (Runs 6b, 6c, 6d)
echo "Submitting 6b_splicing_downstream.sh (Includes 6b, 6c, 6d)..."
JOB1=$(sbatch scripts/6b_splicing_downstream.sh)
JOB_ID_1=$(echo $JOB1 | awk '{print $4}')
echo "  -> Job ID: $JOB_ID_1"

# 2. Extended Microexon Analysis (Runs 6e, 6f, 6g)
# We make this dependent on Job 1 to ensure logical order and reduce I/O contention,
# although strict data dependency is minimal.
echo "Submitting 6h_microexon_extended.sh (Includes 6e, 6f, 6g)..."
# Using --dependency=afterok:JOB_ID ensures this runs only after the first job finishes successfully
JOB2=$(sbatch --dependency=afterok:$JOB_ID_1 scripts/6h_microexon_extended.sh)
JOB_ID_2=$(echo $JOB2 | awk '{print $4}')
echo "  -> Job ID: $JOB_ID_2 (Depends on $JOB_ID_1)"

echo "============================================"
echo "Jobs submitted successfully."
echo "Check status with: squeue -u $USER"
echo "============================================"
squeue -u $USER
