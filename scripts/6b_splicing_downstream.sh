#!/bin/bash
#SBATCH --job-name=6b_splicing_downstream
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6b_splicing_downstream.err"
#SBATCH --output="./logs/6b_splicing_downstream.out"

# =============================================================================
# Splicing Downstream Analysis Pipeline
# Project: 90-1239779069
# Scripts: 6b_splicing_downstream.R, 6c_microexon_analysis.R, 6d_splicing_enrichment.R
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SCRIPTS_DIR="${BASE_DIR}/scripts"

echo "============================================"
echo "Splicing Downstream Analysis Pipeline"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/splicing_analysis

# Create log directory if not exists
mkdir -p ${BASE_DIR}/logs

# =============================================================================
# Step 1: Splicing Downstream Analysis (PCA, distributions, directionality)
# =============================================================================
echo ""
echo "============================================"
echo "Step 1: Splicing Downstream Analysis"
echo "============================================"

Rscript ${SCRIPTS_DIR}/6b_splicing_downstream.R

if [ $? -eq 0 ]; then
    echo "Step 1 completed successfully"
else
    echo "Step 1 failed" >&2
    exit 1
fi

# =============================================================================
# Step 2: Microexon Analysis
# =============================================================================
echo ""
echo "============================================"
echo "Step 2: Microexon Analysis"
echo "============================================"

Rscript ${SCRIPTS_DIR}/6c_microexon_analysis.R

if [ $? -eq 0 ]; then
    echo "Step 2 completed successfully"
else
    echo "Step 2 failed" >&2
    exit 1
fi

# =============================================================================
# Step 3: GO/GSEA Enrichment Analysis
# =============================================================================
echo ""
echo "============================================"
echo "Step 3: GO/GSEA Enrichment Analysis"
echo "============================================"

Rscript ${SCRIPTS_DIR}/6d_splicing_enrichment.R

if [ $? -eq 0 ]; then
    echo "Step 3 completed successfully"
else
    echo "Step 3 failed" >&2
    exit 1
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================"
echo "Pipeline Complete!"
echo "============================================"
echo ""
echo "Output directories:"
echo "  - results/06_splicing_downstream/  (PCA, distributions, cross-plots)"
echo "  - results/07_microexon_analysis/   (Microexon-stratified analysis)"
echo "  - results/08_enrichment_analysis/  (GO/GSEA results)"
echo ""
echo "End time: $(date)"
echo "============================================"
