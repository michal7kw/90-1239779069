#!/bin/bash
#SBATCH --job-name=de_splicing_overlap
#SBATCH --output=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/de_splicing_overlap_%j.log
#SBATCH --error=/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/logs/de_splicing_overlap_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=workq

# =============================================================================
# GO/GSEA Enrichment for Genes with BOTH DE and Differential Splicing
# Project: 90-1239779069
# =============================================================================

echo "============================================"
echo "DE + Splicing Overlap Enrichment Analysis"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda/etc/profile.d/conda.sh
conda activate splicing_analysis

# Check if required packages are available
echo ""
echo "Checking R packages..."
Rscript -e "
suppressPackageStartupMessages({
    required <- c('clusterProfiler', 'org.Mm.eg.db', 'AnnotationDbi',
                  'enrichplot', 'ggplot2', 'dplyr', 'tidyr', 'DOSE', 'VennDiagram')
    missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) {
        cat('Missing packages:', paste(missing, collapse = ', '), '\n')
        quit(status = 1)
    } else {
        cat('All required packages are available\n')
    }
})
"

if [ $? -ne 0 ]; then
    echo "ERROR: Missing required R packages"
    exit 1
fi

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069

# Create logs directory if needed
mkdir -p logs

# Run the R script
echo ""
echo "Running enrichment analysis..."
Rscript scripts/6i_de_splicing_overlap_enrichment.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "============================================"
    echo "Analysis completed successfully!"
    echo "Results in: results/13_de_splicing_overlap_enrichment/"
    echo "End time: $(date)"
    echo "============================================"
else
    echo ""
    echo "ERROR: Analysis failed!"
    exit 1
fi
