#!/bin/bash
#SBATCH --job-name=7_sashimi
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/7_sashimi.err"
#SBATCH --output="./logs/7_sashimi.out"

# =============================================================================
# Sashimi Plot Generation
# Project: 90-1239779069
# Generates sashimi plots for top DE/DS genes and user-specified genes
#
# Features:
#   - Aggregated: Averages replicates per group (-A mean)
#   - Strict filtering: High filter (-M 100) shows only major isoforms
#   - Dual output: Both PDF and PNG formats
#   - Clean visualization with larger dimensions
# =============================================================================

set -euo pipefail
export PYTHONNOUSERSITE=1

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
BAM_DIR="${BASE_DIR}/results/02_aligned"
DESEQ_DIR="${BASE_DIR}/results/04_deseq2"
SPLICING_DIR="${BASE_DIR}/results/05_splicing"
OUTPUT_DIR="${BASE_DIR}/results/06_sashimi"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# User-specified genes file (create this file with gene symbols, one per line)
USER_GENES="${BASE_DIR}/config/genes_of_interest.txt"

# Plot settings
MIN_JUNCTION_COV=100    # Strict filtering - show only major isoforms
PLOT_HEIGHT=5
PLOT_WIDTH=18
PNG_RESOLUTION=300

echo "============================================"
echo "Sashimi Plot Generation"
echo "Start time: $(date)"
echo "============================================"
echo ""
echo "Settings:"
echo "  - Aggregation: mean (replicates averaged per group)"
echo "  - Min junction coverage: ${MIN_JUNCTION_COV}"
echo "  - Output formats: PDF + PNG (${PNG_RESOLUTION} DPI)"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi

# Create output directories
mkdir -p ${OUTPUT_DIR}/{top_DE,top_splicing,user_genes}

# Create BAM configuration file for ggsashimi
# Using numbered prefixes for proper group ordering in plots
cat > ${OUTPUT_DIR}/bam_config.tsv << EOF
Parental_1	${BAM_DIR}/1/1_Aligned.sortedByCoord.out.bam	1_Parental
Parental_2	${BAM_DIR}/2/2_Aligned.sortedByCoord.out.bam	1_Parental
Parental_3	${BAM_DIR}/3/3_Aligned.sortedByCoord.out.bam	1_Parental
Neg_1	${BAM_DIR}/4/4_Aligned.sortedByCoord.out.bam	2_Neg
Neg_2	${BAM_DIR}/5/5_Aligned.sortedByCoord.out.bam	2_Neg
Neg_3	${BAM_DIR}/6/6_Aligned.sortedByCoord.out.bam	2_Neg
Pos_1	${BAM_DIR}/7/7_Aligned.sortedByCoord.out.bam	3_Pos
Pos_2	${BAM_DIR}/8/8_Aligned.sortedByCoord.out.bam	3_Pos
Pos_3	${BAM_DIR}/9/9_Aligned.sortedByCoord.out.bam	3_Pos
KO_1	${BAM_DIR}/13/13_Aligned.sortedByCoord.out.bam	4_KO
KO_2	${BAM_DIR}/14/14_Aligned.sortedByCoord.out.bam	4_KO
KO_3	${BAM_DIR}/15/15_Aligned.sortedByCoord.out.bam	4_KO
EOF


# Export environment variables for Python script
export BASE_DIR OUTPUT_DIR DESEQ_DIR SPLICING_DIR GTF_FILE USER_GENES
export MIN_JUNCTION_COV PLOT_HEIGHT PLOT_WIDTH PNG_RESOLUTION

# Run the Python script
python3 7_sashimi.py

echo ""
echo "============================================"
echo "Sashimi plot generation complete!"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Settings used:"
echo "  - Aggregation: mean (replicates averaged per group)"
echo "  - Min junction coverage: ${MIN_JUNCTION_COV}"
echo "  - Plot dimensions: ${PLOT_WIDTH}x${PLOT_HEIGHT}"
echo "  - PNG resolution: ${PNG_RESOLUTION} DPI"
echo ""
echo "Subdirectories:"
echo "  top_DE/       - Top differentially expressed genes (PDF + PNG)"
echo "  top_splicing/ - Top differential splicing events (PDF + PNG)"
echo "  user_genes/   - User-specified genes (PDF + PNG)"
echo ""
echo "To add custom genes, create/edit:"
echo "  ${USER_GENES}"
echo "(one gene symbol per line)"
echo ""
echo "End time: $(date)"
echo "============================================"
