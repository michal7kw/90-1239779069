#!/bin/bash
#SBATCH --job-name=7c_sashimi
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/7c_sashimi.err"
#SBATCH --output="./logs/7c_sashimi.out"

# =============================================================================
# Minimal/Clean Sashimi Plot Generation
# Even cleaner version with:
#   - Higher junction threshold (100) to show only major splicing events
#   - SVG output for vector graphics (scalable)
#   - PNG output for easy viewing
# =============================================================================

set -euo pipefail
export PYTHONNOUSERSITE=1

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
BAM_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/results/06_sashimi/minimal"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

echo "============================================"
echo "Minimal Sashimi Plot Generation"
echo "Start time: $(date)"
echo "============================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi

mkdir -p ${OUTPUT_DIR}

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

generate_clean_plot() {
    local GENE_NAME=$1
    local REGION=$2
    local MIN_COV=${3:-100}  # Default 100, can override

    echo "Generating clean plot for: ${GENE_NAME} (${REGION}), min_cov=${MIN_COV}"

    # PDF version
    ggsashimi.py \
        -b ${OUTPUT_DIR}/bam_config.tsv \
        -c ${REGION} \
        -g ${GTF_FILE} \
        -o ${OUTPUT_DIR}/${GENE_NAME} \
        -M ${MIN_COV} \
        -C 3 \
        -O 3 \
        -A mean \
        --alpha 0.5 \
        --height 5 \
        --width 20 \
        --base-size 16 \
        --ann-height 4 \
        --shrink \
        -F pdf \
        2>&1 || echo "  Warning: PDF generation had issues"

    # PNG version (high resolution)
    ggsashimi.py \
        -b ${OUTPUT_DIR}/bam_config.tsv \
        -c ${REGION} \
        -g ${GTF_FILE} \
        -o ${OUTPUT_DIR}/${GENE_NAME}_hires \
        -M ${MIN_COV} \
        -C 3 \
        -O 3 \
        -A mean \
        --alpha 0.5 \
        --height 5 \
        --width 20 \
        --base-size 16 \
        --ann-height 4 \
        --shrink \
        -F png \
        -R 300 \
        2>&1 || echo "  Warning: PNG generation had issues"
}

# Generate plots - using higher thresholds for cleaner output
generate_clean_plot "Col6a5" "chr9:105732277-105838842" 100
generate_clean_plot "Atf3" "chr1:190901493-190951236" 100
generate_clean_plot "Arid5a" "chr1:36345814-36364110" 50
generate_clean_plot "Klhl24" "chr16:19915292-19948971" 50
generate_clean_plot "Chac1" "chr2:119180710-119185862" 50

echo ""
echo "============================================"
echo "Minimal sashimi plots complete!"
echo "Output: ${OUTPUT_DIR}"
echo "End time: $(date)"
echo "============================================"
