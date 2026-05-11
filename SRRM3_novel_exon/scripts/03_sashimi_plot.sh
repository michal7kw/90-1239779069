#!/bin/bash
# Sashimi plot for novel SRRM3 exon
# Location: SRRM3_novel_exon/scripts/03_sashimi_plot.sh

set -e

ANALYSIS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/SRRM3_novel_exon"
BAM_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/results/02_aligned"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"
OUTPUT_DIR="${ANALYSIS_DIR}/results"
LOG_FILE="${ANALYSIS_DIR}/logs/03_sashimi_plot.log"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${ANALYSIS_DIR}/logs

exec > >(tee -a ${LOG_FILE}) 2>&1
echo "=== Sashimi plot generation started: $(date) ==="

# Create BAM config file for ggsashimi
cat > ${ANALYSIS_DIR}/bam_config.tsv << EOF
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

echo "Created BAM config file: ${ANALYSIS_DIR}/bam_config.tsv"

# Check if ggsashimi environment exists
GGSASHIMI_ENV="/beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/ggsashimi"
if [[ ! -d "${GGSASHIMI_ENV}" ]]; then
    echo "ERROR: ggsashimi conda environment not found at ${GGSASHIMI_ENV}"
    echo "Please create the environment first."
    exit 1
fi

# Activate conda and ggsashimi environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate ${GGSASHIMI_ENV}

echo "Generating sashimi plot for novel exon region..."
echo "Region: chr5:135897500-135902500"

# Generate sashimi plot
# -M 10: minimum junction count to show
# -C 3: color palette (3 = Set1)
# -O 3: overlay (3 = mean)
# --shrink: shrink introns
ggsashimi.py \
    -b ${ANALYSIS_DIR}/bam_config.tsv \
    -c chr5:135897500-135902500 \
    -g ${GTF_FILE} \
    -M 10 \
    -C 3 \
    -O 3 \
    -A mean \
    --alpha 0.5 \
    --height 3 \
    --width 10 \
    --shrink \
    -o ${OUTPUT_DIR}/sashimi_novel_exon

echo ""

# Also generate a zoomed-in view around the novel exon
echo "Generating zoomed sashimi plot..."
echo "Region: chr5:135898000-135902000"

ggsashimi.py \
    -b ${ANALYSIS_DIR}/bam_config.tsv \
    -c chr5:135898000-135902000 \
    -g ${GTF_FILE} \
    -M 5 \
    -C 3 \
    -O 3 \
    -A mean \
    --alpha 0.5 \
    --height 3 \
    --width 10 \
    --shrink \
    -o ${OUTPUT_DIR}/sashimi_novel_exon_zoomed

echo ""
echo "=== Sashimi plot generation completed: $(date) ==="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/sashimi_novel_exon.pdf"
echo "  - ${OUTPUT_DIR}/sashimi_novel_exon_zoomed.pdf"
