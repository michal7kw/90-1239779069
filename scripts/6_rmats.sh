#!/bin/bash
#SBATCH --job-name=6_rmats
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6_rmats.err"
#SBATCH --output="./logs/6_rmats.out"

# =============================================================================
# rMATS Differential Splicing Analysis
# Project: 90-1239779069
# Reference: Parental
# Comparisons: All pairwise
# =============================================================================

set -euo pipefail

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
BAM_DIR="${BASE_DIR}/results/02_aligned"
OUTPUT_DIR="${BASE_DIR}/results/05_splicing"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"
TMP_DIR="${OUTPUT_DIR}/tmp"

echo "============================================"
echo "rMATS Differential Splicing Analysis"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rmats

# Create directories
mkdir -p ${OUTPUT_DIR}/{config,Neg_vs_Parental,Pos_vs_Parental,KO_vs_Parental,Pos_vs_Neg,KO_vs_Neg,KO_vs_Pos}
mkdir -p ${TMP_DIR}

# Create BAM file lists for each group
echo "Creating BAM file lists..."

# Parental (samples 1, 2, 3)
echo "${BAM_DIR}/1/1_Aligned.sortedByCoord.out.bam,${BAM_DIR}/2/2_Aligned.sortedByCoord.out.bam,${BAM_DIR}/3/3_Aligned.sortedByCoord.out.bam" > ${OUTPUT_DIR}/config/Parental.txt

# Neg (samples 4, 5, 6)
echo "${BAM_DIR}/4/4_Aligned.sortedByCoord.out.bam,${BAM_DIR}/5/5_Aligned.sortedByCoord.out.bam,${BAM_DIR}/6/6_Aligned.sortedByCoord.out.bam" > ${OUTPUT_DIR}/config/Neg.txt

# Pos (samples 7, 8, 9)
echo "${BAM_DIR}/7/7_Aligned.sortedByCoord.out.bam,${BAM_DIR}/8/8_Aligned.sortedByCoord.out.bam,${BAM_DIR}/9/9_Aligned.sortedByCoord.out.bam" > ${OUTPUT_DIR}/config/Pos.txt

# KO (samples 13, 14, 15)
echo "${BAM_DIR}/13/13_Aligned.sortedByCoord.out.bam,${BAM_DIR}/14/14_Aligned.sortedByCoord.out.bam,${BAM_DIR}/15/15_Aligned.sortedByCoord.out.bam" > ${OUTPUT_DIR}/config/KO.txt

# Function to run rMATS comparison
run_rmats() {
    local GROUP1=$1
    local GROUP2=$2
    local COMPARISON="${GROUP1}_vs_${GROUP2}"
    local OUT="${OUTPUT_DIR}/${COMPARISON}"

    echo ""
    echo "============================================"
    echo "Running: ${COMPARISON}"
    echo "============================================"

    rmats.py \
        --b1 ${OUTPUT_DIR}/config/${GROUP1}.txt \
        --b2 ${OUTPUT_DIR}/config/${GROUP2}.txt \
        --gtf ${GTF_FILE} \
        -t paired \
        --readLength 150 \
        --nthread 16 \
        --od ${OUT} \
        --tmp ${TMP_DIR}/${COMPARISON} \
        --statoff

    # Run statistical analysis separately (allows for more memory)
    rmats.py \
        --b1 ${OUTPUT_DIR}/config/${GROUP1}.txt \
        --b2 ${OUTPUT_DIR}/config/${GROUP2}.txt \
        --gtf ${GTF_FILE} \
        -t paired \
        --readLength 150 \
        --nthread 16 \
        --od ${OUT} \
        --tmp ${TMP_DIR}/${COMPARISON} \
        --task stat

    echo "Completed: ${COMPARISON}"
}

# Run all pairwise comparisons (with Parental as reference)
echo ""
echo "Running pairwise comparisons..."

# Comparisons against Parental (reference)
run_rmats "Neg" "Parental"
run_rmats "Pos" "Parental"
run_rmats "KO" "Parental"

# Other pairwise comparisons
run_rmats "Pos" "Neg"
run_rmats "KO" "Neg"
run_rmats "KO" "Pos"

# Generate summary
echo ""
echo "============================================"
echo "Generating summary..."
echo "============================================"

# Create summary of significant splicing events
cat > ${OUTPUT_DIR}/summarize_splicing.py << 'PYEOF'
#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path

output_dir = os.environ.get('OUTPUT_DIR', '.')
comparisons = ['Neg_vs_Parental', 'Pos_vs_Parental', 'KO_vs_Parental',
               'Pos_vs_Neg', 'KO_vs_Neg', 'KO_vs_Pos']
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

summary = []
for comp in comparisons:
    for event in event_types:
        jc_file = Path(output_dir) / comp / f'{event}.MATS.JC.txt'
        if jc_file.exists():
            df = pd.read_csv(jc_file, sep='\t')
            sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)]
            summary.append({
                'Comparison': comp,
                'Event_Type': event,
                'Total_Events': len(df),
                'Significant_FDR005_dPSI01': len(sig)
            })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(Path(output_dir) / 'splicing_summary.csv', index=False)
print(summary_df.to_string(index=False))
PYEOF

export OUTPUT_DIR
python3 ${OUTPUT_DIR}/summarize_splicing.py

# Cleanup temp files
rm -rf ${TMP_DIR}

echo ""
echo "============================================"
echo "rMATS analysis complete!"
echo "Output directory: ${OUTPUT_DIR}"
echo "End time: $(date)"
echo "============================================"
