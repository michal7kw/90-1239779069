#!/bin/bash
#SBATCH --job-name=6h_microexon_extended
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6h_microexon_extended.err"
#SBATCH --output="./logs/6h_microexon_extended.out"

# =============================================================================
# Extended Microexon Analysis Pipeline
# Project: 90-1239779069
# Scripts: 6e, 6f, 6g - Functional annotation, expression correlation, comparative
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SCRIPTS_DIR="${BASE_DIR}/scripts"
OUTPUT_DIR="${BASE_DIR}/results/09_microexon_extended"

echo "============================================"
echo "Extended Microexon Analysis Pipeline"
echo "Start time: $(date)"
echo "============================================"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/splicing_analysis

# Create directories
mkdir -p ${BASE_DIR}/logs
mkdir -p ${OUTPUT_DIR}/functional_annotation
mkdir -p ${OUTPUT_DIR}/expression_correlation
mkdir -p ${OUTPUT_DIR}/comparative_analysis
mkdir -p ${OUTPUT_DIR}/summary

# =============================================================================
# Step 1: Functional Annotation (5'UTR, CDS, 3'UTR classification)
# =============================================================================
echo ""
echo "============================================"
echo "Step 1: Functional Annotation"
echo "============================================"

Rscript ${SCRIPTS_DIR}/6e_microexon_functional_annotation.R 2>&1 | tee ${BASE_DIR}/logs/6e_functional_annotation.log

if [ $? -eq 0 ]; then
    echo "Step 1 completed successfully"
else
    echo "Step 1 failed" >&2
    exit 1
fi

# =============================================================================
# Step 2: Expression Correlation (dPSI vs log2FC)
# =============================================================================
echo ""
echo "============================================"
echo "Step 2: Expression Correlation"
echo "============================================"

Rscript ${SCRIPTS_DIR}/6f_microexon_expression_correlation.R 2>&1 | tee ${BASE_DIR}/logs/6f_expression_correlation.log

if [ $? -eq 0 ]; then
    echo "Step 2 completed successfully"
else
    echo "Step 2 failed" >&2
    exit 1
fi

# =============================================================================
# Step 3: Comparative Analysis (Common/unique with directionality)
# =============================================================================
echo ""
echo "============================================"
echo "Step 3: Comparative Analysis"
echo "============================================"

Rscript ${SCRIPTS_DIR}/6g_microexon_comparative_analysis.R 2>&1 | tee ${BASE_DIR}/logs/6g_comparative_analysis.log

if [ $? -eq 0 ]; then
    echo "Step 3 completed successfully"
else
    echo "Step 3 failed" >&2
    exit 1
fi

# =============================================================================
# Generate Summary Report
# =============================================================================
echo ""
echo "============================================"
echo "Generating Summary"
echo "============================================"

# Create summary file
SUMMARY_FILE="${OUTPUT_DIR}/summary/extended_microexon_summary.txt"

cat > ${SUMMARY_FILE} << 'EOF'
==============================================================================
EXTENDED MICROEXON ANALYSIS SUMMARY
Project: 90-1239779069
==============================================================================

ANALYSIS OVERVIEW
-----------------
This analysis extends the microexon characterization with:
1. Functional annotation (5'UTR, CDS, 3'UTR classification)
2. Expression-splicing correlation analysis
3. Comprehensive comparative analysis across SRRM3 perturbations

EXPERIMENTAL GROUPS
-------------------
- Parental (samples 1,2,3): Normal SRRM3 - Baseline
- Neg (samples 4,5,6): Knockdown - Partial loss of SRRM3
- Pos (samples 7,8,9): Overexpression - Gain of SRRM3
- KO (samples 13,14,15): Knockout - Complete loss of SRRM3

BIOLOGICAL EXPECTATIONS
-----------------------
- Neg & KO (loss of function): Decreased microexon inclusion (negative dPSI)
- Pos (gain of function): Increased microexon inclusion (positive dPSI)
- Neg vs KO: Same direction, KO effect should be stronger
- Neg/KO vs Pos: Opposite directions

OUTPUT FILES
------------
results/09_microexon_extended/
├── functional_annotation/
│   ├── microexon_genomic_annotation.csv  - Full annotation table
│   ├── position_summary_stats.csv        - Statistics by position
│   ├── biotype_summary_stats.csv         - Statistics by gene biotype
│   └── *.pdf                             - Distribution plots
├── expression_correlation/
│   ├── *_correlation.csv                 - Correlation statistics per comparison
│   ├── correlation_summary.csv           - Combined correlations
│   ├── correlation_heatmap.pdf           - Heatmap visualization
│   └── *_scatter_*.pdf                   - Scatter plots
├── comparative_analysis/
│   ├── dPSI_comparison_matrix.csv        - Wide-format dPSI matrix
│   ├── common_microexons_all.csv         - Events significant in all comparisons
│   ├── unique_*.csv                      - Comparison-specific events
│   ├── directionality_patterns.csv       - Direction pattern counts
│   ├── concordance_analysis.csv          - Pairwise concordance
│   └── *.pdf                             - UpSet, Venn, flow diagrams
└── summary/
    └── extended_microexon_summary.txt    - This file

THRESHOLDS USED
---------------
- Splicing significance: FDR < 0.05, |dPSI| > 0.1
- Expression significance: padj < 0.05, |log2FC| > 1
- Microexon definition: 0-30bp
- Small exon definition: 30-50bp

==============================================================================
EOF

echo "Summary written to: ${SUMMARY_FILE}"

# =============================================================================
# Pipeline Complete
# =============================================================================
echo ""
echo "============================================"
echo "Extended Microexon Analysis Complete!"
echo "============================================"
echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Sub-directories:"
echo "  - functional_annotation/: Genomic position and biotype analysis"
echo "  - expression_correlation/: dPSI vs log2FC correlation"
echo "  - comparative_analysis/: Common/unique events and concordance"
echo "  - summary/: Summary statistics"
echo ""
echo "End time: $(date)"
echo "============================================"
