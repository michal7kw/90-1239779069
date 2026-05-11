#!/bin/bash
# Extract junction reads for novel SRRM3 exon analysis
# Separate analysis directory - does NOT modify existing pipeline

set -e

ANALYSIS_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/SRRM3_novel_exon"
BAM_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/results/02_aligned"
OUTPUT_DIR="${ANALYSIS_DIR}/results"
LOG_FILE="${ANALYSIS_DIR}/logs/01_extract_junctions.log"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${ANALYSIS_DIR}/logs

exec > >(tee -a ${LOG_FILE}) 2>&1
echo "=== Junction extraction started: $(date) ==="

# Target junction coordinates (1-based, using exact splice site positions)
# Novel exon: chr5:135,898,574-135,898,652
#
# From annotation:
# Upstream annotated exon ends at: 135,898,148
# Novel exon starts at: 135,898,574 and ends at: 135,898,652
# Downstream annotated exon starts at: 135,901,934
#
# Junction coordinates (intron boundaries, 1-based):
# Inclusion junction 1: 135,898,148 -> 135,898,574 (upstream end -> novel start)
# Inclusion junction 2: 135,898,652 -> 135,901,934 (novel end -> downstream start)
# Skipping junction: 135,898,148 -> 135,901,934 (upstream end -> downstream start)

echo "Target region: chr5:135897000-135903000"
echo "Novel exon: chr5:135898574-135898652"
echo ""

# Create Python script for parsing junctions
cat > ${ANALYSIS_DIR}/scripts/parse_junctions.py << 'PYEOF'
#!/usr/bin/env python3
import sys
import re

def parse_cigar(cigar, pos):
    """Parse CIGAR string and extract splice junctions"""
    junctions = []
    current_pos = pos

    # Parse CIGAR operations
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    for match in pattern.finditer(cigar):
        length = int(match.group(1))
        op = match.group(2)

        if op in ('M', 'D', '=', 'X'):  # Consume reference
            current_pos += length
        elif op == 'N':  # Splice junction
            # Junction: from (current_pos - 1) to (current_pos + length)
            # In 1-based coords: junction spans from current_pos to current_pos + length - 1
            junc_start = current_pos  # Last base of upstream exon (1-based)
            junc_end = current_pos + length  # First base of downstream exon (1-based)
            junctions.append((junc_start, junc_end))
            current_pos += length
        elif op in ('I', 'S', 'H', 'P'):  # Don't consume reference
            pass

    return junctions

# Read SAM input from stdin
for line in sys.stdin:
    if line.startswith('@'):
        continue
    fields = line.strip().split('\t')
    if len(fields) < 6:
        continue

    chrom = fields[2]
    pos = int(fields[3])  # 1-based position
    cigar = fields[5]

    if 'N' not in cigar:
        continue

    junctions = parse_cigar(cigar, pos)
    for junc_start, junc_end in junctions:
        print(f"{chrom}\t{junc_start}\t{junc_end}")
PYEOF

chmod +x ${ANALYSIS_DIR}/scripts/parse_junctions.py

# Initialize summary file
echo -e "sample\tgroup\tinc1_reads\tinc2_reads\tskip_reads\tnovel_coverage" > ${OUTPUT_DIR}/junction_counts_summary.csv

# Define sample groups
declare -A sample_groups
sample_groups[1]="Parental"
sample_groups[2]="Parental"
sample_groups[3]="Parental"
sample_groups[4]="Neg"
sample_groups[5]="Neg"
sample_groups[6]="Neg"
sample_groups[7]="Pos"
sample_groups[8]="Pos"
sample_groups[9]="Pos"
sample_groups[13]="KO"
sample_groups[14]="KO"
sample_groups[15]="KO"

# Process each sample
for sample in 1 2 3 4 5 6 7 8 9 13 14 15; do
    BAM="${BAM_DIR}/${sample}/${sample}_Aligned.sortedByCoord.out.bam"
    GROUP="${sample_groups[$sample]}"

    echo "Processing sample ${sample} (${GROUP})..."

    if [[ ! -f "${BAM}" ]]; then
        echo "  WARNING: BAM file not found: ${BAM}"
        continue
    fi

    # Extract all splice junctions in the region
    samtools view ${BAM} chr5:135897000-135903000 | \
        python3 ${ANALYSIS_DIR}/scripts/parse_junctions.py | \
        sort | uniq -c | \
        awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${OUTPUT_DIR}/sample_${sample}_junctions.txt

    # Count specific junctions with tolerance
    # Inclusion junction 1: 135898148 -> 135898574 (allow ±3bp tolerance)
    INC1=$(awk '$3 >= 135898145 && $3 <= 135898151 && $4 >= 135898571 && $4 <= 135898577 {sum += $1} END {print sum+0}' \
        ${OUTPUT_DIR}/sample_${sample}_junctions.txt)

    # Inclusion junction 2: 135898652 -> 135901934 (allow ±3bp tolerance)
    INC2=$(awk '$3 >= 135898649 && $3 <= 135898655 && $4 >= 135901931 && $4 <= 135901937 {sum += $1} END {print sum+0}' \
        ${OUTPUT_DIR}/sample_${sample}_junctions.txt)

    # Skipping junction: 135898148 -> 135901934 (allow ±3bp tolerance)
    SKIP=$(awk '$3 >= 135898145 && $3 <= 135898151 && $4 >= 135901931 && $4 <= 135901937 {sum += $1} END {print sum+0}' \
        ${OUTPUT_DIR}/sample_${sample}_junctions.txt)

    # Count reads overlapping the novel exon
    COVERAGE=$(samtools view -c ${BAM} chr5:135898574-135898652)

    echo "  Inc1 (upstream->novel): ${INC1}"
    echo "  Inc2 (novel->downstream): ${INC2}"
    echo "  Skip (upstream->downstream): ${SKIP}"
    echo "  Novel exon coverage: ${COVERAGE}"

    # Add to summary
    echo -e "${sample}\t${GROUP}\t${INC1}\t${INC2}\t${SKIP}\t${COVERAGE}" >> ${OUTPUT_DIR}/junction_counts_summary.csv
done

echo ""
echo "=== Junction extraction completed: $(date) ==="
echo "Output files:"
echo "  - ${OUTPUT_DIR}/junction_counts_summary.csv"
echo "  - ${OUTPUT_DIR}/sample_*_junctions.txt"
