#!/bin/bash
#SBATCH --job-name=8_igv_prep
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8_igv_prep.err"
#SBATCH --output="./logs/8_igv_prep.out"

# =============================================================================
# Prepare files for IGV visualization of differential splicing
# Creates:
#   1. BED files for significant splicing events (from rMATS)
#   2. BigWig coverage files (from BAM)
#   3. IGV session XML file
# =============================================================================

set -euo pipefail

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
BAM_DIR="${BASE_DIR}/results/02_aligned"
SPLICING_DIR="${BASE_DIR}/results/05_splicing"
OUTPUT_DIR="${BASE_DIR}/results/07_igv"
CHROM_SIZES="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/star/chrNameLength.txt"

echo "============================================"
echo "Preparing IGV visualization files"
echo "Start time: $(date)"
echo "============================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rmats

mkdir -p ${OUTPUT_DIR}/{bed,bigwig}

# =============================================================================
# 1. Create BED files from rMATS significant splicing events
# =============================================================================
echo ""
echo "=== Creating BED files for significant splicing events ==="

python3 << 'PYEOF'
import os
import pandas as pd
from pathlib import Path

base_dir = os.environ.get('BASE_DIR', '/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069')
splicing_dir = Path(base_dir) / 'results/05_splicing'
output_dir = Path(base_dir) / 'results/07_igv/bed'

comparisons = ['Neg_vs_Parental', 'Pos_vs_Parental', 'KO_vs_Parental',
               'Pos_vs_Neg', 'KO_vs_Neg', 'KO_vs_Pos']
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

# Color scheme for different event types (RGB)
colors = {
    'SE': '255,0,0',      # Red - Skipped Exon
    'A5SS': '0,255,0',    # Green - Alt 5' splice site
    'A3SS': '0,0,255',    # Blue - Alt 3' splice site
    'MXE': '255,165,0',   # Orange - Mutually exclusive exons
    'RI': '128,0,128'     # Purple - Retained intron
}

for comp in comparisons:
    all_events = []

    for event in event_types:
        jc_file = splicing_dir / comp / f'{event}.MATS.JC.txt'
        if not jc_file.exists():
            continue

        try:
            df = pd.read_csv(jc_file, sep='\t')
            # Filter significant events: FDR < 0.05 and |dPSI| > 0.1
            sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)].copy()

            if len(sig) == 0:
                continue

            for _, row in sig.iterrows():
                chrom = row['chr']
                # Get coordinates based on event type
                if event == 'SE':
                    start = int(row['exonStart_0base'])
                    end = int(row['exonEnd'])
                elif event == 'RI':
                    start = int(row['riExonStart_0base'])
                    end = int(row['riExonEnd'])
                elif event in ['A5SS', 'A3SS']:
                    start = min(int(row['longExonStart_0base']), int(row['shortES']))
                    end = max(int(row['longExonEnd']), int(row['shortEE']))
                elif event == 'MXE':
                    start = min(int(row['1stExonStart_0base']), int(row['2ndExonStart_0base']))
                    end = max(int(row['1stExonEnd']), int(row['2ndExonEnd']))

                gene = row.get('geneSymbol', row.get('GeneID', 'Unknown'))
                if isinstance(gene, str):
                    gene = gene.replace('"', '')

                dpsi = row['IncLevelDifference']
                fdr = row['FDR']

                # BED format: chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb
                name = f"{gene}_{event}_dPSI={dpsi:.2f}"
                score = min(1000, int(abs(dpsi) * 1000))  # Scale dPSI to score
                strand = row.get('strand', '.')

                all_events.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'score': score,
                    'strand': strand,
                    'thickStart': start,
                    'thickEnd': end,
                    'itemRgb': colors[event],
                    'event_type': event,
                    'dpsi': dpsi,
                    'fdr': fdr
                })
        except Exception as e:
            print(f"  Error processing {jc_file}: {e}")

    if all_events:
        # Write combined BED file for this comparison
        bed_df = pd.DataFrame(all_events)
        bed_file = output_dir / f'{comp}_significant_splicing.bed'

        # Write BED header
        with open(bed_file, 'w') as f:
            f.write(f'track name="{comp}_splicing" description="Significant splicing events (FDR<0.05, |dPSI|>0.1)" itemRgb="On"\n')

        # Append data
        bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'itemRgb']].to_csv(
            bed_file, sep='\t', header=False, index=False, mode='a'
        )

        print(f"  {comp}: {len(all_events)} significant events written to BED")

        # Also write separate BED files per event type
        for event in event_types:
            event_df = bed_df[bed_df['event_type'] == event]
            if len(event_df) > 0:
                event_file = output_dir / f'{comp}_{event}.bed'
                with open(event_file, 'w') as f:
                    f.write(f'track name="{comp}_{event}" description="{event} events" color="{colors[event].replace(",", " ")}"\n')
                event_df[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
                    event_file, sep='\t', header=False, index=False, mode='a'
                )

print("\nBED files created successfully!")
PYEOF

export BASE_DIR

# =============================================================================
# 2. Create BigWig files from BAM (for smoother coverage visualization)
# =============================================================================
echo ""
echo "=== Creating BigWig coverage files ==="

# Check if chromosome sizes file exists, if not create from BAM header
if [ ! -f "${CHROM_SIZES}" ]; then
    echo "Creating chromosome sizes file..."
    samtools view -H ${BAM_DIR}/1/1_Aligned.sortedByCoord.out.bam | \
        grep "^@SQ" | cut -f2,3 | sed 's/SN://;s/LN:/\t/' > ${OUTPUT_DIR}/chrom.sizes
    CHROM_SIZES="${OUTPUT_DIR}/chrom.sizes"
fi

# Create BigWig for each sample
declare -A SAMPLES=(
    ["1"]="Parental_1"
    ["2"]="Parental_2"
    ["3"]="Parental_3"
    ["4"]="Neg_1"
    ["5"]="Neg_2"
    ["6"]="Neg_3"
    ["7"]="Pos_1"
    ["8"]="Pos_2"
    ["9"]="Pos_3"
    ["13"]="KO_1"
    ["14"]="KO_2"
    ["15"]="KO_3"
)

for sample_id in "${!SAMPLES[@]}"; do
    sample_name="${SAMPLES[$sample_id]}"
    bam_file="${BAM_DIR}/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"
    bigwig_file="${OUTPUT_DIR}/bigwig/${sample_name}.bw"

    if [ -f "${bam_file}" ] && [ ! -f "${bigwig_file}" ]; then
        echo "  Creating BigWig for ${sample_name}..."
        bamCoverage -b ${bam_file} -o ${bigwig_file} \
            --binSize 10 \
            --normalizeUsing RPKM \
            --numberOfProcessors 8 \
            2>/dev/null || echo "    Warning: bamCoverage failed for ${sample_name}"
    else
        echo "  Skipping ${sample_name} (already exists or BAM not found)"
    fi
done

# =============================================================================
# 3. Create IGV session file
# =============================================================================
echo ""
echo "=== Creating IGV session file ==="

cat > ${OUTPUT_DIR}/splicing_analysis.xml << 'XMLEOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">
    <Resources>
XMLEOF

# Add BigWig resources
for sample_name in Parental_1 Parental_2 Parental_3 Neg_1 Neg_2 Neg_3 Pos_1 Pos_2 Pos_3 KO_1 KO_2 KO_3; do
    echo "        <Resource path=\"bigwig/${sample_name}.bw\"/>" >> ${OUTPUT_DIR}/splicing_analysis.xml
done

# Add BAM resources
for sample_id in 1 2 3 4 5 6 7 8 9 13 14 15; do
    echo "        <Resource path=\"${BAM_DIR}/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam\"/>" >> ${OUTPUT_DIR}/splicing_analysis.xml
done

# Add BED resources
for comp in Neg_vs_Parental Pos_vs_Parental KO_vs_Parental; do
    echo "        <Resource path=\"bed/${comp}_significant_splicing.bed\"/>" >> ${OUTPUT_DIR}/splicing_analysis.xml
done

cat >> ${OUTPUT_DIR}/splicing_analysis.xml << 'XMLEOF'
    </Resources>
    <Panel name="DataPanel" height="600">
        <!-- Coverage tracks will be loaded here -->
    </Panel>
    <Panel name="FeaturePanel" height="200">
        <!-- BED tracks will be loaded here -->
    </Panel>
    <PanelLayout dividerFractions="0.75"/>
</Session>
XMLEOF

# =============================================================================
# 4. Create a regions of interest file for top splicing events
# =============================================================================
echo ""
echo "=== Creating regions of interest file ==="

python3 << 'PYEOF2'
import os
import pandas as pd
from pathlib import Path

base_dir = os.environ.get('BASE_DIR', '/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069')
splicing_dir = Path(base_dir) / 'results/05_splicing'
output_dir = Path(base_dir) / 'results/07_igv'

# Collect top events across all comparisons
all_top_events = []

comparisons = ['Neg_vs_Parental', 'Pos_vs_Parental', 'KO_vs_Parental']
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

for comp in comparisons:
    for event in event_types:
        jc_file = splicing_dir / comp / f'{event}.MATS.JC.txt'
        if not jc_file.exists():
            continue

        try:
            df = pd.read_csv(jc_file, sep='\t')
            sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)]

            # Get top 3 by absolute dPSI
            top = sig.nlargest(3, 'IncLevelDifference', keep='first')

            for _, row in top.iterrows():
                chrom = row['chr']
                if event == 'SE':
                    start = int(row['upstreamES']) - 500
                    end = int(row['downstreamEE']) + 500
                elif event == 'RI':
                    start = int(row['riExonStart_0base']) - 500
                    end = int(row['riExonEnd']) + 500
                else:
                    start = int(row.get('longExonStart_0base', row.get('1stExonStart_0base', 0))) - 500
                    end = int(row.get('longExonEnd', row.get('2ndExonEnd', 0))) + 500

                gene = str(row.get('geneSymbol', row.get('GeneID', 'Unknown'))).replace('"', '')
                dpsi = row['IncLevelDifference']

                all_top_events.append({
                    'region': f"{chrom}:{max(1,start)}-{end}",
                    'description': f"{gene}_{event}_{comp}_dPSI={dpsi:.2f}"
                })
        except Exception as e:
            pass

# Write regions file (IGV batch format)
with open(output_dir / 'top_splicing_regions.txt', 'w') as f:
    f.write("# Top differential splicing regions for IGV\n")
    f.write("# Use: View > Go to regions in file\n")
    for evt in all_top_events[:50]:  # Top 50 regions
        f.write(f"{evt['region']}\t{evt['description']}\n")

# Write as BED for region navigation
with open(output_dir / 'bed/top_regions.bed', 'w') as f:
    f.write('track name="Top_Splicing_Regions" description="Top differential splicing events"\n')
    for evt in all_top_events[:50]:
        region = evt['region']
        chrom, coords = region.split(':')
        start, end = coords.split('-')
        f.write(f"{chrom}\t{start}\t{end}\t{evt['description']}\n")

print(f"Created regions file with {min(50, len(all_top_events))} top events")
PYEOF2

echo ""
echo "============================================"
echo "IGV preparation complete!"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Files created:"
echo "  - bed/              : BED files for significant splicing events"
echo "  - bigwig/           : Coverage tracks for each sample"
echo "  - splicing_analysis.xml : IGV session file"
echo "  - top_splicing_regions.txt : Quick navigation regions"
echo ""
echo "To use in IGV:"
echo "  1. Open IGV and set genome to mm39 (GRCm39)"
echo "  2. File > Load Session > splicing_analysis.xml"
echo "  3. Or manually load BED/BigWig files"
echo "  4. Right-click on BAM track > Sashimi Plot"
echo ""
echo "End time: $(date)"
echo "============================================"
