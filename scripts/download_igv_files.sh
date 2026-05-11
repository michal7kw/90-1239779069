#!/bin/bash
# =============================================================================
# IGV Visualization Files Download Helper
# =============================================================================
# Run this script on your LOCAL machine to download files from the server
#
# Usage:
#   ./download_igv_files.sh [tier] [server] [local_dest]
#
# Tiers:
#   1 or minimal   - Pre-rendered sashimi only (~36 MB)
#   2 or coverage  - BigWig + BED files (~420 MB)
#   3 or sashimi   - Full subset with BAM for sashimi (~10 GB)
#   4 or full      - Complete visualization package
#
# Example:
#   ./download_igv_files.sh coverage user@server.com ./igv_data
# =============================================================================

set -e

# Configuration - UPDATE THESE
SERVER="${2:-user@your-server.com}"
REMOTE_BASE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/results"
LOCAL_DEST="${3:-./igv_splicing_data}"

TIER="${1:-2}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=============================================="
echo "IGV Visualization Files Download"
echo "=============================================="
echo ""

# Create local destination
mkdir -p "${LOCAL_DEST}"
cd "${LOCAL_DEST}"

case "$TIER" in
    1|minimal)
        echo -e "${GREEN}Tier 1: Minimal - Pre-rendered Sashimi Plots${NC}"
        echo "Downloading ~36 MB..."
        echo ""

        mkdir -p sashimi_plots
        scp -r "${SERVER}:${REMOTE_BASE}/06_sashimi/top_splicing/" sashimi_plots/

        echo ""
        echo -e "${GREEN}Download complete!${NC}"
        echo "Files in: ${LOCAL_DEST}/sashimi_plots/"
        echo "Open PDF/PNG files directly - no IGV needed"
        ;;

    2|coverage)
        echo -e "${GREEN}Tier 2: Coverage - BigWig + BED files${NC}"
        echo "Downloading ~420 MB..."
        echo ""

        mkdir -p coverage annotations navigation

        echo "Downloading BigWig files..."
        scp "${SERVER}:${REMOTE_BASE}/08_igv_subset/bigwig/*.bw" coverage/

        echo "Downloading BED annotations..."
        scp "${SERVER}:${REMOTE_BASE}/07_igv/bed/*.bed" annotations/

        echo "Downloading navigation file..."
        scp "${SERVER}:${REMOTE_BASE}/07_igv/top_splicing_regions.txt" navigation/

        # Create simple session file
        cat > session_coverage.xml << 'EOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <Resource path="coverage/Parental_1_subset.bw"/>
        <Resource path="coverage/Parental_2_subset.bw"/>
        <Resource path="coverage/Parental_3_subset.bw"/>
        <Resource path="coverage/Neg_1_subset.bw"/>
        <Resource path="coverage/Neg_2_subset.bw"/>
        <Resource path="coverage/Neg_3_subset.bw"/>
        <Resource path="coverage/Pos_1_subset.bw"/>
        <Resource path="coverage/Pos_2_subset.bw"/>
        <Resource path="coverage/Pos_3_subset.bw"/>
        <Resource path="coverage/KO_1_subset.bw"/>
        <Resource path="coverage/KO_2_subset.bw"/>
        <Resource path="coverage/KO_3_subset.bw"/>
        <Resource path="annotations/KO_vs_Parental_significant_splicing.bed"/>
        <Resource path="annotations/Neg_vs_Parental_significant_splicing.bed"/>
        <Resource path="annotations/Pos_vs_Parental_significant_splicing.bed"/>
        <Resource path="annotations/top_regions.bed"/>
    </Resources>
</Session>
EOF

        echo ""
        echo -e "${GREEN}Download complete!${NC}"
        echo "Files in: ${LOCAL_DEST}/"
        echo ""
        echo "To use in IGV:"
        echo "  1. Open IGV, set genome to mm39"
        echo "  2. File → Load Session → session_coverage.xml"
        echo "  3. Regions → Region Navigator → Load navigation/top_splicing_regions.txt"
        ;;

    3|sashimi)
        echo -e "${GREEN}Tier 3: Full Interactive - BigWig + BAM for Sashimi${NC}"
        echo "Downloading ~10 GB (this will take a while)..."
        echo ""

        mkdir -p coverage alignments annotations navigation

        echo "Downloading BigWig files..."
        scp "${SERVER}:${REMOTE_BASE}/08_igv_subset/bigwig/*.bw" coverage/

        echo "Downloading BAM files (large)..."
        scp "${SERVER}:${REMOTE_BASE}/08_igv_subset/bam/*.bam" alignments/
        scp "${SERVER}:${REMOTE_BASE}/08_igv_subset/bam/*.bam.bai" alignments/

        echo "Downloading BED annotations..."
        scp "${SERVER}:${REMOTE_BASE}/07_igv/bed/*.bed" annotations/

        echo "Downloading navigation file..."
        scp "${SERVER}:${REMOTE_BASE}/07_igv/top_splicing_regions.txt" navigation/

        # Create full session file
        cat > session_full.xml << 'EOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm39" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <!-- Coverage tracks -->
        <Resource path="coverage/Parental_1_subset.bw"/>
        <Resource path="coverage/Parental_2_subset.bw"/>
        <Resource path="coverage/Parental_3_subset.bw"/>
        <Resource path="coverage/Neg_1_subset.bw"/>
        <Resource path="coverage/Neg_2_subset.bw"/>
        <Resource path="coverage/Neg_3_subset.bw"/>
        <Resource path="coverage/Pos_1_subset.bw"/>
        <Resource path="coverage/Pos_2_subset.bw"/>
        <Resource path="coverage/Pos_3_subset.bw"/>
        <Resource path="coverage/KO_1_subset.bw"/>
        <Resource path="coverage/KO_2_subset.bw"/>
        <Resource path="coverage/KO_3_subset.bw"/>
        <!-- BAM files for sashimi -->
        <Resource path="alignments/Parental_1_subset.bam"/>
        <Resource path="alignments/Parental_2_subset.bam"/>
        <Resource path="alignments/Parental_3_subset.bam"/>
        <Resource path="alignments/Neg_1_subset.bam"/>
        <Resource path="alignments/Neg_2_subset.bam"/>
        <Resource path="alignments/Neg_3_subset.bam"/>
        <Resource path="alignments/Pos_1_subset.bam"/>
        <Resource path="alignments/Pos_2_subset.bam"/>
        <Resource path="alignments/Pos_3_subset.bam"/>
        <Resource path="alignments/KO_1_subset.bam"/>
        <Resource path="alignments/KO_2_subset.bam"/>
        <Resource path="alignments/KO_3_subset.bam"/>
        <!-- Annotations -->
        <Resource path="annotations/KO_vs_Parental_significant_splicing.bed"/>
        <Resource path="annotations/Neg_vs_Parental_significant_splicing.bed"/>
        <Resource path="annotations/Pos_vs_Parental_significant_splicing.bed"/>
        <Resource path="annotations/top_regions.bed"/>
    </Resources>
</Session>
EOF

        echo ""
        echo -e "${GREEN}Download complete!${NC}"
        echo "Files in: ${LOCAL_DEST}/"
        echo ""
        echo "To use in IGV:"
        echo "  1. Open IGV, set genome to mm39"
        echo "  2. File → Load Session → session_full.xml"
        echo "  3. For sashimi: Right-click BAM track → Sashimi Plot"
        ;;

    4|full)
        echo -e "${YELLOW}Tier 4: Complete Package${NC}"
        echo "Downloading full visualization package..."
        echo ""

        scp "${SERVER}:${REMOTE_BASE}/08_igv_subset/igv_subset_package.tar.gz" .

        echo "Extracting..."
        tar -xzf igv_subset_package.tar.gz

        echo ""
        echo -e "${GREEN}Download complete!${NC}"
        echo "Files in: ${LOCAL_DEST}/"
        echo ""
        echo "To use: Open splicing_subset_session.xml in IGV (genome: mm39)"
        ;;

    *)
        echo -e "${RED}Unknown tier: $TIER${NC}"
        echo ""
        echo "Available tiers:"
        echo "  1 or minimal  - Pre-rendered sashimi only (~36 MB)"
        echo "  2 or coverage - BigWig + BED files (~420 MB)"
        echo "  3 or sashimi  - Full subset with BAM (~10 GB)"
        echo "  4 or full     - Complete package from tarball"
        exit 1
        ;;
esac

echo ""
echo "=============================================="
