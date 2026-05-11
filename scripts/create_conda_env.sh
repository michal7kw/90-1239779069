#!/bin/bash
#SBATCH --job-name=create_conda
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=./logs/create_conda.out
#SBATCH --error=./logs/create_conda.err

# Script to create the conda environment for splicing analysis
# Usage: sbatch create_conda_env.sh

echo "Creating 'splicing_analysis' conda environment using mamba..."

set -e

# Create environment with all required packages
# r-base is included to ensure R is installed
# Channels: conda-forge, bioconda
# Packages based on script dependencies:
# - Core R: dplyr, ggplot2, tidyr, grid (base), etc.
# - Visualization: pheatmap, RColorBrewer, ggrepel, patchwork, UpSetR, VennDiagram, ggalluvial
# - Bioconductor: AnnotationDbi, clusterProfiler, DOSE, enrichplot, GenomicFeatures, GenomicRanges, org.Mm.eg.db, rtracklayer

mamba create -y -n splicing_analysis \
    -c conda-forge -c bioconda \
    r-base \
    r-dplyr \
    r-ggplot2 \
    r-ggrepel \
    r-patchwork \
    r-pheatmap \
    r-rcolorbrewer \
    r-tidyr \
    r-upsetr \
    r-venndiagram \
    r-ggalluvial \
    bioconductor-annotationdbi \
    bioconductor-clusterprofiler \
    bioconductor-dose \
    bioconductor-enrichplot \
    bioconductor-genomicfeatures \
    bioconductor-genomicranges \
    bioconductor-org.mm.eg.db \
    bioconductor-rtracklayer

echo "Environment 'splicing_analysis' created."
echo "Activate it with: conda activate splicing_analysis"
