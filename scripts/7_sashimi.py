#!/usr/bin/env python3
"""
Generate sashimi plots for top DE genes, top splicing events, and user-specified genes.

Features:
  - Aggregated replicates per group (-A mean)
  - Strict junction filtering (-M 100)
  - Dual output: PDF + PNG
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
import re

# Configuration from environment
BASE_DIR = os.environ.get('BASE_DIR')
OUTPUT_DIR = os.environ.get('OUTPUT_DIR')
DESEQ_DIR = os.environ.get('DESEQ_DIR')
SPLICING_DIR = os.environ.get('SPLICING_DIR')
GTF_FILE = os.environ.get('GTF_FILE')
USER_GENES_FILE = os.environ.get('USER_GENES')
BAM_CONFIG = f"{OUTPUT_DIR}/bam_config.tsv"

# Plot settings from environment
MIN_JUNCTION_COV = os.environ.get('MIN_JUNCTION_COV', '100')
PLOT_HEIGHT = os.environ.get('PLOT_HEIGHT', '5')
PLOT_WIDTH = os.environ.get('PLOT_WIDTH', '18')
PNG_RESOLUTION = os.environ.get('PNG_RESOLUTION', '300')

def parse_gtf_for_gene(gtf_file, gene_name):
    """Extract genomic coordinates for a gene from GTF file."""
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue
            attrs = fields[8]
            gene_match = re.search(r'gene_name "([^"]+)"', attrs)
            if gene_match and gene_match.group(1) == gene_name:
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                # Add flanking region
                start = max(1, start - 1000)
                end = end + 1000
                return f"{chrom}:{start}-{end}", strand
    return None, None

def parse_gtf_for_gene_id(gtf_file, gene_id):
    """Extract genomic coordinates for a gene from GTF file using gene_id."""
    gene_id_clean = gene_id.split('.')[0]
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue
            attrs = fields[8]
            id_match = re.search(r'gene_id "([^"]+)"', attrs)
            name_match = re.search(r'gene_name "([^"]+)"', attrs)
            if id_match:
                gtf_gene_id = id_match.group(1).split('.')[0]
                if gtf_gene_id == gene_id_clean:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    gene_name = name_match.group(1) if name_match else gene_id
                    start = max(1, start - 1000)
                    end = end + 1000
                    return f"{chrom}:{start}-{end}", gene_name
    return None, None

def run_ggsashimi(region, gene_name, output_subdir):
    """Run ggsashimi for a specific region - generates both PDF and PNG."""
    output_prefix = f"{OUTPUT_DIR}/{output_subdir}/{gene_name}"

    # Base command arguments (shared between PDF and PNG)
    base_args = [
        "ggsashimi.py",
        "-b", BAM_CONFIG,
        "-c", region,
        "-g", GTF_FILE,
        "-M", MIN_JUNCTION_COV,      # Strict filtering
        "-C", "3",                    # Color by group
        "-O", "3",                    # Overlay mode
        "-A", "mean",                 # Aggregate replicates by mean
        "--alpha", "0.5",
        "--height", PLOT_HEIGHT,
        "--width", PLOT_WIDTH,
        "--base-size", "14",
        "--ann-height", "3",
        "--shrink"
    ]

    print(f"  Generating: {gene_name} ({region})")
    success = False

    # Generate PDF
    try:
        cmd_pdf = base_args + ["-o", output_prefix, "-F", "pdf"]
        result = subprocess.run(cmd_pdf, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            success = True
        else:
            print(f"    Warning (PDF): {result.stderr[:200]}")
    except subprocess.TimeoutExpired:
        print(f"    Timeout for {gene_name} (PDF)")
    except Exception as e:
        print(f"    Error (PDF): {str(e)[:100]}")

    # Generate PNG (high resolution)
    try:
        cmd_png = base_args + ["-o", f"{output_prefix}_hires", "-F", "png", "-R", PNG_RESOLUTION]
        result = subprocess.run(cmd_png, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"    Warning (PNG): {result.stderr[:200]}")
    except subprocess.TimeoutExpired:
        print(f"    Timeout for {gene_name} (PNG)")
    except Exception as e:
        print(f"    Error (PNG): {str(e)[:100]}")

    return success

# Get top DE genes
print("\n=== Processing Top DE Genes ===")
de_genes_processed = set()
comparisons = ['Neg_vs_Parental', 'Pos_vs_Parental', 'KO_vs_Parental',
               'Pos_vs_Neg', 'KO_vs_Neg', 'KO_vs_Pos']

for comp in comparisons:
    sig_file = Path(DESEQ_DIR) / f"{comp}_significant.csv"
    if sig_file.exists():
        df = pd.read_csv(sig_file)
        # Get top 5 genes by adjusted p-value for each comparison
        top_genes = df.nsmallest(5, 'padj')['gene_id'].tolist()
        for gene_id in top_genes:
            if gene_id in de_genes_processed:
                continue
            region, gene_name = parse_gtf_for_gene_id(GTF_FILE, gene_id)
            if region:
                success = run_ggsashimi(region, gene_name, "top_DE")
                if success:
                    de_genes_processed.add(gene_id)
            if len(de_genes_processed) >= 20:
                break
    if len(de_genes_processed) >= 20:
        break

print(f"Generated {len(de_genes_processed)} DE gene plots")

# Get top splicing events
print("\n=== Processing Top Splicing Events ===")
splicing_genes_processed = set()
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

for comp in comparisons:
    for event in event_types:
        jc_file = Path(SPLICING_DIR) / comp / f'{event}.MATS.JC.txt'
        if jc_file.exists():
            try:
                df = pd.read_csv(jc_file, sep='\t')
                # Filter for significant events
                sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)]
                if len(sig) == 0:
                    continue
                # Get top 3 events
                top_events = sig.nsmallest(3, 'FDR')
                for _, row in top_events.iterrows():
                    gene_name = row['GeneID'] if 'GeneID' in row else row.get('geneSymbol', 'Unknown')
                    if gene_name in splicing_genes_processed:
                        continue
                    chrom = row['chr']
                    start = min(row.get('exonStart_0base', row.get('upstreamES', 0)),
                               row.get('upstreamES', row.get('exonStart_0base', 0)))
                    end = max(row.get('exonEnd', row.get('downstreamEE', 0)),
                             row.get('downstreamEE', row.get('exonEnd', 0)))
                    if start and end:
                        start = max(1, int(start) - 500)
                        end = int(end) + 500
                        region = f"{chrom}:{start}-{end}"
                        success = run_ggsashimi(region, f"{gene_name}_{event}", "top_splicing")
                        if success:
                            splicing_genes_processed.add(gene_name)
                    if len(splicing_genes_processed) >= 15:
                        break
            except Exception as e:
                print(f"  Error processing {jc_file}: {e}")
        if len(splicing_genes_processed) >= 15:
            break
    if len(splicing_genes_processed) >= 15:
        break

print(f"Generated {len(splicing_genes_processed)} splicing event plots")

# Process user-specified genes
print("\n=== Processing User-Specified Genes ===")
if os.path.exists(USER_GENES_FILE):
    with open(USER_GENES_FILE, 'r') as f:
        user_genes = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    for gene_name in user_genes:
        region, strand = parse_gtf_for_gene(GTF_FILE, gene_name)
        if region:
            run_ggsashimi(region, gene_name, "user_genes")
        else:
            print(f"  Gene not found in GTF: {gene_name}")
    print(f"Processed {len(user_genes)} user-specified genes")
else:
    print(f"User genes file not found: {USER_GENES_FILE}")
    print("Create this file with gene symbols (one per line) to generate custom plots")

print("\n=== Sashimi plot generation complete ===")
