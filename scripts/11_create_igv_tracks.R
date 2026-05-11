#!/usr/bin/env Rscript
# =============================================================================
# Script: 11_create_igv_tracks.R
# Purpose: Create BED files for IGV visualization of significant splicing events
# Output: results/microexons_specific/
#   - {Comparison}_significant_SE.bed - Deduplicated skipped exon events (strongest |dPSI| per exon)
#   - {Comparison}_significant_SE.csv - All events with flanking exon context
#   - {Comparison}_MICROEXONS.bed - Deduplicated microexon events (≤30bp)
#   - {Comparison}_MICROEXONS.csv - All microexon events with flanking exon context
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(dplyr)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/microexons_specific")

# Significance thresholds
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
MICROEXON_MAX_LENGTH <- 30  # bp

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Main comparisons (vs Parental baseline)
COMPARISONS <- c(
    "KO_vs_Parental",
    "Neg_vs_Parental",
    "Pos_vs_Parental"
)

# -----------------------------------------------------------------------------
# Function: Convert dPSI to BED score (0-1000 scale)
# Higher |dPSI| = lower score (darker color in IGV)
# -----------------------------------------------------------------------------
dpsi_to_score <- function(dpsi) {
    # Scale |dPSI| (0-1) to score (0-1000)
    # Invert so that stronger effects have lower scores (darker in IGV)
    score <- round((1 - abs(dpsi)) * 250)
    score <- pmax(0, pmin(1000, score))  # Clamp to valid range
    return(score)
}

# -----------------------------------------------------------------------------
# Function: Create BED file from rMATS SE output
# -----------------------------------------------------------------------------
create_bed_from_se <- function(comparison, filter_microexons = FALSE) {
    # Read SE (Skipped Exon) file
    se_file <- file.path(SPLICING_DIR, comparison, "SE.MATS.JC.txt")

    if (!file.exists(se_file)) {
        warning(sprintf("File not found: %s", se_file))
        return(NULL)
    }

    # Read data
    se <- read.delim(se_file, stringsAsFactors = FALSE)

    # Filter for significant events
    se_sig <- se %>%
        filter(FDR < FDR_THRESHOLD, abs(IncLevelDifference) >= DPSI_THRESHOLD)

    if (nrow(se_sig) == 0) {
        warning(sprintf("No significant events for %s", comparison))
        return(NULL)
    }

    # Calculate exon length and event_id (unique triplet identifier)
    se_sig <- se_sig %>%
        mutate(
            exon_length = exonEnd - exonStart_0base,
            event_id = sprintf("%s:%d-%d|up:%d-%d|down:%d-%d",
                               chr, exonStart_0base, exonEnd,
                               upstreamES, upstreamEE, downstreamES, downstreamEE)
        )

    # Filter for microexons if requested
    if (filter_microexons) {
        se_sig <- se_sig %>%
            filter(exon_length <= MICROEXON_MAX_LENGTH)

        if (nrow(se_sig) == 0) {
            warning(sprintf("No microexons for %s", comparison))
            return(NULL)
        }
    }

    # CSV: full event context (all events, including duplicates per exon region)
    csv_df <- se_sig %>%
        select(event_id, geneSymbol, chr, strand, exonStart_0base, exonEnd, exon_length,
               upstreamES, upstreamEE, downstreamES, downstreamEE,
               FDR, IncLevelDifference) %>%
        arrange(chr, exonStart_0base)

    # BED: deduplicated — keep strongest |dPSI| per unique exon region
    bed <- se_sig %>%
        group_by(chr, exonStart_0base, exonEnd, strand, geneSymbol) %>%
        slice_max(abs(IncLevelDifference), n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(
            chrom = chr,
            chromStart = exonStart_0base,
            chromEnd = exonEnd,
            name = if (filter_microexons) {
                sprintf('"%s"_%dbp_dPSI=%s', geneSymbol, exon_length, IncLevelDifference)
            } else {
                sprintf('"%s"_dPSI=%s', geneSymbol, IncLevelDifference)
            },
            score = dpsi_to_score(IncLevelDifference),
            strand_col = strand
        ) %>%
        select(chrom, chromStart, chromEnd, name, score, strand_col) %>%
        arrange(chrom, chromStart)

    return(list(csv = csv_df, bed = bed))
}

# -----------------------------------------------------------------------------
# Main Processing
# -----------------------------------------------------------------------------
cat("=" , rep("=", 69), "\n", sep = "")
cat("Creating IGV BED Tracks for Splicing Events\n")
cat("=" , rep("=", 69), "\n", sep = "")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("FDR threshold: %s\n", FDR_THRESHOLD))
cat(sprintf("|dPSI| threshold: %s\n", DPSI_THRESHOLD))
cat(sprintf("Microexon max length: %d bp\n", MICROEXON_MAX_LENGTH))
cat("\n")

for (comparison in COMPARISONS) {
    cat(sprintf("Processing %s...\n", comparison))

    # 1. All significant SE events
    result_all <- create_bed_from_se(comparison, filter_microexons = FALSE)
    if (!is.null(result_all)) {
        # Write deduplicated BED
        bed_file <- file.path(OUTPUT_DIR, sprintf("%s_significant_SE.bed", comparison))
        write.table(result_all$bed, bed_file,
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        # Write full-context CSV
        csv_file <- file.path(OUTPUT_DIR, sprintf("%s_significant_SE.csv", comparison))
        write.csv(result_all$csv, csv_file, row.names = FALSE)
        cat(sprintf("  %s_significant_SE.bed: %d events (deduplicated from %d)\n",
                    comparison, nrow(result_all$bed), nrow(result_all$csv)))
    }

    # 2. Microexons only
    result_micro <- create_bed_from_se(comparison, filter_microexons = TRUE)
    if (!is.null(result_micro)) {
        # Write deduplicated BED
        bed_file <- file.path(OUTPUT_DIR, sprintf("%s_MICROEXONS.bed", comparison))
        write.table(result_micro$bed, bed_file,
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        # Write full-context CSV
        csv_file <- file.path(OUTPUT_DIR, sprintf("%s_MICROEXONS.csv", comparison))
        write.csv(result_micro$csv, csv_file, row.names = FALSE)
        cat(sprintf("  %s_MICROEXONS.bed: %d microexons (deduplicated from %d)\n",
                    comparison, nrow(result_micro$bed), nrow(result_micro$csv)))
    }
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n")
cat("=" , rep("=", 69), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 69), "\n", sep = "")

# List created files
cat("\nCreated files:\n")
files <- list.files(OUTPUT_DIR, pattern = "\\.(bed|csv)$", full.names = TRUE)
for (f in files) {
    info <- file.info(f)
    cat(sprintf("  %s (%s)\n", basename(f), format(info$size, big.mark = ",")))
}

cat("\n")
cat("BED file format (deduplicated — strongest |dPSI| per exon region):\n")
cat("  Column 1: Chromosome\n")
cat("  Column 2: Start position (0-based)\n")
cat("  Column 3: End position\n")
cat("  Column 4: Name (gene symbol + dPSI)\n")
cat("  Column 5: Score (0-1000, lower = stronger effect)\n")
cat("  Column 6: Strand (+/-)\n")

cat("\n")
cat("CSV file columns (all events with flanking exon context):\n")
cat("  event_id: Unique triplet identifier (skipped|upstream|downstream)\n")
cat("  geneSymbol, chr, strand, exonStart_0base, exonEnd, exon_length\n")
cat("  upstreamES, upstreamEE: Upstream flanking exon coordinates\n")
cat("  downstreamES, downstreamEE: Downstream flanking exon coordinates\n")
cat("  FDR, IncLevelDifference: Significance and effect size\n")

cat("\n")
cat("To load in IGV:\n")
cat("  File -> Load from File -> Select .bed file\n")
cat("\n")
cat("Done!\n")
