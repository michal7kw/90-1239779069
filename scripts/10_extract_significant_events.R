#!/usr/bin/env Rscript
# =============================================================================
# Script: 10_extract_significant_events.R
# Purpose: Extract and consolidate all significant splicing events from rMATS
# Author: Recreated based on output file structure
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})
# Note: Using base R read.delim() and write.csv() - no readr needed

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/10_significant_events")

# Significance thresholds
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define comparisons and event types
COMPARISONS <- c(
    "KO_vs_Parental",
    "Neg_vs_Parental",
    "Pos_vs_Parental",
    "KO_vs_Neg",
    "KO_vs_Pos",
    "Pos_vs_Neg"
)

EVENT_TYPES <- list(
    SE = list(
        file = "SE.MATS.JC.txt",
        desc = "Skipped Exon",
        coord_cols = c("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
                       "downstreamES", "downstreamEE")
    ),
    A3SS = list(
        file = "A3SS.MATS.JC.txt",
        desc = "Alternative 3' Splice Site",
        coord_cols = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE",
                       "flankingES", "flankingEE")
    ),
    A5SS = list(
        file = "A5SS.MATS.JC.txt",
        desc = "Alternative 5' Splice Site",
        coord_cols = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE",
                       "flankingES", "flankingEE")
    ),
    MXE = list(
        file = "MXE.MATS.JC.txt",
        desc = "Mutually Exclusive Exons",
        coord_cols = c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
                       "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES",
                       "downstreamEE")
    ),
    RI = list(
        file = "RI.MATS.JC.txt",
        desc = "Retained Intron",
        coord_cols = c("riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE",
                       "downstreamES", "downstreamEE")
    )
)

# -----------------------------------------------------------------------------
# Function: Read and filter rMATS output
# -----------------------------------------------------------------------------
read_rmats_significant <- function(comparison, event_type, event_info) {
    file_path <- file.path(SPLICING_DIR, comparison, event_info$file)

    if (!file.exists(file_path)) {
        message(sprintf("  File not found: %s", file_path))
        return(NULL)
    }

    # Read rMATS output
    df <- read.delim(file_path, stringsAsFactors = FALSE)

    # Filter for significant events
    df_sig <- df %>%
        filter(FDR < FDR_THRESHOLD, abs(IncLevelDifference) >= DPSI_THRESHOLD)

    if (nrow(df_sig) == 0) {
        return(NULL)
    }

    # Add metadata columns
    df_sig <- df_sig %>%
        mutate(
            Comparison = comparison,
            EventType = event_type,
            EventDescription = event_info$desc,
            DeltaPSI = IncLevelDifference,
            PSI_Group1 = IncLevel1,
            PSI_Group2 = IncLevel2
        )

    # Standardize column names
    df_sig <- df_sig %>%
        rename(
            GeneSymbol = geneSymbol,
            Chromosome = chr
        )

    # Select and reorder columns
    # Common columns first, then event-specific coordinate columns
    common_cols <- c("Comparison", "EventType", "EventDescription", "GeneSymbol",
                     "GeneID", "Chromosome", "strand", "FDR", "DeltaPSI",
                     "PSI_Group1", "PSI_Group2")

    # Get all coordinate columns that exist in this event type
    # Note: R adds 'X' prefix to column names starting with numbers
    all_coord_cols <- c("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
                        "downstreamES", "downstreamEE", "longExonStart_0base",
                        "longExonEnd", "shortES", "shortEE", "flankingES",
                        "flankingEE", "X1stExonStart_0base", "X1stExonEnd",
                        "X2ndExonStart_0base", "X2ndExonEnd", "riExonStart_0base",
                        "riExonEnd")

    coord_cols_present <- intersect(all_coord_cols, names(df_sig))

    # Select available columns
    available_cols <- intersect(c(common_cols, coord_cols_present), names(df_sig))
    df_sig <- df_sig %>% select(all_of(available_cols))

    # Rename strand column
    if ("strand" %in% names(df_sig)) {
        df_sig <- df_sig %>% rename(Strand = strand)
    }

    return(df_sig)
}

# -----------------------------------------------------------------------------
# Main Processing
# -----------------------------------------------------------------------------
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("Extracting Significant Splicing Events from rMATS\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat(sprintf("FDR threshold: %s\n", FDR_THRESHOLD))
cat(sprintf("|dPSI| threshold: %s\n", DPSI_THRESHOLD))
cat("\n")

# Collect all significant events
all_events <- list()
comparison_events <- list()

for (comparison in COMPARISONS) {
    cat(sprintf("Processing %s...\n", comparison))
    comp_events <- list()

    for (event_type in names(EVENT_TYPES)) {
        event_info <- EVENT_TYPES[[event_type]]
        df <- read_rmats_significant(comparison, event_type, event_info)

        if (!is.null(df) && nrow(df) > 0) {
            cat(sprintf("  %s: %d significant events\n", event_type, nrow(df)))
            comp_events[[event_type]] <- df
            all_events[[paste(comparison, event_type, sep = "_")]] <- df
        }
    }

    # Combine events for this comparison
    if (length(comp_events) > 0) {
        comparison_events[[comparison]] <- bind_rows(comp_events)
    }
}

# -----------------------------------------------------------------------------
# Create Output Files
# -----------------------------------------------------------------------------
cat("\n")
cat("Creating output files...\n")

# 1. All significant events (complete data)
all_events_df <- bind_rows(all_events)
all_events_df <- all_events_df %>%
    arrange(Comparison, EventType, FDR)

write.csv(all_events_df,
          file.path(OUTPUT_DIR, "all_significant_splicing_events.csv"),
          row.names = FALSE)
cat(sprintf("  all_significant_splicing_events.csv: %d events\n", nrow(all_events_df)))

# 2. Simplified version (key columns only)
simple_df <- all_events_df %>%
    select(Comparison, EventType, GeneSymbol, Chromosome, Strand, FDR, DeltaPSI) %>%
    arrange(Comparison, desc(abs(DeltaPSI)))

write.csv(simple_df,
          file.path(OUTPUT_DIR, "significant_events_simple.csv"),
          row.names = FALSE)
cat(sprintf("  significant_events_simple.csv: %d events\n", nrow(simple_df)))

# 3. Per-comparison files
for (comparison in names(comparison_events)) {
    df <- comparison_events[[comparison]]
    filename <- sprintf("%s_significant_events.csv", comparison)
    write.csv(df, file.path(OUTPUT_DIR, filename), row.names = FALSE)
    cat(sprintf("  %s: %d events\n", filename, nrow(df)))
}

# 4. Gene-level summary
genes_summary <- all_events_df %>%
    group_by(GeneSymbol) %>%
    summarise(
        Comparisons = paste(unique(Comparison), collapse = ","),
        EventTypes = paste(unique(EventType), collapse = ","),
        TotalEvents = n(),
        Mean_dPSI = round(mean(DeltaPSI), 3),
        Min_dPSI = round(min(DeltaPSI), 3),
        Max_dPSI = round(max(DeltaPSI), 3),
        .groups = "drop"
    ) %>%
    arrange(desc(TotalEvents))

write.csv(genes_summary,
          file.path(OUTPUT_DIR, "genes_with_splicing_events.csv"),
          row.names = FALSE)
cat(sprintf("  genes_with_splicing_events.csv: %d genes\n", nrow(genes_summary)))

# 5. Summary by comparison and event type
summary_stats <- all_events_df %>%
    group_by(Comparison, EventType) %>%
    summarise(
        Count = n(),
        Mean_dPSI = round(mean(DeltaPSI), 3),
        Min_dPSI = round(min(DeltaPSI), 3),
        Max_dPSI = round(max(DeltaPSI), 3),
        .groups = "drop"
    ) %>%
    arrange(Comparison, EventType)

write.csv(summary_stats,
          file.path(OUTPUT_DIR, "summary_by_comparison_eventtype.csv"),
          row.names = FALSE)
cat(sprintf("  summary_by_comparison_eventtype.csv: %d rows\n", nrow(summary_stats)))

# -----------------------------------------------------------------------------
# 6. Create IGV-ready coordinate files
# -----------------------------------------------------------------------------
cat("\nCreating IGV coordinate files...\n")

# Function to add IGV coordinates for each event type separately
add_igv_coords <- function(df, event_type) {
    if (event_type == "SE") {
        df <- df %>%
            mutate(
                event_start = exonStart_0base,
                event_end = exonEnd,
                ExonLength = exonEnd - exonStart_0base
            )
    } else if (event_type %in% c("A3SS", "A5SS")) {
        df <- df %>%
            mutate(
                event_start = longExonStart_0base,
                event_end = longExonEnd,
                ExonLength = NA_real_
            )
    } else if (event_type == "MXE") {
        # R adds 'X' prefix to column names starting with numbers
        df <- df %>%
            mutate(
                event_start = X1stExonStart_0base,
                event_end = X2ndExonEnd,
                ExonLength = NA_real_
            )
    } else if (event_type == "RI") {
        df <- df %>%
            mutate(
                event_start = riExonStart_0base,
                event_end = riExonEnd,
                ExonLength = NA_real_
            )
    }

    df %>%
        mutate(
            IGV_coordinate = paste0(Chromosome, ":",
                                    pmax(0, event_start - 500), "-",
                                    event_end + 500),
            IGV_narrow = paste0(Chromosome, ":",
                                pmax(0, event_start - 100), "-",
                                event_end + 100)
        )
}

# Process each event type separately and combine
get_event_coords <- function(df) {
    event_types <- unique(df$EventType)
    result_list <- list()

    for (et in event_types) {
        et_df <- df %>% filter(EventType == et)
        if (nrow(et_df) > 0) {
            result_list[[et]] <- add_igv_coords(et_df, et)
        }
    }

    bind_rows(result_list)
}

# Add IGV coordinates to all events
all_events_igv <- get_event_coords(all_events_df)

# 6a. Full IGV coordinates file (TSV for easy viewing)
igv_full <- all_events_igv %>%
    select(
        Gene = GeneSymbol,
        Comparison,
        EventType,
        dPSI = DeltaPSI,
        FDR,
        Chromosome,
        Start = event_start,
        End = event_end,
        ExonLength,
        IGV_coordinate,
        IGV_narrow
    ) %>%
    arrange(Comparison, desc(abs(dPSI))) %>%
    mutate(
        dPSI = round(dPSI, 3),
        FDR = signif(FDR, 3)
    )

write.table(igv_full,
            file.path(OUTPUT_DIR, "IGV_coordinates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  IGV_coordinates.tsv: %d events\n", nrow(igv_full)))

# 6b. Simple copy-paste file (just gene label and coordinate)
igv_simple <- igv_full %>%
    mutate(
        Label = paste0(Gene, " | ", EventType, " | dPSI=", dPSI, " | ", Comparison)
    ) %>%
    select(Label, IGV_coordinate) %>%
    distinct()

write.table(igv_simple,
            file.path(OUTPUT_DIR, "IGV_coordinates_simple.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
cat(sprintf("  IGV_coordinates_simple.txt: %d unique coordinates\n", nrow(igv_simple)))

# 6c. Microexon-specific file (SE events ≤30bp)
igv_microexons <- igv_full %>%
    filter(EventType == "SE", !is.na(ExonLength), ExonLength <= 30) %>%
    arrange(Comparison, desc(abs(dPSI)))

write.table(igv_microexons,
            file.path(OUTPUT_DIR, "IGV_coordinates_microexons.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  IGV_coordinates_microexons.tsv: %d microexon events\n", nrow(igv_microexons)))

# 6d. Top hits file (|dPSI| >= 0.3)
igv_top <- igv_full %>%
    filter(abs(dPSI) >= 0.3) %>%
    arrange(desc(abs(dPSI)))

write.table(igv_top,
            file.path(OUTPUT_DIR, "IGV_coordinates_top_hits.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  IGV_coordinates_top_hits.tsv: %d events with |dPSI| >= 0.3\n", nrow(igv_top)))

# 6e. Per-comparison IGV files for easy navigation
for (comp in unique(igv_full$Comparison)) {
    comp_igv <- igv_full %>%
        filter(Comparison == comp) %>%
        arrange(desc(abs(dPSI)))

    # Simple format for this comparison
    comp_simple <- comp_igv %>%
        mutate(Label = paste0(Gene, " | ", EventType, " | dPSI=", dPSI)) %>%
        select(Label, IGV_coordinate)

    write.table(comp_simple,
                file.path(OUTPUT_DIR, paste0("IGV_", comp, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}
cat(sprintf("  IGV_{Comparison}.txt: 6 per-comparison files\n"))

# -----------------------------------------------------------------------------
# Print Summary
# -----------------------------------------------------------------------------
cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("SUMMARY\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")

cat("\nTotal significant events by comparison:\n")
all_events_df %>%
    group_by(Comparison) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count)) %>%
    print()

cat("\nTotal significant events by event type:\n")
all_events_df %>%
    group_by(EventType) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count)) %>%
    print()

cat("\nTop 10 genes with most splicing events:\n")
genes_summary %>%
    head(10) %>%
    select(GeneSymbol, TotalEvents, EventTypes) %>%
    print()

cat("\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("Done!\n")
