#!/usr/bin/env Rscript
# =============================================================================
# 14b_neg_exclusive_genes.R
# Stage 1: Find genes with microexon splicing events significant in
# Neg_vs_Parental but NOT in Pos_vs_Parental. Outputs the gene list and
# associated event details for downstream analysis (14c).
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/14_todo_analysis/task2_neg_exclusive_microexon")

FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
MIN_AVG_READS <- 10
MICROEXON_MAX <- 30

EVENT_TYPES <- c("SE", "A3SS", "A5SS", "MXE", "RI")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Utility Functions
# =============================================================================
sum_counts <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(0)
    sum(as.numeric(unlist(strsplit(as.character(x), ","))), na.rm = TRUE)
}

load_rmats_results <- function(comparison, event_type) {
    file_path <- file.path(SPLICING_DIR, comparison, paste0(event_type, ".MATS.JC.txt"))
    if (!file.exists(file_path)) {
        warning(paste("File not found:", file_path))
        return(NULL)
    }
    df <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    dup_cols <- duplicated(colnames(df))
    if (any(dup_cols)) df <- df[, !dup_cols]
    df$event_type <- event_type
    df$comparison <- comparison
    if (all(c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2") %in% colnames(df))) {
        df$avg_reads <- (
            sapply(df$IJC_SAMPLE_1, sum_counts) +
            sapply(df$SJC_SAMPLE_1, sum_counts) +
            sapply(df$IJC_SAMPLE_2, sum_counts) +
            sapply(df$SJC_SAMPLE_2, sum_counts)
        ) / 6
        df <- df[df$avg_reads >= MIN_AVG_READS, ]
    }
    return(df)
}

create_event_key <- function(df, event_type) {
    switch(event_type,
        SE = paste(df$chr, df$exonStart_0base, df$exonEnd,
                   df$upstreamEE, df$downstreamES, "SE", sep = ":"),
        A3SS = paste(df$chr, df$longExonStart_0base, df$longExonEnd,
                     df$shortES, df$shortEE, df$flankingES, df$flankingEE, "A3SS", sep = ":"),
        A5SS = paste(df$chr, df$longExonStart_0base, df$longExonEnd,
                     df$shortES, df$shortEE, df$flankingES, df$flankingEE, "A5SS", sep = ":"),
        MXE = paste(df$chr, df$`1stExonStart_0base`, df$`1stExonEnd`,
                    df$`2ndExonStart_0base`, df$`2ndExonEnd`,
                    df$upstreamEE, df$downstreamES, "MXE", sep = ":"),
        RI = paste(df$chr, df$riExonStart_0base, df$riExonEnd,
                   df$upstreamEE, df$downstreamES, "RI", sep = ":"),
        stop(paste("Unknown event type:", event_type))
    )
}

compute_exon_length <- function(df, event_type) {
    switch(event_type,
        SE = df$exonEnd - df$exonStart_0base,
        A3SS = df$longExonEnd - df$longExonStart_0base,
        A5SS = df$longExonEnd - df$longExonStart_0base,
        MXE = df$`1stExonEnd` - df$`1stExonStart_0base`,
        RI = df$riExonEnd - df$riExonStart_0base,
        stop(paste("Unknown event type:", event_type))
    )
}

classify_size <- function(lengths) {
    factor(
        case_when(
            lengths <= MICROEXON_MAX ~ "Microexon (0-30bp)",
            lengths <= 50 ~ "Small (30-50bp)",
            TRUE ~ "Regular (>50bp)"
        ),
        levels = c("Microexon (0-30bp)", "Small (30-50bp)", "Regular (>50bp)")
    )
}

# =============================================================================
# Main: Load all events, find neg-exclusive microexon genes
# =============================================================================
cat("========================================\n")
cat("14b: Neg-Exclusive Microexon Gene Extraction\n")
cat("========================================\n")

load_all_events <- function(comparison) {
    all_events <- list()
    for (et in EVENT_TYPES) {
        df <- load_rmats_results(comparison, et)
        if (!is.null(df) && nrow(df) > 0) {
            df$event_key <- create_event_key(df, et)
            df$exon_length <- compute_exon_length(df, et)
            df$size_category <- classify_size(df$exon_length)
            df$significant <- df$FDR < FDR_THRESHOLD & abs(df$IncLevelDifference) > DPSI_THRESHOLD
            df$dPSI <- df$IncLevelDifference
            all_events[[et]] <- df
        }
    }
    bind_rows(all_events)
}

neg_all <- load_all_events("Neg_vs_Parental")
pos_all <- load_all_events("Pos_vs_Parental")

neg_sig <- neg_all %>% filter(significant)
pos_sig <- pos_all %>% filter(significant)

cat(sprintf("  Neg significant: %d events\n", nrow(neg_sig)))
cat(sprintf("  Pos significant: %d events\n", nrow(pos_sig)))

# Significant SE microexons in each comparison
neg_micro <- neg_sig %>% filter(event_type == "SE", exon_length <= MICROEXON_MAX)
pos_micro <- pos_sig %>% filter(event_type == "SE", exon_length <= MICROEXON_MAX)

cat(sprintf("  Significant SE microexons: Neg=%d, Pos=%d\n",
            nrow(neg_micro), nrow(pos_micro)))

# Neg-exclusive microexon event keys
neg_excl_keys <- setdiff(neg_micro$event_key, pos_micro$event_key)
cat(sprintf("  Neg-exclusive microexon events: %d\n", length(neg_excl_keys)))

if (length(neg_excl_keys) == 0) {
    cat("  No Neg-exclusive microexon events found. Nothing to output.\n")
    quit(status = 0)
}

neg_excl_events <- neg_micro %>% filter(event_key %in% neg_excl_keys)
neg_excl_genes <- unique(neg_excl_events$GeneID)
cat(sprintf("  Unique genes with neg-exclusive microexons: %d\n", length(neg_excl_genes)))

# =============================================================================
# Save outputs
# =============================================================================

# 1. Neg-exclusive microexon events (full detail)
neg_excl_detail <- neg_excl_events %>%
    select(event_key, event_type, GeneID, geneSymbol, chr, strand,
           exonStart_0base, exonEnd, exon_length, dPSI, FDR) %>%
    arrange(geneSymbol, exonStart_0base)
write.csv(neg_excl_detail,
          file.path(OUTPUT_DIR, "neg_exclusive_microexon_events.csv"), row.names = FALSE)

# 2. Gene list (compact, for consumption by 14c)
gene_summary <- neg_excl_detail %>%
    group_by(GeneID, geneSymbol) %>%
    summarize(
        n_neg_excl_microexons = n(),
        microexon_coords = paste(paste0(chr, ":", exonStart_0base, "-", exonEnd), collapse = "; "),
        mean_dPSI = round(mean(dPSI), 3),
        dPSI_values = paste(round(dPSI, 3), collapse = "; "),
        .groups = "drop"
    ) %>%
    arrange(desc(n_neg_excl_microexons))
write.csv(gene_summary,
          file.path(OUTPUT_DIR, "neg_exclusive_microexon_genes.csv"), row.names = FALSE)

# 3. Also save the gene IDs as a simple text file for easy consumption
writeLines(neg_excl_genes, file.path(OUTPUT_DIR, "neg_exclusive_gene_ids.txt"))

# 4. Save all Neg events for these genes (for Neg-Pos scatter in 14c)
neg_events_for_genes <- neg_all %>%
    filter(GeneID %in% neg_excl_genes) %>%
    select(event_key, event_type, GeneID, geneSymbol, exon_length,
           size_category, dPSI, FDR, significant)
write.csv(neg_events_for_genes,
          file.path(OUTPUT_DIR, "neg_all_events_for_genes.csv"), row.names = FALSE)

cat("\n  Output files:\n")
cat(sprintf("    %s/neg_exclusive_microexon_events.csv (%d events)\n",
            OUTPUT_DIR, nrow(neg_excl_detail)))
cat(sprintf("    %s/neg_exclusive_microexon_genes.csv (%d genes)\n",
            OUTPUT_DIR, nrow(gene_summary)))
cat(sprintf("    %s/neg_exclusive_gene_ids.txt (%d gene IDs)\n",
            OUTPUT_DIR, length(neg_excl_genes)))
cat(sprintf("    %s/neg_all_events_for_genes.csv (%d events)\n",
            OUTPUT_DIR, nrow(neg_events_for_genes)))

cat("\n14b complete.\n")
