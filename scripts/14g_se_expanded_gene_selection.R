#!/usr/bin/env Rscript
# =============================================================================
# 14g_se_expanded_gene_selection.R
#
# SE-ONLY version of 14g_expanded_gene_selection.R
#
# Expanded gene selection for 3 cumulative groups, each using Neg_vs_Parental
# significant SE (Skipped Exon) events only that are neg_only OR shared (both).
#
# Events are classified by SE exon_length:
#   Group 1: Genes with microexon SE (0-30bp) significant in Neg
#   Group 2: Group 1 + genes with small exon SE (30-50bp) significant in Neg
#   Group 3: Group 2 + genes with regular SE (>50bp) significant in Neg
#
# "neg_only + both" = significant in Neg_vs_Parental regardless of Pos status
#
# Output directory: results/14_todo_analysis/task5_expanded_groups_se/
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
OUTPUT_BASE <- file.path(BASE_DIR, "results/14_todo_analysis/task5_expanded_groups_se")

FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
MIN_AVG_READS <- 10
MICROEXON_MAX <- 30
SMALL_MAX <- 50

EVENT_TYPES <- c("SE", "A3SS", "A5SS", "MXE", "RI")

# Group definitions (SE events only)
GROUPS <- list(
    group1_neg_microexon = list(
        label = "Microexon SE (0-30bp)",
        description = "Genes with microexon SE events significant in Neg (neg_only + both)",
        max_size = MICROEXON_MAX
    ),
    group2_neg_micro_small = list(
        label = "Microexon + Small SE (0-50bp)",
        description = "Genes with microexon or small exon SE events significant in Neg",
        max_size = SMALL_MAX
    ),
    group3_neg_all_sizes = list(
        label = "All SE Sizes",
        description = "Genes with any SE event significant in Neg (all sizes)",
        max_size = Inf
    )
)

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
            lengths <= SMALL_MAX ~ "Small (30-50bp)",
            TRUE ~ "Regular (>50bp)"
        ),
        levels = c("Microexon (0-30bp)", "Small (30-50bp)", "Regular (>50bp)")
    )
}

# =============================================================================
# Main: Load all events from both comparisons
# =============================================================================
cat("========================================================\n")
cat("14g_se: Expanded Gene Selection (3 Cumulative Groups)\n")
cat("  *** SE events only ***\n")
cat("========================================================\n\n")

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

cat("Loading Neg_vs_Parental events...\n")
neg_all <- load_all_events("Neg_vs_Parental")
cat("Loading Pos_vs_Parental events...\n")
pos_all <- load_all_events("Pos_vs_Parental")

neg_sig <- neg_all %>% filter(significant)
pos_sig <- pos_all %>% filter(significant)

cat(sprintf("  Neg total events: %d, significant: %d\n", nrow(neg_all), nrow(neg_sig)))
cat(sprintf("  Pos total events: %d, significant: %d\n", nrow(pos_all), nrow(pos_sig)))

# Filter to SE events only for group selection
neg_sig_se <- neg_sig %>% filter(event_type == "SE")
cat(sprintf("\n  Neg significant SE events: %d\n", nrow(neg_sig_se)))
cat(sprintf("    Microexon (0-30bp): %d\n", sum(neg_sig_se$exon_length <= MICROEXON_MAX)))
cat(sprintf("    Small (30-50bp):    %d\n", sum(neg_sig_se$exon_length > MICROEXON_MAX & neg_sig_se$exon_length <= SMALL_MAX)))
cat(sprintf("    Regular (>50bp):    %d\n", sum(neg_sig_se$exon_length > SMALL_MAX)))

# =============================================================================
# Process each group
# =============================================================================

for (group_name in names(GROUPS)) {
    group <- GROUPS[[group_name]]
    group_dir <- file.path(OUTPUT_BASE, group_name)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

    cat(sprintf("\n--- %s: %s ---\n", group_name, group$label))
    cat(sprintf("  %s\n", group$description))

    # Select SE events for this group (Neg significant, within size cutoff)
    if (is.infinite(group$max_size)) {
        selection_events <- neg_sig_se
    } else {
        selection_events <- neg_sig_se %>% filter(exon_length <= group$max_size)
    }

    selection_genes <- unique(selection_events$GeneID)
    cat(sprintf("  Selection SE events: %d\n", nrow(selection_events)))
    cat(sprintf("  Unique genes: %d\n", length(selection_genes)))

    if (length(selection_genes) == 0) {
        cat("  WARNING: No genes found for this group. Skipping.\n")
        next
    }

    # --- 1. Save selection events ---
    selection_detail <- selection_events %>%
        select(event_key, event_type, GeneID, geneSymbol, chr, strand,
               exonStart_0base, exonEnd, exon_length, size_category, dPSI, FDR) %>%
        arrange(geneSymbol, exonStart_0base)
    write.csv(selection_detail,
              file.path(group_dir, "selection_events.csv"), row.names = FALSE)

    # --- 2. Gene list with summary ---
    gene_summary <- selection_detail %>%
        group_by(GeneID, geneSymbol) %>%
        summarize(
            n_selection_events = n(),
            size_categories = paste(sort(unique(as.character(size_category))), collapse = "; "),
            event_coords = paste(paste0(chr, ":", exonStart_0base, "-", exonEnd), collapse = "; "),
            mean_dPSI = round(mean(dPSI), 3),
            .groups = "drop"
        ) %>%
        arrange(desc(n_selection_events))
    write.csv(gene_summary,
              file.path(group_dir, "gene_list.csv"), row.names = FALSE)

    # --- 3. Gene IDs text file ---
    writeLines(selection_genes, file.path(group_dir, "gene_ids.txt"))

    # --- 4. ALL Neg events for these genes (all event types, all sizes) ---
    neg_events_for_genes <- neg_all %>%
        filter(GeneID %in% selection_genes) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR, significant)
    write.csv(neg_events_for_genes,
              file.path(group_dir, "neg_all_events_for_genes.csv"), row.names = FALSE)
    cat(sprintf("  Neg events (all types) for these genes: %d\n", nrow(neg_events_for_genes)))

    # --- 5. ALL Pos events for these genes ---
    pos_events_for_genes <- pos_all %>%
        filter(GeneID %in% selection_genes) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR, significant)
    write.csv(pos_events_for_genes,
              file.path(group_dir, "pos_all_events_for_genes.csv"), row.names = FALSE)
    cat(sprintf("  Pos events (all types) for these genes: %d\n", nrow(pos_events_for_genes)))

    # --- 6. Per-gene Pos summary ---
    gene_pos_summary <- pos_events_for_genes %>%
        group_by(GeneID, geneSymbol) %>%
        summarize(
            n_events_pos = n(),
            n_sig_pos = sum(significant),
            n_event_types = n_distinct(event_type),
            event_types = paste(sort(unique(event_type)), collapse = ","),
            n_microexon = sum(size_category == "Microexon (0-30bp)" & event_type == "SE"),
            mean_abs_dPSI = round(mean(abs(dPSI), na.rm = TRUE), 3),
            .groups = "drop"
        )

    gene_pos_summary <- gene_pos_summary %>%
        left_join(
            gene_summary %>% select(GeneID, n_selection_events, mean_dPSI),
            by = "GeneID"
        ) %>%
        rename(selection_mean_dPSI = mean_dPSI) %>%
        arrange(desc(n_selection_events))

    write.csv(gene_pos_summary,
              file.path(group_dir, "gene_pos_summary.csv"), row.names = FALSE)

    cat(sprintf("  Gene summary saved: %d genes\n", nrow(gene_pos_summary)))
    cat(sprintf("  Genes with sig Pos events: %d (%.1f%%)\n",
                sum(gene_pos_summary$n_sig_pos > 0),
                100 * sum(gene_pos_summary$n_sig_pos > 0) / nrow(gene_pos_summary)))
}

# =============================================================================
# Cross-group comparison summary
# =============================================================================
cat("\n========================================================\n")
cat("  CROSS-GROUP SUMMARY (SE only)\n")
cat("========================================================\n\n")

for (group_name in names(GROUPS)) {
    group <- GROUPS[[group_name]]
    group_dir <- file.path(OUTPUT_BASE, group_name)
    gene_ids_file <- file.path(group_dir, "gene_ids.txt")
    if (file.exists(gene_ids_file)) {
        n_genes <- length(readLines(gene_ids_file))
        sel_events <- read.csv(file.path(group_dir, "selection_events.csv"))
        pos_summ <- read.csv(file.path(group_dir, "gene_pos_summary.csv"))
        cat(sprintf("  %-30s: %4d genes, %5d selection SE events, %d (%.0f%%) with sig Pos\n",
                    group$label, n_genes, nrow(sel_events),
                    sum(pos_summ$n_sig_pos > 0),
                    100 * sum(pos_summ$n_sig_pos > 0) / n_genes))
    }
}

cat(sprintf("\n  Output base: %s\n", OUTPUT_BASE))
cat("\n14g_se complete.\n")
