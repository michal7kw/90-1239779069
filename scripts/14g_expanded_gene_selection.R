#!/usr/bin/env Rscript
# =============================================================================
# 14g_expanded_gene_selection.R
#
# Expanded gene selection for 3 cumulative groups, each using Neg_vs_Parental
# significant events (ALL event types: SE, A3SS, A5SS, MXE, RI) that are
# neg_only OR shared (both).
#
# Events are classified by computed exon_length, matching the Venn diagram
# (14a) size classification:
#   Group 1: Genes with microexon events (0-30bp) significant in Neg
#   Group 2: Group 1 + genes with small exon events (30-50bp) significant in Neg
#   Group 3: Group 2 + genes with regular events (>50bp) significant in Neg
#
# "neg_only + both" = significant in Neg_vs_Parental regardless of Pos status
# (i.e. all Neg-significant events, whether or not Pos is also significant)
#
# For each group, saves the same file structure as 14b+14c so that 14h can
# run 14f-style characterization downstream.
#
# Output directory: results/14_todo_analysis/task5_expanded_groups/
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
OUTPUT_BASE <- file.path(BASE_DIR, "results/14_todo_analysis/task5_expanded_groups")

FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
MIN_AVG_READS <- 10
MICROEXON_MAX <- 30
SMALL_MAX <- 50

EVENT_TYPES <- c("SE", "A3SS", "A5SS", "MXE", "RI")

# Group definitions (all event types, classified by exon_length like Venn diagram)
GROUPS <- list(
    group1_neg_microexon = list(
        label = "Microexon (0-30bp)",
        description = "Genes with microexon events (all types) significant in Neg (neg_only + both)",
        max_size = MICROEXON_MAX
    ),
    group2_neg_micro_small = list(
        label = "Microexon + Small (0-50bp)",
        description = "Genes with microexon or small exon events (all types) significant in Neg",
        max_size = SMALL_MAX
    ),
    group3_neg_all_sizes = list(
        label = "All sizes",
        description = "Genes with any event (all types, all sizes) significant in Neg",
        max_size = Inf
    )
)

# =============================================================================
# Utility Functions (same as 14a/14b)
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

# Generic coordinate extraction (same columns used for exon_length)
compute_exon_coords <- function(df, event_type) {
    switch(event_type,
        SE = data.frame(exon_start = df$exonStart_0base, exon_end = df$exonEnd),
        A3SS = data.frame(exon_start = df$longExonStart_0base, exon_end = df$longExonEnd),
        A5SS = data.frame(exon_start = df$longExonStart_0base, exon_end = df$longExonEnd),
        MXE = data.frame(exon_start = df$`1stExonStart_0base`, exon_end = df$`1stExonEnd`),
        RI = data.frame(exon_start = df$riExonStart_0base, exon_end = df$riExonEnd),
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
cat("14g: Expanded Gene Selection (3 Cumulative Groups)\n")
cat("  All event types included (matching Venn diagram)\n")
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
            # Add generic coordinate columns for output
            coords <- compute_exon_coords(df, et)
            df$exon_start <- coords$exon_start
            df$exon_end <- coords$exon_end
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

# All Neg significant events (the pool for all 3 groups) - ALL event types
cat(sprintf("\n  Neg significant events (all types): %d\n", nrow(neg_sig)))
cat(sprintf("    By event type: SE=%d, A3SS=%d, A5SS=%d, MXE=%d, RI=%d\n",
            sum(neg_sig$event_type == "SE"), sum(neg_sig$event_type == "A3SS"),
            sum(neg_sig$event_type == "A5SS"), sum(neg_sig$event_type == "MXE"),
            sum(neg_sig$event_type == "RI")))
cat(sprintf("    By size (all types):\n"))
cat(sprintf("      Microexon (0-30bp): %d\n", sum(neg_sig$exon_length <= MICROEXON_MAX)))
cat(sprintf("      Small (30-50bp):    %d\n", sum(neg_sig$exon_length > MICROEXON_MAX & neg_sig$exon_length <= SMALL_MAX)))
cat(sprintf("      Regular (>50bp):    %d\n", sum(neg_sig$exon_length > SMALL_MAX)))

# =============================================================================
# Process each group
# =============================================================================

for (group_name in names(GROUPS)) {
    group <- GROUPS[[group_name]]
    group_dir <- file.path(OUTPUT_BASE, group_name)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

    cat(sprintf("\n--- %s: %s ---\n", group_name, group$label))
    cat(sprintf("  %s\n", group$description))

    # Select events for this group (all Neg significant, within size cutoff)
    if (is.infinite(group$max_size)) {
        selection_events <- neg_sig
    } else {
        selection_events <- neg_sig %>% filter(exon_length <= group$max_size)
    }

    selection_genes <- unique(selection_events$GeneID)
    cat(sprintf("  Selection events (all types): %d\n", nrow(selection_events)))
    cat(sprintf("  Unique genes: %d\n", length(selection_genes)))

    if (length(selection_genes) == 0) {
        cat("  WARNING: No genes found for this group. Skipping.\n")
        next
    }

    # Event type breakdown of selection
    et_counts <- selection_events %>% count(event_type)
    cat("  Event type breakdown:\n")
    for (i in seq_len(nrow(et_counts))) {
        cat(sprintf("    %s: %d\n", et_counts$event_type[i], et_counts$n[i]))
    }

    # --- 1. Save selection events (the events driving gene selection) ---
    selection_detail <- selection_events %>%
        select(event_key, event_type, GeneID, geneSymbol, chr, strand,
               exon_start, exon_end, exon_length, size_category, dPSI, FDR) %>%
        arrange(geneSymbol, exon_start)
    write.csv(selection_detail,
              file.path(group_dir, "selection_events.csv"), row.names = FALSE)

    # --- 2. Gene list with summary ---
    gene_summary <- selection_detail %>%
        group_by(GeneID, geneSymbol) %>%
        summarize(
            n_selection_events = n(),
            event_types = paste(sort(unique(event_type)), collapse = "; "),
            size_categories = paste(sort(unique(as.character(size_category))), collapse = "; "),
            event_coords = paste(paste0(chr, ":", exon_start, "-", exon_end), collapse = "; "),
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
            n_microexon = sum(size_category == "Microexon (0-30bp)"),
            mean_abs_dPSI = round(mean(abs(dPSI), na.rm = TRUE), 3),
            .groups = "drop"
        )

    # Merge with selection gene info
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
cat("  CROSS-GROUP SUMMARY\n")
cat("========================================================\n\n")

for (group_name in names(GROUPS)) {
    group <- GROUPS[[group_name]]
    group_dir <- file.path(OUTPUT_BASE, group_name)
    gene_ids_file <- file.path(group_dir, "gene_ids.txt")
    if (file.exists(gene_ids_file)) {
        n_genes <- length(readLines(gene_ids_file))
        sel_events <- read.csv(file.path(group_dir, "selection_events.csv"))
        pos_summ <- read.csv(file.path(group_dir, "gene_pos_summary.csv"))
        cat(sprintf("  %-25s: %4d genes, %5d selection events, %d (%.0f%%) with sig Pos\n",
                    group$label, n_genes, nrow(sel_events),
                    sum(pos_summ$n_sig_pos > 0),
                    100 * sum(pos_summ$n_sig_pos > 0) / n_genes))
    }
}

cat(sprintf("\n  Output base: %s\n", OUTPUT_BASE))
cat("\n14g complete.\n")
