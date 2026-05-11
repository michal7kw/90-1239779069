#!/usr/bin/env Rscript
# =============================================================================
# 14c_neg_exclusive_composition.R
# Stage 2: For the neg-exclusive microexon genes identified by 14b, load ALL
# their splicing events from Pos_vs_Parental and characterize the composition
# (event types, size distribution, dPSI, Neg-vs-Pos scatter).
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
INPUT_DIR <- file.path(BASE_DIR, "results/14_todo_analysis/task2_neg_exclusive_microexon")
OUTPUT_DIR <- INPUT_DIR  # same directory, different files

FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
MIN_AVG_READS <- 10
MICROEXON_MAX <- 30

EVENT_TYPES <- c("SE", "A3SS", "A5SS", "MXE", "RI")

event_colors <- c(SE = "#1f77b4", A5SS = "#ff7f0e", A3SS = "#2ca02c",
                  MXE = "#d62728", RI = "#9467bd")
size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")

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
# Load gene list from 14b
# =============================================================================
cat("========================================\n")
cat("14c: Neg-Exclusive Microexon Gene Composition\n")
cat("========================================\n")

gene_id_file <- file.path(INPUT_DIR, "neg_exclusive_gene_ids.txt")
if (!file.exists(gene_id_file)) {
    stop("Input not found: ", gene_id_file,
         "\n  Run 14b_neg_exclusive_genes.R first.")
}
neg_excl_genes <- readLines(gene_id_file)
cat(sprintf("  Loaded %d neg-exclusive microexon gene IDs from 14b\n", length(neg_excl_genes)))

# Also load neg events for these genes (produced by 14b, for scatter plot)
neg_events_file <- file.path(INPUT_DIR, "neg_all_events_for_genes.csv")
if (file.exists(neg_events_file)) {
    neg_events <- read.csv(neg_events_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded %d Neg events for these genes\n", nrow(neg_events)))
} else {
    neg_events <- NULL
    cat("  Warning: neg_all_events_for_genes.csv not found, scatter plot will be skipped\n")
}

# Load gene info for summary
gene_info_file <- file.path(INPUT_DIR, "neg_exclusive_microexon_genes.csv")
gene_info <- if (file.exists(gene_info_file)) {
    read.csv(gene_info_file, stringsAsFactors = FALSE)
} else NULL

# =============================================================================
# Load ALL Pos_vs_Parental events for these genes
# =============================================================================
cat("\n  Loading Pos_vs_Parental events for neg-exclusive genes...\n")

pos_gene_events <- list()
for (et in EVENT_TYPES) {
    df <- load_rmats_results("Pos_vs_Parental", et)
    if (!is.null(df) && nrow(df) > 0) {
        df <- df %>% filter(GeneID %in% neg_excl_genes)
        if (nrow(df) > 0) {
            df$event_key <- create_event_key(df, et)
            df$exon_length <- compute_exon_length(df, et)
            df$size_category <- classify_size(df$exon_length)
            df$significant <- df$FDR < FDR_THRESHOLD & abs(df$IncLevelDifference) > DPSI_THRESHOLD
            df$dPSI <- df$IncLevelDifference
            pos_gene_events[[et]] <- df
        }
    }
}
pos_gene_events <- bind_rows(pos_gene_events)

cat(sprintf("  Total Pos events for these genes: %d (across %d genes)\n",
            nrow(pos_gene_events), length(unique(pos_gene_events$GeneID))))

# Save all Pos events
pos_export <- pos_gene_events %>%
    select(event_key, event_type, GeneID, geneSymbol, exon_length,
           size_category, dPSI, FDR, significant)
write.csv(pos_export,
          file.path(OUTPUT_DIR, "neg_exclusive_genes_pos_events.csv"), row.names = FALSE)

# =============================================================================
# Event type distribution
# =============================================================================
cat("\n  Generating plots...\n")

et_counts <- pos_gene_events %>%
    count(event_type) %>%
    mutate(pct = round(100 * n / sum(n), 1))

p_et <- ggplot(et_counts, aes(x = reorder(event_type, -n), y = n, fill = event_type)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(n, "\n(", pct, "%)")), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = event_colors) +
    labs(title = "Event Type Distribution in Pos_vs_Parental",
         subtitle = paste("Genes with Neg-exclusive microexons (n =",
                          length(neg_excl_genes), "genes)"),
         x = "Event Type", y = "Number of Events") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none") +
    ylim(0, max(et_counts$n) * 1.15)
ggsave(file.path(OUTPUT_DIR, "event_type_distribution.pdf"), p_et, width = 7, height = 5)

# =============================================================================
# Size distribution (SE events only)
# =============================================================================
pos_se <- pos_gene_events %>% filter(event_type == "SE")
if (nrow(pos_se) > 0) {
    sz_counts <- pos_se %>%
        count(size_category) %>%
        mutate(pct = round(100 * n / sum(n), 1))

    p_sz <- ggplot(sz_counts, aes(x = size_category, y = n, fill = size_category)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = paste0(n, "\n(", pct, "%)")), vjust = -0.3, size = 3.5) +
        scale_fill_manual(values = size_colors) +
        labs(title = "SE Event Size Distribution in Pos_vs_Parental",
             subtitle = "SE events in genes with Neg-exclusive microexons",
             x = "Size Category", y = "Number of SE Events") +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none") +
        ylim(0, max(sz_counts$n) * 1.15)
    ggsave(file.path(OUTPUT_DIR, "size_distribution.pdf"), p_sz, width = 7, height = 5)
}

# =============================================================================
# dPSI distribution colored by significance
# =============================================================================
p_dpsi <- ggplot(pos_gene_events, aes(x = dPSI, fill = significant)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = c(-DPSI_THRESHOLD, DPSI_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#95a5a6"),
                      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
                      name = "FDR < 0.05 &\n|dPSI| > 0.1") +
    labs(title = "dPSI Distribution in Pos_vs_Parental",
         subtitle = "Events in genes with Neg-exclusive microexons",
         x = "IncLevelDifference (dPSI)", y = "Density") +
    theme_minimal(base_size = 13)
ggsave(file.path(OUTPUT_DIR, "dpsi_distribution.pdf"), p_dpsi, width = 7, height = 5)

# =============================================================================
# dPSI scatter: Neg vs Pos for matched events
# =============================================================================
if (!is.null(neg_events)) {
    neg_for_join <- neg_events %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR, significant) %>%
        rename(dPSI_neg = dPSI, FDR_neg = FDR, sig_neg = significant)

    pos_for_join <- pos_gene_events %>%
        select(event_key, dPSI, FDR, significant) %>%
        rename(dPSI_pos = dPSI, FDR_pos = FDR, sig_pos = significant)

    matched <- inner_join(neg_for_join, pos_for_join, by = "event_key")

    if (nrow(matched) > 0) {
        matched$sig_label <- case_when(
            matched$sig_neg & matched$sig_pos ~ "Both",
            matched$sig_neg ~ "Neg only",
            matched$sig_pos ~ "Pos only",
            TRUE ~ "Neither"
        )

        corr_val <- cor(matched$dPSI_neg, matched$dPSI_pos, use = "complete.obs")

        p_scatter <- ggplot(matched, aes(x = dPSI_neg, y = dPSI_pos, color = sig_label)) +
            geom_point(alpha = 0.6, size = 1.5) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            geom_hline(yintercept = 0, color = "grey80") +
            geom_vline(xintercept = 0, color = "grey80") +
            scale_color_manual(values = c("Both" = "#e74c3c", "Neg only" = "#377EB8",
                                           "Pos only" = "#4DAF4A", "Neither" = "#bdc3c7"),
                               name = "Significant in") +
            labs(title = "dPSI Correlation: Neg vs Pos for Neg-Exclusive Microexon Genes",
                 subtitle = sprintf("n = %d matched events, r = %.3f", nrow(matched), corr_val),
                 x = "dPSI (Neg_vs_Parental)", y = "dPSI (Pos_vs_Parental)") +
            theme_minimal(base_size = 13) +
            coord_fixed()
        ggsave(file.path(OUTPUT_DIR, "neg_vs_pos_scatter.pdf"), p_scatter, width = 7, height = 7)

        cat(sprintf("  Scatter: %d matched events, r = %.3f\n", nrow(matched), corr_val))
    }
}

# =============================================================================
# Per-gene summary
# =============================================================================
gene_summary_pos <- pos_gene_events %>%
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

# Merge with gene info from 14b if available
if (!is.null(gene_info)) {
    gene_summary_pos <- gene_summary_pos %>%
        left_join(
            gene_info %>% select(GeneID, n_neg_excl_microexons, mean_dPSI),
            by = "GeneID"
        ) %>%
        rename(neg_excl_mean_dPSI = mean_dPSI) %>%
        arrange(desc(n_neg_excl_microexons))
}

write.csv(gene_summary_pos,
          file.path(OUTPUT_DIR, "neg_exclusive_genes_pos_summary.csv"), row.names = FALSE)

cat(sprintf("\n  Gene summary: %d genes\n", nrow(gene_summary_pos)))
cat(sprintf("  Mean Pos events per gene: %.1f\n", mean(gene_summary_pos$n_events_pos)))
cat(sprintf("  Mean Pos significant per gene: %.1f\n", mean(gene_summary_pos$n_sig_pos)))

cat("\n14c complete.\n")
