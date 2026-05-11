#!/usr/bin/env Rscript
# =============================================================================
# 14a_venn_diagrams.R
# Venn diagrams of splicing event overlap between Neg_vs_Parental and
# Pos_vs_Parental, stratified by exon size (Microexon/Small/Regular).
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(VennDiagram)
    library(grid)
})

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/14_todo_analysis/task1_venn_diagrams")

FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
MIN_AVG_READS <- 10
MICROEXON_MAX <- 30
SMALL_EXON_MAX <- 50

EVENT_TYPES <- c("SE", "A3SS", "A5SS", "MXE", "RI")

size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")

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
    # rMATS has duplicate 'ID' column (positions 1 and 12); deduplicate
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
            lengths <= SMALL_EXON_MAX ~ "Small (30-50bp)",
            TRUE ~ "Regular (>50bp)"
        ),
        levels = c("Microexon (0-30bp)", "Small (30-50bp)", "Regular (>50bp)")
    )
}

# =============================================================================
# Load all events for both comparisons
# =============================================================================
cat("========================================\n")
cat("14a: Venn Diagrams — Splicing Event Overlap\n")
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

cat(sprintf("  Neg_vs_Parental: %d events loaded, %d significant\n",
            nrow(neg_all), sum(neg_all$significant)))
cat(sprintf("  Pos_vs_Parental: %d events loaded, %d significant\n",
            nrow(pos_all), sum(pos_all$significant)))

neg_sig <- neg_all %>% filter(significant)
pos_sig <- pos_all %>% filter(significant)

# =============================================================================
# Draw Venn diagrams
# =============================================================================
draw_size_venn <- function(size_label, neg_keys, pos_keys, filename) {
    neg_only <- length(setdiff(neg_keys, pos_keys))
    pos_only <- length(setdiff(pos_keys, neg_keys))
    both <- length(intersect(neg_keys, pos_keys))

    if (neg_only + pos_only + both == 0) {
        cat(sprintf("  Skipping %s Venn: no events\n", size_label))
        return(data.frame(size = size_label, neg_only = 0, pos_only = 0, both = 0))
    }

    cat(sprintf("  %s: Neg-only=%d, Pos-only=%d, Both=%d\n",
                size_label, neg_only, pos_only, both))

    pdf(file.path(OUTPUT_DIR, filename), width = 6, height = 5)
    grid.newpage()
    draw.pairwise.venn(
        area1 = neg_only + both,
        area2 = pos_only + both,
        cross.area = both,
        category = c("Neg vs Parental", "Pos vs Parental"),
        fill = c("#377EB8", "#4DAF4A"),
        alpha = 0.5,
        cat.pos = c(-30, 30),
        cat.dist = 0.05,
        cat.cex = 1.1,
        cex = 1.3,
        fontfamily = "sans",
        cat.fontfamily = "sans"
    )
    grid.text(size_label, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()

    data.frame(size = size_label, neg_only = neg_only, pos_only = pos_only, both = both)
}

# Per-size Venns
summary_rows <- list()
for (sz in levels(neg_sig$size_category)) {
    neg_keys <- neg_sig %>% filter(size_category == sz) %>% pull(event_key)
    pos_keys <- pos_sig %>% filter(size_category == sz) %>% pull(event_key)
    safe_name <- gsub("[^a-zA-Z0-9]", "_", tolower(gsub(" \\(.*", "", sz)))
    summary_rows[[sz]] <- draw_size_venn(sz, neg_keys, pos_keys,
                                          paste0("venn_", safe_name, ".pdf"))
}

# All events combined
summary_rows[["All"]] <- draw_size_venn("All Events",
                                         neg_sig$event_key, pos_sig$event_key,
                                         "venn_all_events.pdf")

# Summary CSV
summary_df <- bind_rows(summary_rows)
summary_df$total_neg <- summary_df$neg_only + summary_df$both
summary_df$total_pos <- summary_df$pos_only + summary_df$both
summary_df$jaccard <- round(summary_df$both /
    (summary_df$neg_only + summary_df$pos_only + summary_df$both), 3)
write.csv(summary_df, file.path(OUTPUT_DIR, "venn_overlap_summary.csv"), row.names = FALSE)
cat("\n  Overlap summary:\n")
print(summary_df)

# =============================================================================
# Combined multi-panel figure
# =============================================================================
pdf(file.path(OUTPUT_DIR, "venn_combined.pdf"), width = 15, height = 5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
for (i in seq_along(levels(neg_sig$size_category))) {
    sz <- levels(neg_sig$size_category)[i]
    neg_keys <- neg_sig %>% filter(size_category == sz) %>% pull(event_key)
    pos_keys <- pos_sig %>% filter(size_category == sz) %>% pull(event_key)
    neg_only <- length(setdiff(neg_keys, pos_keys))
    pos_only <- length(setdiff(pos_keys, neg_keys))
    both <- length(intersect(neg_keys, pos_keys))
    if (neg_only + pos_only + both == 0) next

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
    venn <- draw.pairwise.venn(
        area1 = neg_only + both,
        area2 = pos_only + both,
        cross.area = both,
        category = c("Neg vs Par", "Pos vs Par"),
        fill = c("#377EB8", "#4DAF4A"),
        alpha = 0.5,
        cat.pos = c(-30, 30),
        cat.dist = 0.05,
        cat.cex = 0.9,
        cex = 1.1,
        fontfamily = "sans",
        cat.fontfamily = "sans",
        ind = FALSE
    )
    grid.draw(venn)
    grid.text(sz, y = 0.95, gp = gpar(fontsize = 11, fontface = "bold"))
    popViewport()
}
dev.off()

# =============================================================================
# Export event lists
# =============================================================================
neg_keys_all <- neg_sig$event_key
pos_keys_all <- pos_sig$event_key
overlap_keys <- intersect(neg_keys_all, pos_keys_all)
neg_only_keys <- setdiff(neg_keys_all, pos_keys_all)
pos_only_keys <- setdiff(pos_keys_all, neg_keys_all)

if (length(overlap_keys) > 0) {
    overlap_neg <- neg_sig %>% filter(event_key %in% overlap_keys) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR) %>%
        rename(dPSI_neg = dPSI, FDR_neg = FDR)
    overlap_pos <- pos_sig %>% filter(event_key %in% overlap_keys) %>%
        select(event_key, dPSI, FDR) %>%
        rename(dPSI_pos = dPSI, FDR_pos = FDR)
    overlap_df <- left_join(overlap_neg, overlap_pos, by = "event_key")
    write.csv(overlap_df, file.path(OUTPUT_DIR, "venn_overlap_events.csv"), row.names = FALSE)
}
if (length(neg_only_keys) > 0) {
    neg_only_df <- neg_sig %>% filter(event_key %in% neg_only_keys) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR)
    write.csv(neg_only_df, file.path(OUTPUT_DIR, "venn_neg_only_events.csv"), row.names = FALSE)
}
if (length(pos_only_keys) > 0) {
    pos_only_df <- pos_sig %>% filter(event_key %in% pos_only_keys) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR)
    write.csv(pos_only_df, file.path(OUTPUT_DIR, "venn_pos_only_events.csv"), row.names = FALSE)
}

cat("\n14a complete.\n")
