#!/usr/bin/env Rscript
# =============================================================================
# 14i_se_negonly_expanded.R
#
# SE-ONLY version of 14i_negonly_expanded.R
#
# Neg-only (exclusive) gene selection expanded to larger exon sizes.
# Uses ONLY neg_only SE events (significant in Neg but NOT in Pos).
#
# Two groups:
#   Group 2b: Neg-only microexon + small exon SE (0-50bp)
#   Group 3b: Neg-only all SE sizes
#
# Combined script: data preparation + 14f-style characterization.
#
# Output: results/14_todo_analysis/task5_expanded_groups_se/group{2b,3b}_*/
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(gridExtra)
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

# Color palettes
event_colors <- c(SE = "#1f77b4", A5SS = "#ff7f0e", A3SS = "#2ca02c",
                  MXE = "#d62728", RI = "#9467bd")
size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")
direction_colors <- c("Inclusion (dPSI > 0)" = "#2ecc71",
                      "Exclusion (dPSI < 0)" = "#e74c3c")
condition_colors <- c("Neg_vs_Parental" = "#3498db",
                      "Pos_vs_Parental" = "#e67e22")

# Group definitions (SE events only)
GROUPS <- list(
    group2b_negonly_micro_small = list(
        label = "Neg-Only Micro+Small SE (0-50bp)",
        description = "Genes with neg-only (not shared) microexon or small exon SE events",
        max_size = SMALL_MAX
    ),
    group3b_negonly_all_sizes = list(
        label = "Neg-Only All SE Sizes",
        description = "Genes with neg-only (not shared) SE events of any size",
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

add_direction <- function(df) {
    df %>%
        mutate(direction = ifelse(dPSI > 0,
                                  "Inclusion (dPSI > 0)",
                                  "Exclusion (dPSI < 0)"))
}

# =============================================================================
# PART 1: Load data and identify neg-only SE events
# =============================================================================
cat("========================================================\n")
cat("14i_se: Neg-Only Expanded Gene Selection + Characterization\n")
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

# Get neg-only SE events (significant in Neg SE, NOT significant in Pos SE)
neg_sig_se <- neg_sig %>% filter(event_type == "SE")
pos_sig_se <- pos_sig %>% filter(event_type == "SE")
neg_only_se <- neg_sig_se %>% filter(!(event_key %in% pos_sig_se$event_key))

cat(sprintf("\n  Neg significant SE events: %d\n", nrow(neg_sig_se)))
cat(sprintf("  Pos significant SE events: %d\n", nrow(pos_sig_se)))
cat(sprintf("  Neg-only SE events (not shared with Pos): %d\n", nrow(neg_only_se)))
cat(sprintf("    Microexon (0-30bp): %d\n", sum(neg_only_se$exon_length <= MICROEXON_MAX)))
cat(sprintf("    Small (30-50bp):    %d\n", sum(neg_only_se$exon_length > MICROEXON_MAX & neg_only_se$exon_length <= SMALL_MAX)))
cat(sprintf("    Regular (>50bp):    %d\n", sum(neg_only_se$exon_length > SMALL_MAX)))

# =============================================================================
# PART 2: Process each group (data prep + characterization)
# =============================================================================

for (group_name in names(GROUPS)) {
    group <- GROUPS[[group_name]]
    group_dir <- file.path(OUTPUT_BASE, group_name)
    output_dir <- file.path(group_dir, "pos_significant_deep")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    cat(sprintf("\n\n========== %s: %s ==========\n", group_name, group$label))
    cat(sprintf("  %s\n", group$description))

    # --- Select neg-only SE events within size cutoff ---
    if (is.infinite(group$max_size)) {
        selection_events <- neg_only_se
    } else {
        selection_events <- neg_only_se %>% filter(exon_length <= group$max_size)
    }

    selection_genes <- unique(selection_events$GeneID)
    cat(sprintf("  Selection SE events: %d\n", nrow(selection_events)))
    cat(sprintf("  Unique genes: %d\n", length(selection_genes)))

    if (length(selection_genes) == 0) {
        cat("  WARNING: No genes found. Skipping.\n")
        next
    }

    # --- Save selection data ---
    selection_detail <- selection_events %>%
        select(event_key, event_type, GeneID, geneSymbol, chr, strand,
               exonStart_0base, exonEnd, exon_length, size_category, dPSI, FDR) %>%
        arrange(geneSymbol, exonStart_0base)
    write.csv(selection_detail, file.path(group_dir, "selection_events.csv"), row.names = FALSE)

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
    write.csv(gene_summary, file.path(group_dir, "gene_list.csv"), row.names = FALSE)
    writeLines(selection_genes, file.path(group_dir, "gene_ids.txt"))

    # All Neg events for these genes (all event types)
    neg_events_for_genes <- neg_all %>%
        filter(GeneID %in% selection_genes) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR, significant)
    write.csv(neg_events_for_genes, file.path(group_dir, "neg_all_events_for_genes.csv"), row.names = FALSE)

    # All Pos events for these genes (all event types)
    pos_events_for_genes <- pos_all %>%
        filter(GeneID %in% selection_genes) %>%
        select(event_key, event_type, GeneID, geneSymbol, exon_length,
               size_category, dPSI, FDR, significant)
    write.csv(pos_events_for_genes, file.path(group_dir, "pos_all_events_for_genes.csv"), row.names = FALSE)

    cat(sprintf("  Neg events (all types) for these genes: %d\n", nrow(neg_events_for_genes)))
    cat(sprintf("  Pos events (all types) for these genes: %d\n", nrow(pos_events_for_genes)))

    # Per-gene Pos summary
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
        ) %>%
        left_join(
            gene_summary %>% select(GeneID, n_selection_events, mean_dPSI),
            by = "GeneID"
        ) %>%
        rename(selection_mean_dPSI = mean_dPSI) %>%
        arrange(desc(n_selection_events))
    write.csv(gene_pos_summary, file.path(group_dir, "gene_pos_summary.csv"), row.names = FALSE)

    # ========================================================================
    # 14f-style characterization
    # ========================================================================
    pos_sig_grp <- pos_events_for_genes %>% filter(significant == TRUE) %>% add_direction()
    neg_sig_grp <- neg_events_for_genes %>% filter(significant == TRUE) %>% add_direction()

    total_genes <- nrow(gene_pos_summary)
    genes_with_sig_pos <- sum(gene_pos_summary$n_sig_pos > 0)
    genes_without_sig_pos <- total_genes - genes_with_sig_pos

    cat(sprintf("  Pos significant events: %d (in %d genes)\n",
                nrow(pos_sig_grp), n_distinct(pos_sig_grp$geneSymbol)))
    cat(sprintf("  Neg significant events: %d (in %d genes)\n",
                nrow(neg_sig_grp), n_distinct(neg_sig_grp$geneSymbol)))
    cat(sprintf("  Genes WITH significant Pos events:    %d (%.1f%%)\n",
                genes_with_sig_pos, 100 * genes_with_sig_pos / total_genes))

    # Gene classification
    gene_class <- gene_pos_summary %>%
        mutate(sig_pos_bin = case_when(
            n_sig_pos == 0 ~ "0", n_sig_pos == 1 ~ "1",
            n_sig_pos == 2 ~ "2", TRUE ~ "3+"
        )) %>%
        mutate(sig_pos_bin = factor(sig_pos_bin, levels = c("0", "1", "2", "3+")))

    gene_class_counts <- gene_class %>%
        count(sig_pos_bin, .drop = FALSE) %>% rename(n_genes = n)

    # Comparison data frames
    neg_type_counts <- neg_sig_grp %>% count(event_type, name = "count") %>% mutate(condition = "Neg_vs_Parental")
    pos_type_counts <- pos_sig_grp %>% count(event_type, name = "count") %>% mutate(condition = "Pos_vs_Parental")
    type_comparison <- bind_rows(neg_type_counts, pos_type_counts) %>%
        complete(event_type = names(event_colors), condition, fill = list(count = 0)) %>%
        mutate(event_type = factor(event_type, levels = names(event_colors)))

    neg_size_counts <- neg_sig_grp %>% filter(event_type == "SE") %>% count(size_category, name = "count") %>% mutate(condition = "Neg_vs_Parental")
    pos_size_counts <- pos_sig_grp %>% filter(event_type == "SE") %>% count(size_category, name = "count") %>% mutate(condition = "Pos_vs_Parental")
    size_comparison <- bind_rows(neg_size_counts, pos_size_counts) %>%
        complete(size_category = names(size_colors), condition, fill = list(count = 0)) %>%
        mutate(size_category = factor(size_category, levels = names(size_colors)))

    neg_dir_counts <- neg_sig_grp %>% count(direction, name = "count") %>% mutate(condition = "Neg_vs_Parental")
    pos_dir_counts <- pos_sig_grp %>% count(direction, name = "count") %>% mutate(condition = "Pos_vs_Parental")
    dir_comparison <- bind_rows(neg_dir_counts, pos_dir_counts) %>%
        complete(direction = names(direction_colors), condition, fill = list(count = 0))

    # Per-gene profile table
    neg_profile <- neg_sig_grp %>%
        group_by(geneSymbol) %>%
        summarise(neg_n_sig = n(), neg_event_types = paste(sort(unique(event_type)), collapse = ","),
                  neg_n_SE = sum(event_type == "SE"),
                  neg_n_microexon = sum(size_category == "Microexon (0-30bp)", na.rm = TRUE),
                  neg_n_inclusion = sum(dPSI > 0), neg_n_exclusion = sum(dPSI < 0),
                  neg_mean_dPSI = round(mean(dPSI), 4),
                  neg_mean_abs_dPSI = round(mean(abs(dPSI)), 4), .groups = "drop")

    pos_profile <- pos_sig_grp %>%
        group_by(geneSymbol) %>%
        summarise(pos_n_sig = n(), pos_event_types = paste(sort(unique(event_type)), collapse = ","),
                  pos_n_SE = sum(event_type == "SE"),
                  pos_n_microexon = sum(size_category == "Microexon (0-30bp)", na.rm = TRUE),
                  pos_n_inclusion = sum(dPSI > 0), pos_n_exclusion = sum(dPSI < 0),
                  pos_mean_dPSI = round(mean(dPSI), 4),
                  pos_mean_abs_dPSI = round(mean(abs(dPSI)), 4), .groups = "drop")

    gene_profile <- gene_pos_summary %>%
        select(GeneID, geneSymbol) %>%
        left_join(neg_profile, by = "geneSymbol") %>%
        left_join(pos_profile, by = "geneSymbol") %>%
        mutate(across(starts_with("neg_n"), ~replace_na(.x, 0L)),
               across(starts_with("pos_n"), ~replace_na(.x, 0L)),
               neg_event_types = replace_na(neg_event_types, ""),
               pos_event_types = replace_na(pos_event_types, ""))

    # --- Generate 6 plots ---
    cat("  Generating plots...\n")

    p1 <- ggplot(gene_class_counts, aes(x = sig_pos_bin, y = n_genes)) +
        geom_col(fill = "#e67e22", width = 0.6) +
        geom_text(aes(label = n_genes), vjust = -0.5, size = 4.5) +
        labs(title = paste0(group$label, " Genes:\nNumber of Significant Pos Events"),
             x = "Number of significant Pos_vs_Parental events", y = "Number of genes") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"))
    ggsave(file.path(output_dir, "gene_classification_bar.pdf"), p1, width = 6, height = 5)

    p2 <- ggplot(type_comparison, aes(x = event_type, y = count, fill = condition)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(aes(label = count), position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = condition_colors, name = "Comparison") +
        labs(title = paste0("Significant Events by Type: Neg vs Pos\n(", group$label, " genes)"),
             x = "Event type", y = "Number of significant events") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top", panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"))
    ggsave(file.path(output_dir, "event_type_comparison.pdf"), p2, width = 7, height = 5)

    p3 <- ggplot(size_comparison, aes(x = size_category, y = count, fill = condition)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(aes(label = count), position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = condition_colors, name = "Comparison") +
        labs(title = paste0("Significant SE Events by Size: Neg vs Pos\n(", group$label, " genes)"),
             subtitle = "Skipped exon events only",
             x = "Size category", y = "Number of significant SE events") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top", panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(angle = 15, hjust = 1))
    ggsave(file.path(output_dir, "size_category_comparison.pdf"), p3, width = 7, height = 5)

    p4 <- ggplot(dir_comparison, aes(x = direction, y = count, fill = condition)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(aes(label = count), position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = condition_colors, name = "Comparison") +
        labs(title = paste0("Significant Events by Direction: Neg vs Pos\n(", group$label, " genes)"),
             x = "Direction", y = "Number of significant events") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top", panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"))
    ggsave(file.path(output_dir, "direction_comparison.pdf"), p4, width = 7, height = 5)

    # Heatmap
    genes_with_pos <- gene_profile %>% filter(pos_n_sig > 0)
    if (nrow(genes_with_pos) > 0) {
        max_hm <- 50
        if (nrow(genes_with_pos) > max_hm) {
            genes_with_pos <- genes_with_pos %>%
                mutate(total = neg_n_sig + pos_n_sig) %>%
                arrange(desc(total)) %>% head(max_hm)
            hm_sub <- sprintf("Top %d genes by total sig events (of %d with sig Pos)",
                              max_hm, sum(gene_profile$pos_n_sig > 0))
        } else {
            hm_sub <- sprintf("Genes with at least 1 significant Pos event (n = %d)",
                              nrow(genes_with_pos))
        }

        hm_data <- bind_rows(
            neg_sig_grp %>% filter(geneSymbol %in% genes_with_pos$geneSymbol) %>% mutate(condition = "Neg_vs_Parental"),
            pos_sig_grp %>% filter(geneSymbol %in% genes_with_pos$geneSymbol) %>% mutate(condition = "Pos_vs_Parental")
        ) %>%
            group_by(geneSymbol, condition, event_type) %>%
            summarise(count = n(), .groups = "drop") %>%
            complete(geneSymbol = unique(genes_with_pos$geneSymbol),
                     condition = c("Neg_vs_Parental", "Pos_vs_Parental"),
                     event_type = names(event_colors), fill = list(count = 0))

        gene_order <- genes_with_pos %>%
            mutate(total = neg_n_sig + pos_n_sig) %>%
            arrange(desc(total)) %>% pull(geneSymbol)
        hm_data <- hm_data %>%
            mutate(geneSymbol = factor(geneSymbol, levels = rev(gene_order)),
                   event_type = factor(event_type, levels = names(event_colors)))

        p5 <- ggplot(hm_data, aes(x = event_type, y = geneSymbol, fill = count)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = ifelse(count > 0, count, "")), size = 3) +
            facet_wrap(~condition) +
            scale_fill_gradient(low = "white", high = "#2c3e50", name = "# Events") +
            labs(title = paste0("Significant Events per Gene: Neg vs Pos\n(", group$label, ")"),
                 subtitle = hm_sub, x = "Event type", y = "") +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(size = 12, face = "bold"),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  strip.text = element_text(face = "bold", size = 11),
                  panel.grid = element_blank())
        ggsave(file.path(output_dir, "gene_event_heatmap.pdf"),
               p5, width = 10, height = max(6, min(nrow(genes_with_pos) * 0.3 + 3, 30)))
    }

    combined <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                             top = grid::textGrob(
                                 paste0(group$label, " Genes: Pos Event Characterization"),
                                 gp = grid::gpar(fontsize = 14, fontface = "bold")))
    ggsave(file.path(output_dir, "combined_panel.pdf"), combined, width = 14, height = 12)
    cat("  Saved 6 PDFs\n")

    # Save CSVs
    write.csv(pos_sig_grp, file.path(output_dir, "pos_significant_events_detail.csv"), row.names = FALSE)
    write.csv(gene_profile, file.path(output_dir, "gene_splicing_profile_comparison.csv"), row.names = FALSE)

    summary_stats <- data.frame(
        metric = c("group_name", "group_label", "total_genes",
                    "genes_with_sig_pos_events", "genes_without_sig_pos_events", "pct_genes_with_sig_pos",
                    "pos_total_sig_events",
                    "pos_sig_SE", "pos_sig_A3SS", "pos_sig_A5SS", "pos_sig_MXE", "pos_sig_RI",
                    "pos_sig_microexon", "pos_sig_small", "pos_sig_regular",
                    "pos_sig_inclusion", "pos_sig_exclusion",
                    "neg_total_sig_events",
                    "neg_sig_SE", "neg_sig_A3SS", "neg_sig_A5SS", "neg_sig_MXE", "neg_sig_RI",
                    "neg_sig_microexon", "neg_sig_small", "neg_sig_regular",
                    "neg_sig_inclusion", "neg_sig_exclusion"),
        value = c(group_name, group$label, total_genes,
                  genes_with_sig_pos, genes_without_sig_pos,
                  round(100 * genes_with_sig_pos / total_genes, 1),
                  nrow(pos_sig_grp),
                  sum(pos_sig_grp$event_type == "SE"), sum(pos_sig_grp$event_type == "A3SS"),
                  sum(pos_sig_grp$event_type == "A5SS"), sum(pos_sig_grp$event_type == "MXE"),
                  sum(pos_sig_grp$event_type == "RI"),
                  sum(pos_sig_grp$event_type == "SE" & pos_sig_grp$size_category == "Microexon (0-30bp)"),
                  sum(pos_sig_grp$event_type == "SE" & pos_sig_grp$size_category == "Small (30-50bp)"),
                  sum(pos_sig_grp$event_type == "SE" & pos_sig_grp$size_category == "Regular (>50bp)"),
                  sum(pos_sig_grp$dPSI > 0), sum(pos_sig_grp$dPSI < 0),
                  nrow(neg_sig_grp),
                  sum(neg_sig_grp$event_type == "SE"), sum(neg_sig_grp$event_type == "A3SS"),
                  sum(neg_sig_grp$event_type == "A5SS"), sum(neg_sig_grp$event_type == "MXE"),
                  sum(neg_sig_grp$event_type == "RI"),
                  sum(neg_sig_grp$event_type == "SE" & neg_sig_grp$size_category == "Microexon (0-30bp)"),
                  sum(neg_sig_grp$event_type == "SE" & neg_sig_grp$size_category == "Small (30-50bp)"),
                  sum(neg_sig_grp$event_type == "SE" & neg_sig_grp$size_category == "Regular (>50bp)"),
                  sum(neg_sig_grp$dPSI > 0), sum(neg_sig_grp$dPSI < 0)),
        stringsAsFactors = FALSE
    )
    write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
    cat("  Saved 3 CSVs\n")

    # Console summary
    cat(sprintf("\n  --- %s SUMMARY ---\n", group$label))
    cat(sprintf("  %d / %d genes (%.0f%%) have significant Pos events\n",
                genes_with_sig_pos, total_genes, 100 * genes_with_sig_pos / total_genes))

    cat("  Event type breakdown (Neg / Pos sig events):\n")
    for (et in names(event_colors)) {
        cat(sprintf("    %-5s: %3d / %3d\n", et,
                    sum(neg_sig_grp$event_type == et), sum(pos_sig_grp$event_type == et)))
    }
    cat("  Size breakdown - SE events only (Neg / Pos):\n")
    for (sc in names(size_colors)) {
        cat(sprintf("    %-20s: %3d / %3d\n", sc,
                    sum(neg_sig_grp$event_type == "SE" & neg_sig_grp$size_category == sc),
                    sum(pos_sig_grp$event_type == "SE" & pos_sig_grp$size_category == sc)))
    }
    cat("  Direction (Neg / Pos):\n")
    cat(sprintf("    Inclusion:  %3d / %3d\n", sum(neg_sig_grp$dPSI > 0), sum(pos_sig_grp$dPSI > 0)))
    cat(sprintf("    Exclusion:  %3d / %3d\n", sum(neg_sig_grp$dPSI < 0), sum(pos_sig_grp$dPSI < 0)))

    cat(sprintf("\n  Output: %s\n", output_dir))
}

cat("\n\n14i_se complete.\n")
