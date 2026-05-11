#!/usr/bin/env Rscript
# =============================================================================
# 14h_expanded_pos_characterization.R
#
# 14f-style deep characterization of significant Pos_vs_Parental events
# for each of the 3 expanded gene groups prepared by 14g.
#
# For each group, produces:
#   6 PDFs: gene_classification_bar, event_type_comparison, size_category_comparison,
#           direction_comparison, gene_event_heatmap, combined_panel
#   3 CSVs: pos_significant_events_detail, gene_splicing_profile_comparison,
#           summary_statistics
#
# Input:  results/14_todo_analysis/task5_expanded_groups/{group}/
# Output: results/14_todo_analysis/task5_expanded_groups/{group}/pos_significant_deep/
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
INPUT_BASE <- file.path(BASE_DIR, "results/14_todo_analysis/task5_expanded_groups")

# Color palettes (matching 14c/14f conventions)
event_colors <- c(SE = "#1f77b4", A5SS = "#ff7f0e", A3SS = "#2ca02c",
                  MXE = "#d62728", RI = "#9467bd")

size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")

direction_colors <- c("Inclusion (dPSI > 0)" = "#2ecc71",
                      "Exclusion (dPSI < 0)" = "#e74c3c")

condition_colors <- c("Neg_vs_Parental" = "#3498db",
                      "Pos_vs_Parental" = "#e67e22")

# Group metadata
GROUP_LABELS <- list(
    group1_neg_microexon = "Microexon (0-30bp)",
    group2_neg_micro_small = "Microexon + Small (0-50bp)",
    group3_neg_all_sizes = "All sizes"
)

# =============================================================================
# Helper: add direction column
# =============================================================================
add_direction <- function(df) {
    df %>%
        mutate(direction = ifelse(dPSI > 0,
                                  "Inclusion (dPSI > 0)",
                                  "Exclusion (dPSI < 0)"))
}

# =============================================================================
# Main: Process each group
# =============================================================================
cat("========================================================\n")
cat("14h: Expanded Pos Characterization (3 Groups)\n")
cat("========================================================\n")

groups <- names(GROUP_LABELS)

for (group_name in groups) {
    group_label <- GROUP_LABELS[[group_name]]
    INPUT_DIR <- file.path(INPUT_BASE, group_name)
    OUTPUT_DIR <- file.path(INPUT_DIR, "pos_significant_deep")
    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

    cat(sprintf("\n\n=== Processing: %s (%s) ===\n", group_name, group_label))

    # ========================================================================
    # Load data
    # ========================================================================
    required_files <- c("pos_all_events_for_genes.csv",
                        "neg_all_events_for_genes.csv",
                        "gene_pos_summary.csv",
                        "selection_events.csv")

    missing <- !file.exists(file.path(INPUT_DIR, required_files))
    if (any(missing)) {
        cat(sprintf("  ERROR: Missing files: %s\n", paste(required_files[missing], collapse = ", ")))
        cat("  Skipping this group.\n")
        next
    }

    pos_events <- read.csv(file.path(INPUT_DIR, "pos_all_events_for_genes.csv"),
                           stringsAsFactors = FALSE)
    neg_events <- read.csv(file.path(INPUT_DIR, "neg_all_events_for_genes.csv"),
                           stringsAsFactors = FALSE)
    pos_summary <- read.csv(file.path(INPUT_DIR, "gene_pos_summary.csv"),
                            stringsAsFactors = FALSE)
    selection_events <- read.csv(file.path(INPUT_DIR, "selection_events.csv"),
                                 stringsAsFactors = FALSE)

    cat(sprintf("  Pos events loaded: %d\n", nrow(pos_events)))
    cat(sprintf("  Neg events loaded: %d\n", nrow(neg_events)))
    cat(sprintf("  Genes in summary:  %d\n", nrow(pos_summary)))
    cat(sprintf("  Selection events:  %d\n", nrow(selection_events)))

    # ========================================================================
    # Filter to significant events
    # ========================================================================
    pos_sig <- pos_events %>% filter(significant == TRUE) %>% add_direction()
    neg_sig <- neg_events %>% filter(significant == TRUE) %>% add_direction()

    cat(sprintf("  Pos significant events: %d (in %d genes)\n",
                nrow(pos_sig), n_distinct(pos_sig$geneSymbol)))
    cat(sprintf("  Neg significant events: %d (in %d genes)\n",
                nrow(neg_sig), n_distinct(neg_sig$geneSymbol)))

    # ========================================================================
    # Headline statistics
    # ========================================================================
    total_genes <- nrow(pos_summary)
    genes_with_sig_pos <- pos_summary %>% filter(n_sig_pos > 0) %>% nrow()
    genes_without_sig_pos <- total_genes - genes_with_sig_pos

    cat(sprintf("  Genes WITH significant Pos events:    %d (%.1f%%)\n",
                genes_with_sig_pos, 100 * genes_with_sig_pos / total_genes))
    cat(sprintf("  Genes WITHOUT significant Pos events: %d (%.1f%%)\n",
                genes_without_sig_pos, 100 * genes_without_sig_pos / total_genes))

    # ========================================================================
    # Gene classification (0/1/2/3+ sig Pos events)
    # ========================================================================
    gene_class <- pos_summary %>%
        mutate(sig_pos_bin = case_when(
            n_sig_pos == 0 ~ "0",
            n_sig_pos == 1 ~ "1",
            n_sig_pos == 2 ~ "2",
            TRUE ~ "3+"
        )) %>%
        mutate(sig_pos_bin = factor(sig_pos_bin, levels = c("0", "1", "2", "3+")))

    gene_class_counts <- gene_class %>%
        count(sig_pos_bin, .drop = FALSE) %>%
        rename(n_genes = n)

    # ========================================================================
    # Comparison data frames
    # ========================================================================

    # Event type comparison
    neg_type_counts <- neg_sig %>%
        count(event_type, name = "count") %>%
        mutate(condition = "Neg_vs_Parental")

    pos_type_counts <- pos_sig %>%
        count(event_type, name = "count") %>%
        mutate(condition = "Pos_vs_Parental")

    type_comparison <- bind_rows(neg_type_counts, pos_type_counts) %>%
        complete(event_type = names(event_colors), condition, fill = list(count = 0)) %>%
        mutate(event_type = factor(event_type, levels = names(event_colors)))

    # Size category comparison (SE events only)
    neg_size_counts <- neg_sig %>%
        filter(event_type == "SE") %>%
        count(size_category, name = "count") %>%
        mutate(condition = "Neg_vs_Parental")

    pos_size_counts <- pos_sig %>%
        filter(event_type == "SE") %>%
        count(size_category, name = "count") %>%
        mutate(condition = "Pos_vs_Parental")

    size_comparison <- bind_rows(neg_size_counts, pos_size_counts) %>%
        complete(size_category = names(size_colors), condition, fill = list(count = 0)) %>%
        mutate(size_category = factor(size_category, levels = names(size_colors)))

    # Direction comparison
    neg_dir_counts <- neg_sig %>%
        count(direction, name = "count") %>%
        mutate(condition = "Neg_vs_Parental")

    pos_dir_counts <- pos_sig %>%
        count(direction, name = "count") %>%
        mutate(condition = "Pos_vs_Parental")

    dir_comparison <- bind_rows(neg_dir_counts, pos_dir_counts) %>%
        complete(direction = names(direction_colors), condition, fill = list(count = 0))

    # ========================================================================
    # Per-gene profile table
    # ========================================================================
    neg_profile <- neg_sig %>%
        group_by(geneSymbol) %>%
        summarise(
            neg_n_sig = n(),
            neg_event_types = paste(sort(unique(event_type)), collapse = ","),
            neg_n_SE = sum(event_type == "SE"),
            neg_n_microexon = sum(size_category == "Microexon (0-30bp)", na.rm = TRUE),
            neg_n_inclusion = sum(dPSI > 0),
            neg_n_exclusion = sum(dPSI < 0),
            neg_mean_dPSI = round(mean(dPSI), 4),
            neg_mean_abs_dPSI = round(mean(abs(dPSI)), 4),
            .groups = "drop"
        )

    pos_profile <- pos_sig %>%
        group_by(geneSymbol) %>%
        summarise(
            pos_n_sig = n(),
            pos_event_types = paste(sort(unique(event_type)), collapse = ","),
            pos_n_SE = sum(event_type == "SE"),
            pos_n_microexon = sum(size_category == "Microexon (0-30bp)", na.rm = TRUE),
            pos_n_inclusion = sum(dPSI > 0),
            pos_n_exclusion = sum(dPSI < 0),
            pos_mean_dPSI = round(mean(dPSI), 4),
            pos_mean_abs_dPSI = round(mean(abs(dPSI)), 4),
            .groups = "drop"
        )

    gene_profile <- pos_summary %>%
        select(GeneID, geneSymbol) %>%
        left_join(neg_profile, by = "geneSymbol") %>%
        left_join(pos_profile, by = "geneSymbol") %>%
        mutate(across(starts_with("neg_n"), ~replace_na(.x, 0L)),
               across(starts_with("pos_n"), ~replace_na(.x, 0L)),
               neg_event_types = replace_na(neg_event_types, ""),
               pos_event_types = replace_na(pos_event_types, ""))

    # ========================================================================
    # Generate 6 plots
    # ========================================================================
    cat("  Generating plots...\n")

    # --- Plot 1: Gene classification barplot ---
    p1 <- ggplot(gene_class_counts, aes(x = sig_pos_bin, y = n_genes)) +
        geom_col(fill = "#e67e22", width = 0.6) +
        geom_text(aes(label = n_genes), vjust = -0.5, size = 4.5) +
        labs(title = paste0(group_label, " Genes:\nNumber of Significant Pos Events"),
             x = "Number of significant Pos_vs_Parental events",
             y = "Number of genes") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"))

    ggsave(file.path(OUTPUT_DIR, "gene_classification_bar.pdf"),
           p1, width = 6, height = 5)

    # --- Plot 2: Event type comparison ---
    p2 <- ggplot(type_comparison, aes(x = event_type, y = count, fill = condition)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(aes(label = count), position = position_dodge(width = 0.7),
                  vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = condition_colors, name = "Comparison") +
        labs(title = paste0("Significant Events by Type: Neg vs Pos\n(", group_label, " genes)"),
             x = "Event type", y = "Number of significant events") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top",
              panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"))

    ggsave(file.path(OUTPUT_DIR, "event_type_comparison.pdf"),
           p2, width = 7, height = 5)

    # --- Plot 3: Size category comparison (SE events only) ---
    p3 <- ggplot(size_comparison, aes(x = size_category, y = count, fill = condition)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(aes(label = count), position = position_dodge(width = 0.7),
                  vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = condition_colors, name = "Comparison") +
        labs(title = paste0("Significant SE Events by Size: Neg vs Pos\n(", group_label, " genes)"),
             subtitle = "Skipped exon events only",
             x = "Size category", y = "Number of significant SE events") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top",
              panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(angle = 15, hjust = 1))

    ggsave(file.path(OUTPUT_DIR, "size_category_comparison.pdf"),
           p3, width = 7, height = 5)

    # --- Plot 4: Direction comparison ---
    p4 <- ggplot(dir_comparison, aes(x = direction, y = count, fill = condition)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(aes(label = count), position = position_dodge(width = 0.7),
                  vjust = -0.5, size = 3.5) +
        scale_fill_manual(values = condition_colors, name = "Comparison") +
        labs(title = paste0("Significant Events by Direction: Neg vs Pos\n(", group_label, " genes)"),
             x = "Direction", y = "Number of significant events") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top",
              panel.grid.major.x = element_blank(),
              plot.title = element_text(size = 12, face = "bold"))

    ggsave(file.path(OUTPUT_DIR, "direction_comparison.pdf"),
           p4, width = 7, height = 5)

    # --- Plot 5: Gene-event heatmap for genes WITH sig Pos events ---
    genes_with_pos <- gene_profile %>% filter(pos_n_sig > 0)

    if (nrow(genes_with_pos) > 0) {
        # Limit to top 50 genes by total events if too many for readability
        max_heatmap_genes <- 50
        if (nrow(genes_with_pos) > max_heatmap_genes) {
            genes_with_pos <- genes_with_pos %>%
                mutate(total = neg_n_sig + pos_n_sig) %>%
                arrange(desc(total)) %>%
                head(max_heatmap_genes)
            heatmap_subtitle <- sprintf(
                "Top %d genes by total sig events (of %d with sig Pos)",
                max_heatmap_genes,
                sum(gene_profile$pos_n_sig > 0))
        } else {
            heatmap_subtitle <- sprintf(
                "Genes with at least 1 significant Pos event (n = %d)",
                nrow(genes_with_pos))
        }

        neg_gene_detail <- neg_sig %>%
            filter(geneSymbol %in% genes_with_pos$geneSymbol) %>%
            mutate(condition = "Neg_vs_Parental")

        pos_gene_detail <- pos_sig %>%
            filter(geneSymbol %in% genes_with_pos$geneSymbol) %>%
            mutate(condition = "Pos_vs_Parental")

        heatmap_data <- bind_rows(neg_gene_detail, pos_gene_detail) %>%
            group_by(geneSymbol, condition, event_type) %>%
            summarise(count = n(), .groups = "drop") %>%
            complete(geneSymbol = unique(genes_with_pos$geneSymbol),
                     condition = c("Neg_vs_Parental", "Pos_vs_Parental"),
                     event_type = names(event_colors),
                     fill = list(count = 0))

        gene_order <- genes_with_pos %>%
            mutate(total = neg_n_sig + pos_n_sig) %>%
            arrange(desc(total)) %>%
            pull(geneSymbol)

        heatmap_data <- heatmap_data %>%
            mutate(geneSymbol = factor(geneSymbol, levels = rev(gene_order)),
                   event_type = factor(event_type, levels = names(event_colors)))

        p5 <- ggplot(heatmap_data, aes(x = event_type, y = geneSymbol, fill = count)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = ifelse(count > 0, count, "")), size = 3) +
            facet_wrap(~condition) +
            scale_fill_gradient(low = "white", high = "#2c3e50", name = "# Events") +
            labs(title = paste0("Significant Events per Gene: Neg vs Pos\n(", group_label, ")"),
                 subtitle = heatmap_subtitle,
                 x = "Event type", y = "") +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(size = 12, face = "bold"),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  strip.text = element_text(face = "bold", size = 11),
                  panel.grid = element_blank())

        heatmap_height <- max(6, min(nrow(genes_with_pos) * 0.3 + 3, 30))
        ggsave(file.path(OUTPUT_DIR, "gene_event_heatmap.pdf"),
               p5, width = 10, height = heatmap_height)
    } else {
        cat("  Skipped gene_event_heatmap.pdf (no genes with sig Pos events)\n")
    }

    # --- Plot 6: Combined panel ---
    combined <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                             top = grid::textGrob(
                                 paste0(group_label, " Genes: Pos Event Characterization"),
                                 gp = grid::gpar(fontsize = 14, fontface = "bold")))

    ggsave(file.path(OUTPUT_DIR, "combined_panel.pdf"),
           combined, width = 14, height = 12)

    cat("  Saved 6 PDFs\n")

    # ========================================================================
    # Save CSV tables
    # ========================================================================

    # Table 1: All significant Pos events
    write.csv(pos_sig, file.path(OUTPUT_DIR, "pos_significant_events_detail.csv"),
              row.names = FALSE)

    # Table 2: Per-gene splicing profile comparison
    write.csv(gene_profile, file.path(OUTPUT_DIR, "gene_splicing_profile_comparison.csv"),
              row.names = FALSE)

    # Table 3: Summary statistics
    summary_stats <- data.frame(
        metric = c(
            "group_name", "group_label",
            "total_genes",
            "genes_with_sig_pos_events",
            "genes_without_sig_pos_events",
            "pct_genes_with_sig_pos",
            # Pos event counts
            "pos_total_sig_events",
            "pos_sig_SE", "pos_sig_A3SS", "pos_sig_A5SS", "pos_sig_MXE", "pos_sig_RI",
            "pos_sig_microexon", "pos_sig_small", "pos_sig_regular",
            "pos_sig_inclusion", "pos_sig_exclusion",
            # Neg event counts
            "neg_total_sig_events",
            "neg_sig_SE", "neg_sig_A3SS", "neg_sig_A5SS", "neg_sig_MXE", "neg_sig_RI",
            "neg_sig_microexon", "neg_sig_small", "neg_sig_regular",
            "neg_sig_inclusion", "neg_sig_exclusion"
        ),
        value = c(
            group_name, group_label,
            total_genes,
            genes_with_sig_pos,
            genes_without_sig_pos,
            round(100 * genes_with_sig_pos / total_genes, 1),
            # Pos
            nrow(pos_sig),
            sum(pos_sig$event_type == "SE"),
            sum(pos_sig$event_type == "A3SS"),
            sum(pos_sig$event_type == "A5SS"),
            sum(pos_sig$event_type == "MXE"),
            sum(pos_sig$event_type == "RI"),
            sum(pos_sig$event_type == "SE" & pos_sig$size_category == "Microexon (0-30bp)"),
            sum(pos_sig$event_type == "SE" & pos_sig$size_category == "Small (30-50bp)"),
            sum(pos_sig$event_type == "SE" & pos_sig$size_category == "Regular (>50bp)"),
            sum(pos_sig$dPSI > 0),
            sum(pos_sig$dPSI < 0),
            # Neg
            nrow(neg_sig),
            sum(neg_sig$event_type == "SE"),
            sum(neg_sig$event_type == "A3SS"),
            sum(neg_sig$event_type == "A5SS"),
            sum(neg_sig$event_type == "MXE"),
            sum(neg_sig$event_type == "RI"),
            sum(neg_sig$event_type == "SE" & neg_sig$size_category == "Microexon (0-30bp)"),
            sum(neg_sig$event_type == "SE" & neg_sig$size_category == "Small (30-50bp)"),
            sum(neg_sig$event_type == "SE" & neg_sig$size_category == "Regular (>50bp)"),
            sum(neg_sig$dPSI > 0),
            sum(neg_sig$dPSI < 0)
        ),
        stringsAsFactors = FALSE
    )

    write.csv(summary_stats, file.path(OUTPUT_DIR, "summary_statistics.csv"),
              row.names = FALSE)

    cat(sprintf("  Saved 3 CSVs\n"))

    # ========================================================================
    # Console summary
    # ========================================================================
    cat(sprintf("\n  --- %s SUMMARY ---\n", group_label))
    cat(sprintf("  %d / %d genes (%.0f%%) have significant Pos events\n",
                genes_with_sig_pos, total_genes,
                100 * genes_with_sig_pos / total_genes))

    cat("  Event type breakdown (Neg / Pos sig events):\n")
    for (et in names(event_colors)) {
        n_neg <- sum(neg_sig$event_type == et)
        n_pos <- sum(pos_sig$event_type == et)
        cat(sprintf("    %-5s: %3d / %3d\n", et, n_neg, n_pos))
    }

    cat("  Size breakdown - SE events only (Neg / Pos):\n")
    for (sc in names(size_colors)) {
        n_neg <- sum(neg_sig$event_type == "SE" & neg_sig$size_category == sc)
        n_pos <- sum(pos_sig$event_type == "SE" & pos_sig$size_category == sc)
        cat(sprintf("    %-20s: %3d / %3d\n", sc, n_neg, n_pos))
    }

    cat("  Direction (Neg / Pos):\n")
    cat(sprintf("    Inclusion:  %3d / %3d\n", sum(neg_sig$dPSI > 0), sum(pos_sig$dPSI > 0)))
    cat(sprintf("    Exclusion:  %3d / %3d\n", sum(neg_sig$dPSI < 0), sum(pos_sig$dPSI < 0)))

    # Top genes
    top_pos_genes <- gene_profile %>%
        filter(pos_n_sig > 0) %>%
        arrange(desc(pos_n_sig)) %>%
        head(10)

    if (nrow(top_pos_genes) > 0) {
        cat("  Top genes by # significant Pos events:\n")
        for (i in seq_len(nrow(top_pos_genes))) {
            g <- top_pos_genes[i, ]
            cat(sprintf("    %-10s: %d Pos sig, %d Neg sig (Pos types: %s)\n",
                        g$geneSymbol, g$pos_n_sig, g$neg_n_sig, g$pos_event_types))
        }
    }

    cat(sprintf("  Output: %s\n", OUTPUT_DIR))
}

# =============================================================================
# Final cross-group comparison
# =============================================================================
cat("\n\n========================================================\n")
cat("  CROSS-GROUP COMPARISON\n")
cat("========================================================\n\n")

cat(sprintf("  %-30s  %6s  %8s  %8s  %8s  %s\n",
            "Group", "Genes", "Neg Sig", "Pos Sig", "w/ Pos", "% w/ Pos"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (group_name in names(GROUP_LABELS)) {
    stats_file <- file.path(INPUT_BASE, group_name, "pos_significant_deep", "summary_statistics.csv")
    if (file.exists(stats_file)) {
        stats <- read.csv(stats_file, stringsAsFactors = FALSE)
        get_val <- function(m) stats$value[stats$metric == m]
        cat(sprintf("  %-30s  %6s  %8s  %8s  %8s  %5s%%\n",
                    get_val("group_label"),
                    get_val("total_genes"),
                    get_val("neg_total_sig_events"),
                    get_val("pos_total_sig_events"),
                    get_val("genes_with_sig_pos_events"),
                    get_val("pct_genes_with_sig_pos")))
    }
}

cat("\n14h complete.\n")
