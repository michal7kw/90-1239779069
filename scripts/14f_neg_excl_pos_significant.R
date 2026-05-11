#!/usr/bin/env Rscript
#
# 14f_neg_excl_pos_significant.R
#
# Deep characterization of significant Pos_vs_Parental events in genes that
# have Neg-exclusive microexon SE events (from 14b/14c pipeline).
#
# Question: Of the 81 Neg-exclusive microexon genes, how many have ANY
# significant splicing events in Pos_vs_Parental? What is the composition
# (event types, sizes, direction) vs the Neg profile?
#
# Input: Pre-computed CSVs from 14b/14c (no raw rMATS re-loading)
# Output: 6 PDFs + 3 CSVs in pos_significant_deep/ subdirectory

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# ============================================================================
# Section 1: Configuration
# ============================================================================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
INPUT_DIR <- file.path(BASE_DIR, "results/14_todo_analysis/task2_neg_exclusive_microexon")
OUTPUT_DIR <- file.path(INPUT_DIR, "pos_significant_deep")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Color palettes (matching 14c conventions)
event_colors <- c(SE = "#1f77b4", A5SS = "#ff7f0e", A3SS = "#2ca02c",
                  MXE = "#d62728", RI = "#9467bd")

size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")

direction_colors <- c("Inclusion (dPSI > 0)" = "#2ecc71",
                      "Exclusion (dPSI < 0)" = "#e74c3c")

condition_colors <- c("Neg_vs_Parental" = "#3498db",
                      "Pos_vs_Parental" = "#e67e22")

# ============================================================================
# Section 2: Load CSVs from 14b/14c
# ============================================================================

cat("=== Loading input data ===\n")

pos_events <- read.csv(file.path(INPUT_DIR, "neg_exclusive_genes_pos_events.csv"),
                       stringsAsFactors = FALSE)
neg_events <- read.csv(file.path(INPUT_DIR, "neg_all_events_for_genes.csv"),
                       stringsAsFactors = FALSE)
pos_summary <- read.csv(file.path(INPUT_DIR, "neg_exclusive_genes_pos_summary.csv"),
                        stringsAsFactors = FALSE)
neg_excl_microexon <- read.csv(file.path(INPUT_DIR, "neg_exclusive_microexon_events.csv"),
                               stringsAsFactors = FALSE)

cat(sprintf("  Pos events loaded: %d\n", nrow(pos_events)))
cat(sprintf("  Neg events loaded: %d\n", nrow(neg_events)))
cat(sprintf("  Genes in summary:  %d\n", nrow(pos_summary)))
cat(sprintf("  Neg-exclusive microexon events: %d\n", nrow(neg_excl_microexon)))

# ============================================================================
# Section 3: Filter to significant events, add direction column
# ============================================================================

cat("\n=== Filtering to significant events ===\n")

add_direction <- function(df) {
  df %>%
    mutate(direction = ifelse(dPSI > 0,
                              "Inclusion (dPSI > 0)",
                              "Exclusion (dPSI < 0)"))
}

pos_sig <- pos_events %>%
  filter(significant == TRUE) %>%
  add_direction()

neg_sig <- neg_events %>%
  filter(significant == TRUE) %>%
  add_direction()

cat(sprintf("  Pos significant events: %d (in %d genes)\n",
            nrow(pos_sig), n_distinct(pos_sig$geneSymbol)))
cat(sprintf("  Neg significant events: %d (in %d genes)\n",
            nrow(neg_sig), n_distinct(neg_sig$geneSymbol)))

# ============================================================================
# Section 4: Headline statistics
# ============================================================================

cat("\n=== Headline Statistics ===\n")

total_genes <- nrow(pos_summary)
genes_with_sig_pos <- pos_summary %>% filter(n_sig_pos > 0) %>% nrow()
genes_without_sig_pos <- total_genes - genes_with_sig_pos

cat(sprintf("  Total Neg-exclusive microexon genes: %d\n", total_genes))
cat(sprintf("  Genes WITH significant Pos events:   %d (%.1f%%)\n",
            genes_with_sig_pos, 100 * genes_with_sig_pos / total_genes))
cat(sprintf("  Genes WITHOUT significant Pos events: %d (%.1f%%)\n",
            genes_without_sig_pos, 100 * genes_without_sig_pos / total_genes))

# ============================================================================
# Section 5: Gene classification (0/1/2/3+ sig Pos events)
# ============================================================================

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

cat("\n  Gene classification by # sig Pos events:\n")
for (i in seq_len(nrow(gene_class_counts))) {
  cat(sprintf("    %s sig events: %d genes\n",
              gene_class_counts$sig_pos_bin[i],
              gene_class_counts$n_genes[i]))
}

# ============================================================================
# Section 6: Comparison data frames (long-format for grouped bar charts)
# ============================================================================

# --- Event type comparison ---
neg_type_counts <- neg_sig %>%
  count(event_type, name = "count") %>%
  mutate(condition = "Neg_vs_Parental")

pos_type_counts <- pos_sig %>%
  count(event_type, name = "count") %>%
  mutate(condition = "Pos_vs_Parental")

type_comparison <- bind_rows(neg_type_counts, pos_type_counts) %>%
  complete(event_type = names(event_colors), condition, fill = list(count = 0)) %>%
  mutate(event_type = factor(event_type, levels = names(event_colors)))

# --- Size category comparison (SE events only) ---
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

# --- Direction comparison ---
neg_dir_counts <- neg_sig %>%
  count(direction, name = "count") %>%
  mutate(condition = "Neg_vs_Parental")

pos_dir_counts <- pos_sig %>%
  count(direction, name = "count") %>%
  mutate(condition = "Pos_vs_Parental")

dir_comparison <- bind_rows(neg_dir_counts, pos_dir_counts) %>%
  complete(direction = names(direction_colors), condition, fill = list(count = 0))

# ============================================================================
# Section 7: Per-gene profile table
# ============================================================================

# Build Neg profile per gene
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

# Build Pos profile per gene
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

# Merge on all 81 genes
gene_profile <- pos_summary %>%
  select(GeneID, geneSymbol) %>%
  left_join(neg_profile, by = "geneSymbol") %>%
  left_join(pos_profile, by = "geneSymbol") %>%
  # Fill NAs for genes without significant events in either condition
  mutate(across(starts_with("neg_n"), ~replace_na(.x, 0L)),
         across(starts_with("pos_n"), ~replace_na(.x, 0L)),
         neg_event_types = replace_na(neg_event_types, ""),
         pos_event_types = replace_na(pos_event_types, ""))

cat(sprintf("\n  Gene profile table: %d genes\n", nrow(gene_profile)))

# ============================================================================
# Section 8: Generate all 6 plots
# ============================================================================

cat("\n=== Generating plots ===\n")

# --- Plot 1: Gene classification barplot ---
p1 <- ggplot(gene_class_counts, aes(x = sig_pos_bin, y = n_genes)) +
  geom_col(fill = "#e67e22", width = 0.6) +
  geom_text(aes(label = n_genes), vjust = -0.5, size = 4.5) +
  labs(title = "Neg-Exclusive Microexon Genes:\nNumber of Significant Pos Events",
       x = "Number of significant Pos_vs_Parental events",
       y = "Number of genes") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "gene_classification_bar.pdf"),
       p1, width = 6, height = 5)
cat("  Saved gene_classification_bar.pdf\n")

# --- Plot 2: Event type comparison ---
p2 <- ggplot(type_comparison, aes(x = event_type, y = count, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = condition_colors, name = "Comparison") +
  labs(title = "Significant Events by Type: Neg vs Pos",
       x = "Event type", y = "Number of significant events") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "event_type_comparison.pdf"),
       p2, width = 7, height = 5)
cat("  Saved event_type_comparison.pdf\n")

# --- Plot 3: Size category comparison (SE events only) ---
p3 <- ggplot(size_comparison, aes(x = size_category, y = count, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = condition_colors, name = "Comparison") +
  labs(title = "Significant SE Events by Size: Neg vs Pos",
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
cat("  Saved size_category_comparison.pdf\n")

# --- Plot 4: Direction comparison ---
p4 <- ggplot(dir_comparison, aes(x = direction, y = count, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = condition_colors, name = "Comparison") +
  labs(title = "Significant Events by Direction: Neg vs Pos",
       x = "Direction", y = "Number of significant events") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "direction_comparison.pdf"),
       p4, width = 7, height = 5)
cat("  Saved direction_comparison.pdf\n")

# --- Plot 5: Gene-event heatmap for genes WITH sig Pos events ---
genes_with_pos <- gene_profile %>% filter(pos_n_sig > 0)

if (nrow(genes_with_pos) > 0) {
  # Build long-format per-gene event detail for heatmap
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

  # Order genes by total sig events (Neg + Pos)
  gene_order <- gene_profile %>%
    filter(pos_n_sig > 0) %>%
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
    labs(title = "Significant Events per Gene: Neg vs Pos",
         subtitle = sprintf("Genes with at least 1 significant Pos event (n = %d)",
                            nrow(genes_with_pos)),
         x = "Event type", y = "") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 11),
          panel.grid = element_blank())

  ggsave(file.path(OUTPUT_DIR, "gene_event_heatmap.pdf"),
         p5, width = 10, height = max(6, nrow(genes_with_pos) * 0.3 + 3))
  cat("  Saved gene_event_heatmap.pdf\n")
} else {
  cat("  Skipped gene_event_heatmap.pdf (no genes with sig Pos events)\n")
}

# --- Plot 6: Combined panel ---
combined <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                         top = grid::textGrob(
                           "Neg-Exclusive Microexon Genes: Pos Event Characterization",
                           gp = grid::gpar(fontsize = 14, fontface = "bold")))

ggsave(file.path(OUTPUT_DIR, "combined_panel.pdf"),
       combined, width = 14, height = 12)
cat("  Saved combined_panel.pdf\n")

# ============================================================================
# Section 9: Save CSV tables
# ============================================================================

cat("\n=== Saving CSV tables ===\n")

# Table 1: All significant Pos events with full detail
write.csv(pos_sig, file.path(OUTPUT_DIR, "pos_significant_events_detail.csv"),
          row.names = FALSE)
cat(sprintf("  pos_significant_events_detail.csv: %d events\n", nrow(pos_sig)))

# Table 2: Per-gene splicing profile comparison
write.csv(gene_profile, file.path(OUTPUT_DIR, "gene_splicing_profile_comparison.csv"),
          row.names = FALSE)
cat(sprintf("  gene_splicing_profile_comparison.csv: %d genes\n", nrow(gene_profile)))

# Table 3: Summary statistics
summary_stats <- data.frame(
  metric = c(
    "total_neg_exclusive_genes",
    "genes_with_sig_pos_events",
    "genes_without_sig_pos_events",
    "pct_genes_with_sig_pos",
    # Pos event counts
    "pos_total_sig_events",
    "pos_sig_SE", "pos_sig_A3SS", "pos_sig_A5SS", "pos_sig_MXE", "pos_sig_RI",
    "pos_sig_microexon", "pos_sig_small", "pos_sig_regular",
    "pos_sig_inclusion", "pos_sig_exclusion",
    # Neg event counts (for comparison)
    "neg_total_sig_events",
    "neg_sig_SE", "neg_sig_A3SS", "neg_sig_A5SS", "neg_sig_MXE", "neg_sig_RI",
    "neg_sig_microexon", "neg_sig_small", "neg_sig_regular",
    "neg_sig_inclusion", "neg_sig_exclusion"
  ),
  value = c(
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
  )
)

write.csv(summary_stats, file.path(OUTPUT_DIR, "summary_statistics.csv"),
          row.names = FALSE)
cat(sprintf("  summary_statistics.csv: %d metrics\n", nrow(summary_stats)))

# ============================================================================
# Section 10: Console summary
# ============================================================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("  SUMMARY: Pos Significant Events in Neg-Exclusive Genes\n")
cat("=", rep("=", 60), "\n", sep = "")
cat(sprintf("\n  %d / %d genes (%.0f%%) have significant Pos events\n",
            genes_with_sig_pos, total_genes,
            100 * genes_with_sig_pos / total_genes))
cat(sprintf("  %d genes have NO significant Pos events\n\n", genes_without_sig_pos))

cat("  Event type breakdown (Neg / Pos sig events):\n")
for (et in names(event_colors)) {
  n_neg <- sum(neg_sig$event_type == et)
  n_pos <- sum(pos_sig$event_type == et)
  cat(sprintf("    %-5s: %3d / %3d\n", et, n_neg, n_pos))
}

cat("\n  Size breakdown - SE events only (Neg / Pos):\n")
for (sc in names(size_colors)) {
  n_neg <- sum(neg_sig$event_type == "SE" & neg_sig$size_category == sc)
  n_pos <- sum(pos_sig$event_type == "SE" & pos_sig$size_category == sc)
  cat(sprintf("    %-20s: %3d / %3d\n", sc, n_neg, n_pos))
}

cat("\n  Direction (Neg / Pos):\n")
n_neg_inc <- sum(neg_sig$dPSI > 0)
n_neg_exc <- sum(neg_sig$dPSI < 0)
n_pos_inc <- sum(pos_sig$dPSI > 0)
n_pos_exc <- sum(pos_sig$dPSI < 0)
cat(sprintf("    Inclusion:  %3d / %3d\n", n_neg_inc, n_pos_inc))
cat(sprintf("    Exclusion:  %3d / %3d\n", n_neg_exc, n_pos_exc))

# Top genes by Pos sig event count
top_pos_genes <- gene_profile %>%
  filter(pos_n_sig > 0) %>%
  arrange(desc(pos_n_sig)) %>%
  head(10)

cat("\n  Top genes by # significant Pos events:\n")
for (i in seq_len(nrow(top_pos_genes))) {
  g <- top_pos_genes[i, ]
  cat(sprintf("    %-10s: %d Pos sig, %d Neg sig (Pos types: %s)\n",
              g$geneSymbol, g$pos_n_sig, g$neg_n_sig, g$pos_event_types))
}

cat("\n  Output saved to:", OUTPUT_DIR, "\n")
cat("  Files: 6 PDFs + 3 CSVs\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("Done!\n")
