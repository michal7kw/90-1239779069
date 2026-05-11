#!/usr/bin/env Rscript
# =============================================================================
# Script: 15b_sensitivity_highconf.R
# Purpose: High-confidence (gene+length) subset analysis and enhanced plots
#          for paper microexon sensitivity cross-reference
#
# Reads the already-matched data from step 15 and produces:
#   - Filtered CSVs for gene+length matches only
#   - Match-type breakdown plots (gene+length vs gene_only)
#   - High-confidence-only detection/distribution plots
#   - Side-by-side comparison: all matches vs gene+length only
#
# Input:  results/15_paper_microexon_sensitivity/ (existing CSVs)
# Output: results/15_paper_microexon_sensitivity/high_confidence/ (new subdir)
#       + 3 new PDFs in the main output directory
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# --- Configuration -----------------------------------------------------------
base_dir    <- getwd()
main_dir    <- file.path(base_dir, "results", "15_paper_microexon_sensitivity")
hc_dir      <- file.path(main_dir, "high_confidence")
comparisons <- c("Neg_vs_Parental", "KO_vs_Parental", "Pos_vs_Parental")

# Tier ordering and colors (same as original script)
tier_order  <- c("HS", "LS", "CS", "CR", "NonResponding", "Others")
tier_colors <- c(
  "HS"            = "#D32F2F",
  "LS"            = "#F57C00",
  "CS"            = "#388E3C",
  "CR"            = "#1976D2",
  "NonResponding" = "#757575",
  "Others"        = "#9E9E9E"
)
match_colors <- c("gene+length" = "steelblue", "gene_only" = "grey60")

dir.create(hc_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== High-Confidence Subset Analysis (gene+length matches) ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Read existing matched data (do NOT recompute)
# =============================================================================
cat("--- Reading existing step 15 output ---\n")

matched_all <- read.csv(file.path(main_dir, "matched_events_by_comparison.csv"),
                        stringsAsFactors = FALSE)

cat("Total matched pairs:", nrow(matched_all), "\n")
cat("  gene+length:", sum(matched_all$match_type == "gene+length"), "\n")
cat("  gene_only:  ", sum(matched_all$match_type == "gene_only"), "\n\n")

# Read paper reference for unmatched computation
paper_table1 <- file.path(base_dir, "data", "Supplementary Table 1.csv")
paper_table8 <- file.path(base_dir, "data", "Supplementary Table 8.csv")

t1 <- read.csv(paper_table1, skip = 2, stringsAsFactors = FALSE) %>%
  select(EVENT, TIER) %>%
  filter(EVENT != "", !is.na(EVENT))

t8 <- read.csv(paper_table8, skip = 3, stringsAsFactors = FALSE) %>%
  select(GENE, EVENT, COORD, LENGTH) %>%
  filter(EVENT != "", !is.na(EVENT))
t8$LENGTH <- as.numeric(t8$LENGTH)

paper <- t8 %>%
  inner_join(t1, by = "EVENT") %>%
  filter(!is.na(LENGTH), !is.na(GENE), GENE != "")

cat("Paper reference:", nrow(paper), "microexons\n\n")

# =============================================================================
# 2. Filter to gene+length matches only
# =============================================================================
cat("--- Filtering to gene+length (high-confidence) matches ---\n")

hc_df <- matched_all %>% filter(match_type == "gene+length")

cat("High-confidence matched pairs:", nrow(hc_df), "\n")
cat("Unique paper events (gene+length):", n_distinct(hc_df$paper_EVENT), "\n\n")

# =============================================================================
# 3. High-confidence sensitivity summary
# =============================================================================
cat("--- Computing high-confidence summaries ---\n")

tier_totals <- paper %>% count(TIER, name = "tier_total")

hc_summary <- hc_df %>%
  group_by(paper_TIER, Comparison) %>%
  summarise(
    n_detected      = n(),
    n_gene_length   = n(),   # all rows are gene+length by construction
    n_gene_only     = 0L,
    mean_dPSI       = mean(DeltaPSI, na.rm = TRUE),
    median_dPSI     = median(DeltaPSI, na.rm = TRUE),
    mean_abs_dPSI   = mean(abs(DeltaPSI), na.rm = TRUE),
    pct_negative_dPSI = mean(DeltaPSI < 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  left_join(tier_totals, by = c("paper_TIER" = "TIER")) %>%
  mutate(detection_rate_pct = round(n_detected / tier_total * 100, 1)) %>%
  arrange(paper_TIER, Comparison)

# dPSI summary for high-conf
hc_tier_dpsi <- hc_df %>%
  group_by(paper_TIER, Comparison) %>%
  summarise(
    n          = n(),
    mean_dPSI  = round(mean(DeltaPSI, na.rm = TRUE), 4),
    median_dPSI = round(median(DeltaPSI, na.rm = TRUE), 4),
    sd_dPSI    = ifelse(n() > 1, round(sd(DeltaPSI, na.rm = TRUE), 4), NA_real_),
    min_dPSI   = round(min(DeltaPSI, na.rm = TRUE), 4),
    max_dPSI   = round(max(DeltaPSI, na.rm = TRUE), 4),
    .groups = "drop"
  )

# --- Flag multi-locus HC events (same paper event → different mouse exons) ---
hc_loci <- hc_df %>%
  mutate(locus = paste(Chromosome, exonStart_0base, exonEnd, sep = ":")) %>%
  group_by(paper_EVENT) %>%
  summarise(n_loci = n_distinct(locus), .groups = "drop") %>%
  filter(n_loci > 1)

if (nrow(hc_loci) > 0) {
  cat("NOTE:", nrow(hc_loci), "paper events match different mouse exons across comparisons\n")
  cat("  (gene has multiple microexons of similar size; per-comparison stats unaffected)\n")
}

# --- Flag mixed match-type events (gene+length in one comp, gene_only in another) ---
mixed_type <- matched_all %>%
  group_by(paper_EVENT) %>%
  summarise(types = paste(sort(unique(match_type)), collapse = ","), .groups = "drop") %>%
  filter(types == "gene_only,gene+length")

if (nrow(mixed_type) > 0) {
  cat("NOTE:", nrow(mixed_type), "paper events have gene+length match in some comparisons",
      "but gene_only in others\n")
  cat("  (microexon not significant in all comparisons; gene_only fallback catches other exons)\n")
}
cat("\n")

# =============================================================================
# 4. Identify unmatched paper events (no gene+length match in ANY comparison)
# =============================================================================
hc_matched_events <- hc_df %>% distinct(paper_EVENT) %>% pull(paper_EVENT)

hc_unmatched <- paper %>%
  filter(!EVENT %in% hc_matched_events) %>%
  select(GENE, EVENT, COORD, LENGTH, TIER)

cat("Paper events with NO gene+length match:", nrow(hc_unmatched),
    "of", nrow(paper), "\n")
cat("Unmatched by tier:\n")
print(table(hc_unmatched$TIER))
cat("\n")

# =============================================================================
# 5. Write CSVs to high_confidence/ subdir
# =============================================================================
cat("--- Writing high-confidence CSVs ---\n")

# Annotate hc_df with multi-locus flag before writing
hc_df_out <- hc_df %>%
  left_join(hc_loci %>% select(paper_EVENT, n_loci), by = "paper_EVENT") %>%
  mutate(multi_locus = !is.na(n_loci)) %>%
  select(-n_loci)

write.csv(hc_df_out, file.path(hc_dir, "highconf_matched_events.csv"),
          row.names = FALSE)
write.csv(hc_summary, file.path(hc_dir, "highconf_sensitivity_summary.csv"),
          row.names = FALSE)
write.csv(hc_tier_dpsi, file.path(hc_dir, "highconf_tier_dpsi_summary.csv"),
          row.names = FALSE)
write.csv(hc_unmatched, file.path(hc_dir, "highconf_unmatched_paper_events.csv"),
          row.names = FALSE)

cat("  Written: high_confidence/highconf_matched_events.csv\n")
cat("  Written: high_confidence/highconf_sensitivity_summary.csv\n")
cat("  Written: high_confidence/highconf_tier_dpsi_summary.csv\n")
cat("  Written: high_confidence/highconf_unmatched_paper_events.csv\n\n")

# =============================================================================
# 6. Plots — A. Match quality breakdown (main output dir, new filenames)
# =============================================================================
cat("--- Generating plots ---\n")

# 6a. Match type breakdown: grouped bar, per tier × comparison
breakdown_data <- matched_all %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order),
         match_type = factor(match_type, levels = c("gene+length", "gene_only"))) %>%
  count(paper_TIER, Comparison, match_type)

p_breakdown <- ggplot(breakdown_data,
                      aes(x = paper_TIER, y = n, fill = match_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  facet_wrap(~Comparison, nrow = 1) +
  scale_fill_manual(values = match_colors, name = "Match Type") +
  labs(
    title = "Match Type Breakdown by Sensitivity Tier",
    subtitle = "gene+length = high-confidence (same gene AND exon size +/-5bp); gene_only = fallback",
    x = "Sensitivity Tier",
    y = "Number of Matched Events"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(main_dir, "match_type_breakdown.pdf"),
       p_breakdown, width = 14, height = 6)
cat("  Saved: match_type_breakdown.pdf\n")

# 6b. Match type pie/donut per tier
# Compute per-tier composition: gene+length, gene_only, unmatched
# For "unmatched" count per tier: paper events not detected in ANY comparison
all_matched_events <- matched_all %>% distinct(paper_EVENT) %>% pull(paper_EVENT)

tier_composition <- paper %>%
  mutate(status = case_when(
    EVENT %in% (matched_all %>% filter(match_type == "gene+length") %>%
                  distinct(paper_EVENT) %>% pull(paper_EVENT)) ~ "gene+length",
    EVENT %in% all_matched_events ~ "gene_only",
    TRUE ~ "unmatched"
  )) %>%
  count(TIER, status) %>%
  mutate(TIER = factor(TIER, levels = tier_order),
         status = factor(status, levels = c("gene+length", "gene_only", "unmatched")))

pie_colors <- c("gene+length" = "steelblue", "gene_only" = "grey60",
                "unmatched" = "#FFCDD2")

p_pie <- ggplot(tier_composition, aes(x = 2, y = n, fill = status)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_wrap(~TIER, nrow = 1) +
  scale_fill_manual(values = pie_colors, name = "Match Status") +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(
    title = "Match Composition by Sensitivity Tier",
    subtitle = "Per-tier breakdown: gene+length (high-conf), gene_only (fallback), unmatched"
  ) +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(main_dir, "match_type_pie.pdf"),
       p_pie, width = 14, height = 5)
cat("  Saved: match_type_pie.pdf\n")

# =============================================================================
# 7. Plots — B. High-confidence-only plots (in high_confidence/ subdir)
# =============================================================================

# 7a. High-conf detection rate bar chart
hc_plot_data <- hc_summary %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order))

p_hc_detection <- ggplot(hc_plot_data,
                         aes(x = paper_TIER, y = detection_rate_pct,
                             fill = paper_TIER)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(n_detected, "/", tier_total)),
            vjust = -0.5, size = 3) +
  facet_wrap(~Comparison, nrow = 1) +
  scale_fill_manual(values = tier_colors, guide = "none") +
  labs(
    title = "Microexon Detection Rate - Gene+Length Matches Only",
    subtitle = "High-confidence matches: same gene AND exon length within +/-5bp",
    x = "Sensitivity Tier",
    y = "Detection Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(hc_dir, "highconf_detection_rate.pdf"),
       p_hc_detection, width = 12, height = 6)
cat("  Saved: high_confidence/highconf_detection_rate.pdf\n")

# 7b. High-conf detection heatmap
hc_heatmap_data <- hc_summary %>%
  select(paper_TIER, Comparison, detection_rate_pct) %>%
  mutate(paper_TIER = factor(paper_TIER, levels = rev(tier_order)),
         Comparison = factor(Comparison, levels = comparisons))

p_hc_heat <- ggplot(hc_heatmap_data,
                    aes(x = Comparison, y = paper_TIER,
                        fill = detection_rate_pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(detection_rate_pct, 0), "%")),
            size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#FFCC80", high = "#D32F2F",
                       midpoint = 15, name = "Detection\nRate (%)") +
  labs(
    title = "Detection Rate Heatmap - Gene+Length Matches Only",
    subtitle = "TIER x Comparison (high-confidence subset)",
    x = "Comparison",
    y = "Sensitivity Tier"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(hc_dir, "highconf_detection_heatmap.pdf"),
       p_hc_heat, width = 8, height = 5)
cat("  Saved: high_confidence/highconf_detection_heatmap.pdf\n")

# 7c. High-conf dPSI distribution
hc_dpsi_data <- hc_df %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order))

if (nrow(hc_dpsi_data) > 0) {
  p_hc_dpsi <- ggplot(hc_dpsi_data,
                      aes(x = paper_TIER, y = DeltaPSI, fill = paper_TIER)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~Comparison, nrow = 1) +
    scale_fill_manual(values = tier_colors, guide = "none") +
    labs(
      title = "Delta PSI Distribution - Gene+Length Matches Only",
      subtitle = "High-confidence subset: each point is a matched significant event",
      x = "Sensitivity Tier",
      y = "Delta PSI (rMATS)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(hc_dir, "highconf_dpsi_distribution.pdf"),
         p_hc_dpsi, width = 12, height = 6)
  cat("  Saved: high_confidence/highconf_dpsi_distribution.pdf\n")
}

# 7d. Length comparison scatter: paper LENGTH vs our exon_length
hc_length_data <- hc_df %>%
  filter(!is.na(exon_length), !is.na(paper_LENGTH)) %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order))

if (nrow(hc_length_data) > 0) {
  length_range <- range(c(hc_length_data$paper_LENGTH,
                          hc_length_data$exon_length), na.rm = TRUE)
  axis_min <- max(0, length_range[1] - 5)
  axis_max <- length_range[2] + 5

  p_length <- ggplot(hc_length_data,
                     aes(x = paper_LENGTH, y = exon_length, color = paper_TIER)) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey40") +
    geom_abline(slope = 1, intercept = 5, linetype = "dashed",
                color = "grey70", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = -5, linetype = "dashed",
                color = "grey70", linewidth = 0.4) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = tier_colors, name = "Tier") +
    coord_fixed(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
    annotate("text", x = axis_max - 2, y = axis_max - 8,
             label = "+/-5bp tolerance", color = "grey50", size = 3,
             fontface = "italic") +
    labs(
      title = "Exon Length Concordance: Paper vs Our Data",
      subtitle = "Gene+length matches only - points should cluster on the diagonal",
      x = "Paper Exon Length (bp, human)",
      y = "Our Exon Length (bp, mouse)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(hc_dir, "highconf_length_comparison.pdf"),
         p_length, width = 8, height = 7)
  cat("  Saved: high_confidence/highconf_length_comparison.pdf\n")
}

# =============================================================================
# 8. Plot — C. Comparison: all matches vs gene+length only (main output dir)
# =============================================================================

# Overall unique-event detection rates per tier (across all comparisons)
overall_all <- matched_all %>%
  distinct(paper_EVENT, .keep_all = TRUE) %>%
  count(paper_TIER, name = "n_detected") %>%
  left_join(tier_totals, by = c("paper_TIER" = "TIER")) %>%
  mutate(detection_pct = round(n_detected / tier_total * 100, 1),
         strategy = "All matches")

overall_hc <- hc_df %>%
  distinct(paper_EVENT, .keep_all = TRUE) %>%
  count(paper_TIER, name = "n_detected") %>%
  left_join(tier_totals, by = c("paper_TIER" = "TIER")) %>%
  mutate(detection_pct = round(n_detected / tier_total * 100, 1),
         strategy = "Gene+Length only")

# Ensure all tiers present for both strategies (fill missing with 0)
all_tiers_in_data <- unique(c(overall_all$paper_TIER, overall_hc$paper_TIER))
for (tier in all_tiers_in_data) {
  if (!tier %in% overall_hc$paper_TIER) {
    tt <- tier_totals %>% filter(TIER == tier) %>% pull(tier_total)
    overall_hc <- bind_rows(overall_hc,
      tibble(paper_TIER = tier, n_detected = 0, tier_total = tt,
             detection_pct = 0, strategy = "Gene+Length only"))
  }
}

comparison_data <- bind_rows(overall_all, overall_hc) %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order),
         strategy = factor(strategy, levels = c("All matches", "Gene+Length only")))

strategy_colors <- c("All matches" = "#78909C", "Gene+Length only" = "steelblue")

p_comparison <- ggplot(comparison_data,
                       aes(x = paper_TIER, y = detection_pct, fill = strategy)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(n_detected, "/", tier_total)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3) +
  scale_fill_manual(values = strategy_colors, name = "Match Strategy") +
  labs(
    title = "Detection Rate: All Matches vs Gene+Length Only",
    subtitle = "Overall unique paper events detected across all 3 comparisons",
    x = "Sensitivity Tier",
    y = "Detection Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(main_dir, "detection_comparison_all_vs_highconf.pdf"),
       p_comparison, width = 9, height = 6)
cat("  Saved: detection_comparison_all_vs_highconf.pdf\n")

# =============================================================================
# 9. Print comparison summary
# =============================================================================
cat("\n=== COMPARISON SUMMARY: All Matches vs Gene+Length Only ===\n\n")

summary_table <- comparison_data %>%
  select(paper_TIER, strategy, n_detected, detection_pct) %>%
  pivot_wider(names_from = strategy,
              values_from = c(n_detected, detection_pct),
              values_fill = 0) %>%
  left_join(comparison_data %>% distinct(paper_TIER, tier_total),
            by = "paper_TIER") %>%
  mutate(drop_pct = round(`detection_pct_All matches` -
                            `detection_pct_Gene+Length only`, 1))

cat(sprintf("%-15s %6s %6s %6s %6s %8s\n",
            "Tier", "All_n", "HC_n", "All_%", "HC_%", "Drop_%"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (i in seq_len(nrow(summary_table))) {
  row <- summary_table[i, ]
  cat(sprintf("%-15s %6d %6d %6.1f %6.1f %8.1f\n",
              as.character(row$paper_TIER),
              row$`n_detected_All matches`,
              row$`n_detected_Gene+Length only`,
              row$`detection_pct_All matches`,
              row$`detection_pct_Gene+Length only`,
              row$drop_pct))
}

cat("\nTotal unique paper events matched:\n")
cat("  All matches:     ", n_distinct(matched_all$paper_EVENT), "of", nrow(paper), "\n")
cat("  Gene+Length only: ", n_distinct(hc_df$paper_EVENT), "of", nrow(paper), "\n")
cat("  HC unmatched:     ", nrow(hc_unmatched), "\n")

cat("\n=== Output Files ===\n")
cat("Main dir (new):\n")
cat("  match_type_breakdown.pdf\n")
cat("  match_type_pie.pdf\n")
cat("  detection_comparison_all_vs_highconf.pdf\n")
cat("\nHigh-confidence subdir:\n")
cat("  highconf_matched_events.csv\n")
cat("  highconf_sensitivity_summary.csv\n")
cat("  highconf_tier_dpsi_summary.csv\n")
cat("  highconf_unmatched_paper_events.csv\n")
cat("  highconf_detection_rate.pdf\n")
cat("  highconf_detection_heatmap.pdf\n")
cat("  highconf_dpsi_distribution.pdf\n")
cat("  highconf_length_comparison.pdf\n")

cat("\nDone:", format(Sys.time()), "\n")
