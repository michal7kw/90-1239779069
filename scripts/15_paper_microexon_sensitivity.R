#!/usr/bin/env Rscript
# =============================================================================
# Script: 15_paper_microexon_sensitivity.R
# Purpose: Cross-reference paper microexon sensitivity classes (HS, LS, CS, CR,
#          NonResponding) with our SRRM3 significant splicing events
#
# Paper: Nature Structural & Molecular Biology 2025
#   - 234 microexons classified by SRRM4/3 sensitivity tiers
#   - Supplementary Table 1: EVENT + TIER (authoritative tier assignments)
#   - Supplementary Table 8: GENE + EVENT + COORD + LENGTH + TIER + PSI data
#
# Matching strategy: Gene name (case-insensitive) + exon length (±5bp tolerance)
#   Human gene names are uppercase (RAPGEF6), mouse are title-case (Rapgef6)
#   Microexon lengths are highly conserved across species
#
# Output: results/15_paper_microexon_sensitivity/
# =============================================================================

library(dplyr)
library(purrr)
library(ggplot2)

# --- Configuration -----------------------------------------------------------
base_dir <- getwd()
paper_table1 <- file.path(base_dir, "data", "Supplementary Table 1.csv")
paper_table8 <- file.path(base_dir, "data", "Supplementary Table 8.csv")
sig_dir     <- file.path(base_dir, "results", "10_significant_events")
out_dir     <- file.path(base_dir, "results", "15_paper_microexon_sensitivity")

comparisons <- c("Neg_vs_Parental", "KO_vs_Parental", "Pos_vs_Parental")
length_tolerance <- 5  # bp

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Paper Microexon Sensitivity Cross-Reference ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Read and merge paper data
# =============================================================================
cat("--- Reading paper supplementary tables ---\n")

# Table 1: authoritative EVENT -> TIER mapping (skip 2 header rows)
t1 <- read.csv(paper_table1, skip = 2, stringsAsFactors = FALSE) %>%
  select(EVENT, TIER) %>%
  filter(EVENT != "", !is.na(EVENT))

cat("Table 1:", nrow(t1), "events with tier assignments\n")
cat("  Tier distribution:\n")
print(table(t1$TIER))

# Table 8: EVENT -> GENE, LENGTH, COORD, etc. (skip 3 header rows)
t8 <- read.csv(paper_table8, skip = 3, stringsAsFactors = FALSE) %>%
  select(GENE, EVENT, COORD, LENGTH) %>%
  filter(EVENT != "", !is.na(EVENT))

# LENGTH is already numeric in this file
t8$LENGTH <- as.numeric(t8$LENGTH)

cat("Table 8:", nrow(t8), "events with gene/length info\n")

# Merge: use Table 1 TIER as authoritative
paper <- t8 %>%
  inner_join(t1, by = "EVENT") %>%
  filter(!is.na(LENGTH), !is.na(GENE), GENE != "")

cat("Merged:", nrow(paper), "events with GENE + LENGTH + TIER\n\n")

# =============================================================================
# 2. Read our significant splicing events
# =============================================================================
cat("--- Reading significant splicing events ---\n")

compute_exon_length <- function(df) {
  df %>%
    mutate(exon_length = case_when(
      EventType == "SE"   ~ exonEnd - exonStart_0base,
      EventType == "A3SS" ~ longExonEnd - longExonStart_0base,
      EventType == "A5SS" ~ longExonEnd - longExonStart_0base,
      EventType == "MXE"  ~ X1stExonEnd - X1stExonStart_0base,
      EventType == "RI"   ~ riExonEnd - riExonStart_0base,
      TRUE ~ NA_real_
    ))
}

all_events <- map_dfr(comparisons, function(comp) {
  f <- file.path(sig_dir, paste0(comp, "_significant_events.csv"))
  read.csv(f, stringsAsFactors = FALSE) %>%
    compute_exon_length()
})

cat("Total significant events across 3 comparisons:", nrow(all_events), "\n")
for (comp in comparisons) {
  n <- sum(all_events$Comparison == comp)
  cat("  ", comp, ":", n, "events\n")
}
cat("\n")

# =============================================================================
# 3. Match paper events to our events
# =============================================================================
cat("--- Matching paper events to our splicing data ---\n")

# Prepare paper gene names for case-insensitive matching
paper <- paper %>%
  mutate(gene_upper = toupper(GENE))

all_events <- all_events %>%
  mutate(gene_upper = toupper(GeneSymbol))

# For each paper event, find matches in each comparison
matches <- list()

for (i in seq_len(nrow(paper))) {
  p_gene   <- paper$gene_upper[i]
  p_length <- paper$LENGTH[i]
  p_event  <- paper$EVENT[i]
  p_tier   <- paper$TIER[i]
  p_gene_orig <- paper$GENE[i]

  for (comp in comparisons) {
    comp_events <- all_events %>% filter(Comparison == comp)

    # Primary: gene + length match
    gene_length_match <- comp_events %>%
      filter(gene_upper == p_gene,
             !is.na(exon_length),
             abs(exon_length - p_length) <= length_tolerance)

    # Fallback: gene-only match
    gene_only_match <- comp_events %>%
      filter(gene_upper == p_gene)

    if (nrow(gene_length_match) > 0) {
      # Keep best match by FDR, then |dPSI|
      best <- gene_length_match %>%
        arrange(FDR, desc(abs(DeltaPSI))) %>%
        slice(1)
      matches[[length(matches) + 1]] <- best %>%
        mutate(paper_EVENT = p_event,
               paper_GENE = p_gene_orig,
               paper_LENGTH = p_length,
               paper_TIER = p_tier,
               match_type = "gene+length",
               n_candidates = nrow(gene_length_match))
    } else if (nrow(gene_only_match) > 0) {
      best <- gene_only_match %>%
        arrange(FDR, desc(abs(DeltaPSI))) %>%
        slice(1)
      matches[[length(matches) + 1]] <- best %>%
        mutate(paper_EVENT = p_event,
               paper_GENE = p_gene_orig,
               paper_LENGTH = p_length,
               paper_TIER = p_tier,
               match_type = "gene_only",
               n_candidates = nrow(gene_only_match))
    }
    # else: no match for this comparison — tracked in unmatched summary
  }
}

matched_df <- bind_rows(matches)

cat("Total matched event-comparison pairs:", nrow(matched_df), "\n")
cat("  gene+length matches:", sum(matched_df$match_type == "gene+length"), "\n")
cat("  gene_only matches:", sum(matched_df$match_type == "gene_only"), "\n\n")

# =============================================================================
# 4. Identify unmatched paper events
# =============================================================================
matched_events <- matched_df %>%
  distinct(paper_EVENT) %>%
  pull(paper_EVENT)

unmatched <- paper %>%
  filter(!EVENT %in% matched_events) %>%
  select(GENE, EVENT, COORD, LENGTH, TIER)

cat("Paper events with NO match in any comparison:", nrow(unmatched), "of", nrow(paper), "\n")
cat("Unmatched by tier:\n")
print(table(unmatched$TIER))
cat("\n")

# =============================================================================
# 5. Summary statistics
# =============================================================================
cat("--- Generating summaries ---\n")

# 5a. Sensitivity summary: per TIER × comparison
tier_totals <- paper %>% count(TIER, name = "tier_total")

sensitivity_summary <- matched_df %>%
  group_by(paper_TIER, Comparison) %>%
  summarise(
    n_detected = n(),
    n_gene_length = sum(match_type == "gene+length"),
    n_gene_only = sum(match_type == "gene_only"),
    mean_dPSI = mean(DeltaPSI, na.rm = TRUE),
    median_dPSI = median(DeltaPSI, na.rm = TRUE),
    mean_abs_dPSI = mean(abs(DeltaPSI), na.rm = TRUE),
    pct_negative_dPSI = mean(DeltaPSI < 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  left_join(tier_totals, by = c("paper_TIER" = "TIER")) %>%
  mutate(detection_rate_pct = round(n_detected / tier_total * 100, 1)) %>%
  arrange(paper_TIER, Comparison)

# 5b. dPSI summary per TIER × comparison
tier_dpsi <- matched_df %>%
  group_by(paper_TIER, Comparison) %>%
  summarise(
    n = n(),
    mean_dPSI = round(mean(DeltaPSI, na.rm = TRUE), 4),
    median_dPSI = round(median(DeltaPSI, na.rm = TRUE), 4),
    sd_dPSI = round(sd(DeltaPSI, na.rm = TRUE), 4),
    min_dPSI = round(min(DeltaPSI, na.rm = TRUE), 4),
    max_dPSI = round(max(DeltaPSI, na.rm = TRUE), 4),
    .groups = "drop"
  )

# =============================================================================
# 6. Write output files
# =============================================================================
cat("--- Writing output files ---\n")

# Matched events
write.csv(
  matched_df %>%
    select(paper_EVENT, paper_GENE, paper_LENGTH, paper_TIER, match_type,
           Comparison, EventType, GeneSymbol, exon_length, DeltaPSI, FDR,
           Chromosome, Strand, exonStart_0base, exonEnd, n_candidates),
  file.path(out_dir, "matched_events_by_comparison.csv"),
  row.names = FALSE
)

# Unmatched events
write.csv(unmatched, file.path(out_dir, "unmatched_paper_events.csv"),
          row.names = FALSE)

# Sensitivity summary
write.csv(sensitivity_summary,
          file.path(out_dir, "sensitivity_summary.csv"),
          row.names = FALSE)

# dPSI summary
write.csv(tier_dpsi, file.path(out_dir, "tier_dpsi_summary.csv"),
          row.names = FALSE)

cat("  Written: matched_events_by_comparison.csv\n")
cat("  Written: unmatched_paper_events.csv\n")
cat("  Written: sensitivity_summary.csv\n")
cat("  Written: tier_dpsi_summary.csv\n")

# =============================================================================
# 7. Plots
# =============================================================================
cat("\n--- Generating plots ---\n")

# Define tier ordering and colors
tier_order <- c("HS", "LS", "CS", "CR", "NonResponding", "Others")
tier_colors <- c(
  "HS" = "#D32F2F",
  "LS" = "#F57C00",
  "CS" = "#388E3C",
  "CR" = "#1976D2",
  "NonResponding" = "#757575",
  "Others" = "#9E9E9E"
)

# 7a. Detection rate bar plot
plot_data_bar <- sensitivity_summary %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order))

p_detection <- ggplot(plot_data_bar,
                      aes(x = paper_TIER, y = detection_rate_pct, fill = paper_TIER)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(n_detected, "/", tier_total)),
            vjust = -0.5, size = 3) +
  facet_wrap(~Comparison, nrow = 1) +
  scale_fill_manual(values = tier_colors, guide = "none") +
  labs(
    title = "Microexon Detection Rate by Sensitivity Tier",
    subtitle = "Significant splicing events (FDR<0.05, |dPSI|>0.1) matching paper microexons",
    x = "Sensitivity Tier (Paper Classification)",
    y = "Detection Rate (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "detection_rate_by_tier.pdf"),
       p_detection, width = 12, height = 6)
cat("  Saved: detection_rate_by_tier.pdf\n")

# 7b. dPSI distribution by tier (box + jitter)
plot_data_dpsi <- matched_df %>%
  mutate(paper_TIER = factor(paper_TIER, levels = tier_order))

if (nrow(plot_data_dpsi) > 0) {
  p_dpsi <- ggplot(plot_data_dpsi,
                   aes(x = paper_TIER, y = DeltaPSI, fill = paper_TIER)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~Comparison, nrow = 1) +
    scale_fill_manual(values = tier_colors, guide = "none") +
    labs(
      title = "Delta PSI Distribution by Sensitivity Tier",
      subtitle = "Each point is a matched significant splicing event",
      x = "Sensitivity Tier",
      y = "Delta PSI (rMATS)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(out_dir, "dpsi_distribution_by_tier.pdf"),
         p_dpsi, width = 12, height = 6)
  cat("  Saved: dpsi_distribution_by_tier.pdf\n")
}

# 7c. Heatmap: detection rate — TIER × comparison
heatmap_data <- sensitivity_summary %>%
  select(paper_TIER, Comparison, detection_rate_pct) %>%
  mutate(paper_TIER = factor(paper_TIER, levels = rev(tier_order)),
         Comparison = factor(Comparison, levels = comparisons))

if (nrow(heatmap_data) > 0) {
  p_heat <- ggplot(heatmap_data,
                   aes(x = Comparison, y = paper_TIER, fill = detection_rate_pct)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(detection_rate_pct, 0), "%")),
              size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "white", mid = "#FFCC80", high = "#D32F2F",
                         midpoint = 25, name = "Detection\nRate (%)") +
    labs(
      title = "Microexon Detection Rate: Tier × Comparison",
      x = "Comparison",
      y = "Sensitivity Tier"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )

  ggsave(file.path(out_dir, "detection_heatmap.pdf"),
         p_heat, width = 8, height = 5)
  cat("  Saved: detection_heatmap.pdf\n")
}

# =============================================================================
# 8. Print final summary
# =============================================================================
cat("\n=== FINAL SUMMARY ===\n")
cat("Paper microexons:", nrow(paper), "\n")
cat("Matched in at least one comparison:", length(matched_events),
    "(", round(length(matched_events) / nrow(paper) * 100, 1), "%)\n")
cat("Unmatched:", nrow(unmatched), "\n\n")

cat("Detection rates by tier (across all comparisons):\n")
overall_by_tier <- matched_df %>%
  distinct(paper_EVENT, Comparison, .keep_all = TRUE) %>%
  group_by(paper_TIER) %>%
  summarise(
    unique_events_detected = n_distinct(paper_EVENT),
    .groups = "drop"
  ) %>%
  left_join(tier_totals, by = c("paper_TIER" = "TIER")) %>%
  mutate(overall_pct = round(unique_events_detected / tier_total * 100, 1))

print(as.data.frame(overall_by_tier))

cat("\nDone:", format(Sys.time()), "\n")
