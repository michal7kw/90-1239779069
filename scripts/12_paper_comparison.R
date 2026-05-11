#!/usr/bin/env Rscript
# =============================================================================
# 12_paper_comparison.R
# =============================================================================
# Compare SRRM3 RNA-seq splicing data with Gonatopoulos-Pournatzis 2018
# paper (Table S3 SPAR-seq data)
#
# Extended to analyze all three vs-Parental comparisons:
#   - Neg_vs_Parental (partial SRRM3 knockdown)
#   - Pos_vs_Parental (SRRM3 overexpression)
#   - KO_vs_Parental  (SRRM3 knockout)
#
# Input:
#   - data/PSI.csv (dPSI values from paper)
#   - data/event-information.csv (coordinates in mm9)
#   - results/05_splicing/{comparison}/SE.MATS.JC.txt (rMATS output in mm39)
#   - data/mm9ToMm39.over.chain.gz (liftOver chain file)
#
# Output:
#   - results/12_paper_comparison/{comparison}/
# =============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

# =============================================================================
# Configuration
# =============================================================================

base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"

# Input files (static - same for all comparisons)
psi_file <- file.path(base_dir, "data/PSI.csv")
event_info_file <- file.path(base_dir, "data/event-information.csv")
chain_file <- file.path(base_dir, "data/mm9ToMm39.over.chain")

# Base output directory
base_output_dir <- file.path(base_dir, "results/12_paper_comparison")
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# Comparisons to analyze
comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")

# Condition labels for plot annotations
condition_labels <- list(
  Neg_vs_Parental = list(
    short = "Neg",
    full = "Partial KD",
    description = "Partial SRRM3 Knockdown"
  ),
  Pos_vs_Parental = list(
    short = "Pos",
    full = "Overexpression",
    description = "SRRM3 Overexpression"
  ),
  KO_vs_Parental = list(
    short = "KO",
    full = "Knockout",
    description = "SRRM3 Knockout"
  )
)

# Parameters
coordinate_tolerance <- 10  # bp tolerance for coordinate matching

# =============================================================================
# Helper Functions
# =============================================================================

#' Load paper data (PSI values and event coordinates)
#' @return list with event_info, srrm3_psi, srrm4_psi, and merged paper_data
load_paper_data <- function() {
  cat("Loading paper data...\n")

  # Load event coordinates (mm9)
  event_info <- read.csv(event_info_file, stringsAsFactors = FALSE)
  cat(sprintf("  - Loaded %d microexon events from paper\n", nrow(event_info)))

  # Load PSI values (transposed matrix)
  psi_raw <- read.csv(psi_file, stringsAsFactors = FALSE)
  colnames(psi_raw)[1] <- "knockdown"

  # Reshape PSI data to long format
  psi_long <- psi_raw %>%
    pivot_longer(
      cols = -knockdown,
      names_to = "event",
      values_to = "dPSI_percent"
    ) %>%
    mutate(
      # Convert from percentage to 0-1 scale (paper uses percent, rMATS uses 0-1)
      dPSI = dPSI_percent / 100
    )

  # Extract Srrm3 and Srrm4 knockdown data
  srrm3_psi <- psi_long %>%
    filter(knockdown == "Srrm3") %>%
    select(event, srrm3_dPSI = dPSI, srrm3_dPSI_percent = dPSI_percent)

  srrm4_psi <- psi_long %>%
    filter(knockdown == "Srrm4") %>%
    select(event, srrm4_dPSI = dPSI, srrm4_dPSI_percent = dPSI_percent)

  cat(sprintf("  - Srrm3 knockdown: %d events with dPSI values\n",
              sum(!is.na(srrm3_psi$srrm3_dPSI))))
  cat(sprintf("  - Srrm4 knockdown: %d events with dPSI values\n",
              sum(!is.na(srrm4_psi$srrm4_dPSI))))

  # Merge coordinates with dPSI values
  paper_data <- event_info %>%
    left_join(srrm3_psi, by = "event") %>%
    left_join(srrm4_psi, by = "event")

  cat(sprintf("  - Combined paper data: %d events\n", nrow(paper_data)))

  list(
    event_info = event_info,
    srrm3_psi = srrm3_psi,
    srrm4_psi = srrm4_psi,
    paper_data = paper_data
  )
}

#' Perform liftOver from mm9 to mm39 coordinates
#' @param paper_data Data frame with mm9 coordinates
#' @return Data frame with mm39 coordinates
liftover_coordinates <- function(paper_data) {
  cat("LiftOver mm9 -> mm39...\n")

  # Load chain file
  chain <- import.chain(chain_file)
  cat("  - Loaded chain file\n")

  # Create GRanges object from paper data (mm9)
  paper_gr_mm9 <- GRanges(
    seqnames = paper_data$chrom,
    ranges = IRanges(start = paper_data$start, end = paper_data$end),
    strand = paper_data$strand,
    event = paper_data$event,
    gene = paper_data$gene,
    VastDB.ID = paper_data$VastDB.ID,
    srrm3_dPSI = paper_data$srrm3_dPSI,
    srrm4_dPSI = paper_data$srrm4_dPSI
  )

  # Perform liftOver
  paper_gr_mm39_list <- liftOver(paper_gr_mm9, chain)

  # Count successful conversions
  n_converted <- sum(lengths(paper_gr_mm39_list) > 0)
  cat(sprintf("  - Successfully converted: %d / %d events (%.1f%%)\n",
              n_converted, length(paper_gr_mm9), 100 * n_converted / length(paper_gr_mm9)))

  # Extract successfully converted coordinates
  paper_gr_mm39 <- unlist(paper_gr_mm39_list[lengths(paper_gr_mm39_list) > 0])

  # Track failed conversions
  failed_events <- paper_data$event[lengths(paper_gr_mm39_list) == 0]
  if (length(failed_events) > 0) {
    cat(sprintf("  - Failed liftOver for: %s\n", paste(failed_events, collapse = ", ")))
  }

  # Create data frame with mm39 coordinates
  paper_mm39 <- data.frame(
    event = paper_gr_mm39$event,
    gene = paper_gr_mm39$gene,
    VastDB.ID = paper_gr_mm39$VastDB.ID,
    chrom_mm39 = as.character(seqnames(paper_gr_mm39)),
    strand = as.character(strand(paper_gr_mm39)),
    start_mm39 = start(paper_gr_mm39),
    end_mm39 = end(paper_gr_mm39),
    exon_size = width(paper_gr_mm39),
    srrm3_dPSI = paper_gr_mm39$srrm3_dPSI,
    srrm4_dPSI = paper_gr_mm39$srrm4_dPSI,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  - Paper events after liftOver: %d\n", nrow(paper_mm39)))

  paper_mm39
}

#' Load rMATS data for a specific comparison
#' @param comparison Comparison name (e.g., "KO_vs_Parental")
#' @return Data frame with rMATS results
load_rmats_data <- function(comparison) {
  rmats_file <- file.path(base_dir, "results/05_splicing", comparison, "SE.MATS.JC.txt")

  if (!file.exists(rmats_file)) {
    stop(sprintf("rMATS file not found: %s", rmats_file))
  }

  # read.delim with check.names = TRUE will rename duplicate columns (ID -> ID.1)
  rmats <- read.delim(rmats_file, stringsAsFactors = FALSE, check.names = TRUE)

  # The first ID column becomes rmats_ID (the second becomes ID.1 if duplicate)
  rmats <- rmats %>%
    rename(rmats_ID = ID) %>%
    mutate(
      geneSymbol = gsub('"', '', geneSymbol),
      GeneID = gsub('"', '', GeneID),
      exon_size = exonEnd - exonStart_0base,
      exon_start_1based = exonStart_0base + 1
    )

  cat(sprintf("  - Total SE events in rMATS: %d\n", nrow(rmats)))
  cat(sprintf("  - Events with FDR < 0.05: %d\n", sum(rmats$FDR < 0.05, na.rm = TRUE)))
  cat(sprintf("  - Events with |dPSI| > 0.1 and FDR < 0.05: %d\n",
              sum(rmats$FDR < 0.05 & abs(rmats$IncLevelDifference) > 0.1, na.rm = TRUE)))

  rmats
}

#' Find coordinate match in rMATS data
#' @param rmats rMATS data frame
#' @param chrom Chromosome
#' @param start Start position
#' @param end End position
#' @param strand Strand
#' @param tolerance Coordinate tolerance in bp
#' @return List with 'best' (single best match) and 'all' (all candidates), or NULL
find_coordinate_match <- function(rmats, chrom, start, end, strand, tolerance = 10) {
  candidates <- rmats %>%
    filter(chr == chrom, strand == !!strand)

  if (nrow(candidates) == 0) return(NULL)

  matches <- candidates %>%
    filter(
      abs(exon_start_1based - start) <= tolerance,
      abs(exonEnd - end) <= tolerance
    )

  if (nrow(matches) == 0) return(NULL)

  all_matches <- matches

  if (nrow(matches) > 1) {
    # Prefer significant events (FDR < 0.05), then strongest |dPSI|
    sig_matches <- matches %>% filter(FDR < 0.05)
    if (nrow(sig_matches) > 0) {
      matches <- sig_matches %>%
        arrange(desc(abs(IncLevelDifference))) %>%
        slice(1)
    } else {
      matches <- matches %>%
        arrange(desc(abs(IncLevelDifference))) %>%
        slice(1)
    }
  }

  return(list(best = matches, all = all_matches))
}

#' Find gene name match in rMATS data
#' @param rmats rMATS data frame
#' @param gene_name Gene name
#' @param target_exon_size Target exon size
#' @param size_tolerance Size tolerance in bp
#' @return List with 'best' (single best match) and 'all' (all candidates), or NULL
find_gene_match <- function(rmats, gene_name, target_exon_size, size_tolerance = 5) {
  candidates <- rmats %>%
    filter(toupper(geneSymbol) == toupper(gene_name))

  if (nrow(candidates) == 0) return(NULL)

  matches <- candidates %>%
    filter(abs(exon_size - target_exon_size) <= size_tolerance)

  if (nrow(matches) == 0) {
    matches <- candidates %>%
      mutate(size_diff = abs(exon_size - target_exon_size)) %>%
      arrange(size_diff) %>%
      slice(1) %>%
      select(-size_diff)
  }

  all_matches <- matches

  if (nrow(matches) > 1) {
    # Prefer significant events (FDR < 0.05), then strongest |dPSI|
    sig_matches <- matches %>% filter(FDR < 0.05)
    if (nrow(sig_matches) > 0) {
      matches <- sig_matches %>%
        arrange(desc(abs(IncLevelDifference))) %>%
        slice(1)
    } else {
      matches <- matches %>%
        arrange(desc(abs(IncLevelDifference))) %>%
        slice(1)
    }
  }

  return(list(best = matches, all = all_matches))
}

#' Match paper events to rMATS events
#' @param paper_mm39 Paper data with mm39 coordinates
#' @param rmats rMATS data frame
#' @param tolerance Coordinate tolerance in bp
#' @return List with 'comparison_df' and 'multi_match_records'
match_events <- function(paper_mm39, rmats, tolerance = 10) {
  match_results <- list()
  multi_match_records <- list()

  for (i in 1:nrow(paper_mm39)) {
    event_data <- paper_mm39[i, ]

    # Try coordinate matching first
    coord_result <- find_coordinate_match(
      rmats = rmats,
      chrom = event_data$chrom_mm39,
      start = event_data$start_mm39,
      end = event_data$end_mm39,
      strand = event_data$strand,
      tolerance = tolerance
    )

    if (!is.null(coord_result)) {
      best <- coord_result$best
      match_results[[i]] <- list(
        paper_event = event_data$event,
        match_type = "coordinate",
        n_rmats_candidates = nrow(coord_result$all),
        rmats_id = best$rmats_ID[1],
        rmats_gene = best$geneSymbol[1],
        rmats_chr = best$chr[1],
        rmats_start = best$exon_start_1based[1],
        rmats_end = best$exonEnd[1],
        rmats_exon_size = best$exon_size[1],
        rmats_dPSI = best$IncLevelDifference[1],
        rmats_FDR = best$FDR[1],
        rmats_PValue = best$PValue[1]
      )

      # Record multi-match events
      if (nrow(coord_result$all) > 1) {
        for (j in 1:nrow(coord_result$all)) {
          row <- coord_result$all[j, ]
          multi_match_records[[length(multi_match_records) + 1]] <- data.frame(
            paper_event = event_data$event,
            paper_gene = event_data$gene,
            match_type = "coordinate",
            n_candidates = nrow(coord_result$all),
            selected = (row$rmats_ID == best$rmats_ID[1]),
            rmats_ID = row$rmats_ID,
            geneSymbol = row$geneSymbol,
            chr = row$chr,
            exonStart_0base = row$exonStart_0base,
            exonEnd = row$exonEnd,
            strand = row$strand,
            upstreamES = row$upstreamES,
            upstreamEE = row$upstreamEE,
            downstreamES = row$downstreamES,
            downstreamEE = row$downstreamEE,
            IncLevelDifference = row$IncLevelDifference,
            FDR = row$FDR,
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      # Try gene name matching as fallback
      gene_result <- find_gene_match(
        rmats = rmats,
        gene_name = event_data$gene,
        target_exon_size = event_data$exon_size
      )

      if (!is.null(gene_result)) {
        best <- gene_result$best
        match_results[[i]] <- list(
          paper_event = event_data$event,
          match_type = "gene_name",
          n_rmats_candidates = nrow(gene_result$all),
          rmats_id = best$rmats_ID[1],
          rmats_gene = best$geneSymbol[1],
          rmats_chr = best$chr[1],
          rmats_start = best$exon_start_1based[1],
          rmats_end = best$exonEnd[1],
          rmats_exon_size = best$exon_size[1],
          rmats_dPSI = best$IncLevelDifference[1],
          rmats_FDR = best$FDR[1],
          rmats_PValue = best$PValue[1]
        )

        # Record multi-match events
        if (nrow(gene_result$all) > 1) {
          for (j in 1:nrow(gene_result$all)) {
            row <- gene_result$all[j, ]
            multi_match_records[[length(multi_match_records) + 1]] <- data.frame(
              paper_event = event_data$event,
              paper_gene = event_data$gene,
              match_type = "gene_name",
              n_candidates = nrow(gene_result$all),
              selected = (row$rmats_ID == best$rmats_ID[1]),
              rmats_ID = row$rmats_ID,
              geneSymbol = row$geneSymbol,
              chr = row$chr,
              exonStart_0base = row$exonStart_0base,
              exonEnd = row$exonEnd,
              strand = row$strand,
              upstreamES = row$upstreamES,
              upstreamEE = row$upstreamEE,
              downstreamES = row$downstreamES,
              downstreamEE = row$downstreamEE,
              IncLevelDifference = row$IncLevelDifference,
              FDR = row$FDR,
              stringsAsFactors = FALSE
            )
          }
        }
      } else {
        match_results[[i]] <- list(
          paper_event = event_data$event,
          match_type = "unmatched",
          n_rmats_candidates = 0,
          rmats_id = NA,
          rmats_gene = NA,
          rmats_chr = NA,
          rmats_start = NA,
          rmats_end = NA,
          rmats_exon_size = NA,
          rmats_dPSI = NA,
          rmats_FDR = NA,
          rmats_PValue = NA
        )
      }
    }
  }

  # Convert to data frame
  match_df <- bind_rows(lapply(match_results, as.data.frame))

  # Merge with paper data
  comparison_df <- paper_mm39 %>%
    left_join(match_df, by = c("event" = "paper_event"))

  # Combine multi-match records
  if (length(multi_match_records) > 0) {
    multi_match_df <- bind_rows(multi_match_records)
  } else {
    multi_match_df <- data.frame()
  }

  list(comparison_df = comparison_df, multi_match_records = multi_match_df)
}

#' Calculate correlation statistics
#' @param comparison_df Comparison data frame with matched events
#' @param condition_short Short condition name (e.g., "KO")
#' @return List with statistics for Srrm3 and Srrm4 comparisons
calculate_statistics <- function(comparison_df, condition_short) {
  matched_for_stats <- comparison_df %>%
    filter(match_type != "unmatched", !is.na(rmats_dPSI))

  results <- list()

  # Srrm3 comparison
  srrm3_stats <- matched_for_stats %>%
    filter(!is.na(srrm3_dPSI))

  if (nrow(srrm3_stats) >= 3) {
    srrm3_pearson <- cor.test(srrm3_stats$srrm3_dPSI, srrm3_stats$rmats_dPSI,
                              method = "pearson")
    srrm3_spearman <- cor.test(srrm3_stats$srrm3_dPSI, srrm3_stats$rmats_dPSI,
                               method = "spearman")
    srrm3_concordance <- mean(sign(srrm3_stats$srrm3_dPSI) == sign(srrm3_stats$rmats_dPSI),
                              na.rm = TRUE) * 100

    results$srrm3 <- list(
      data = srrm3_stats,
      n = nrow(srrm3_stats),
      pearson = srrm3_pearson,
      spearman = srrm3_spearman,
      concordance = srrm3_concordance
    )
  } else {
    results$srrm3 <- list(
      data = srrm3_stats,
      n = nrow(srrm3_stats),
      pearson = list(estimate = NA, p.value = NA),
      spearman = list(estimate = NA, p.value = NA),
      concordance = NA
    )
  }

  # Srrm4 comparison
  srrm4_stats <- matched_for_stats %>%
    filter(!is.na(srrm4_dPSI))

  if (nrow(srrm4_stats) >= 3) {
    srrm4_pearson <- cor.test(srrm4_stats$srrm4_dPSI, srrm4_stats$rmats_dPSI,
                              method = "pearson")
    srrm4_spearman <- cor.test(srrm4_stats$srrm4_dPSI, srrm4_stats$rmats_dPSI,
                               method = "spearman")
    srrm4_concordance <- mean(sign(srrm4_stats$srrm4_dPSI) == sign(srrm4_stats$rmats_dPSI),
                              na.rm = TRUE) * 100

    results$srrm4 <- list(
      data = srrm4_stats,
      n = nrow(srrm4_stats),
      pearson = srrm4_pearson,
      spearman = srrm4_spearman,
      concordance = srrm4_concordance
    )
  } else {
    results$srrm4 <- list(
      data = srrm4_stats,
      n = nrow(srrm4_stats),
      pearson = list(estimate = NA, p.value = NA),
      spearman = list(estimate = NA, p.value = NA),
      concordance = NA
    )
  }

  # Create statistics summary data frame
  results$summary <- data.frame(
    comparison = c(
      sprintf("Srrm3_KD_vs_SRRM3_%s", condition_short),
      sprintf("Srrm4_KD_vs_SRRM3_%s", condition_short)
    ),
    n_events = c(results$srrm3$n, results$srrm4$n),
    pearson_r = c(results$srrm3$pearson$estimate, results$srrm4$pearson$estimate),
    pearson_p = c(results$srrm3$pearson$p.value, results$srrm4$pearson$p.value),
    spearman_rho = c(results$srrm3$spearman$estimate, results$srrm4$spearman$estimate),
    spearman_p = c(results$srrm3$spearman$p.value, results$srrm4$spearman$p.value),
    direction_concordance_pct = c(results$srrm3$concordance, results$srrm4$concordance)
  )

  results
}

#' Generate all visualizations for a comparison
#' @param comparison_df Comparison data frame
#' @param stats Statistics from calculate_statistics()
#' @param comparison Comparison name
#' @param condition_label Condition label list
#' @param output_dir Output directory
generate_plots <- function(comparison_df, stats, comparison, condition_label, output_dir) {

  # Theme for plots
  theme_publication <- theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  condition_short <- condition_label$short
  condition_desc <- condition_label$description

  plots <- list()

  # Scatter plot: Srrm3 KD vs this study
  if (stats$srrm3$n >= 3 && !is.na(stats$srrm3$pearson$estimate)) {
    srrm3_data <- stats$srrm3$data

    p_srrm3 <- ggplot(srrm3_data, aes(x = srrm3_dPSI, y = rmats_dPSI)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue", fill = "lightblue", alpha = 0.3) +
      geom_point(aes(color = match_type), size = 3, alpha = 0.8) +
      geom_text(aes(label = gene), hjust = -0.1, vjust = -0.5, size = 2.5) +
      scale_color_manual(values = c("coordinate" = "#2166AC", "gene_name" = "#B2182B")) +
      labs(
        title = sprintf("Srrm3 Knockdown (Paper) vs %s (This Study)", condition_desc),
        subtitle = sprintf("Pearson r = %.3f, p = %.4f | Direction concordance = %.1f%%",
                          stats$srrm3$pearson$estimate,
                          stats$srrm3$pearson$p.value,
                          stats$srrm3$concordance),
        x = "Paper Srrm3 KD dPSI",
        y = sprintf("This Study %s dPSI", condition_short),
        color = "Match Type"
      ) +
      coord_fixed(ratio = 1) +
      xlim(c(-1, 0.5)) +
      ylim(c(-1, 0.5)) +
      theme_publication

    ggsave(file.path(output_dir, sprintf("scatter_srrm3_vs_%s.pdf", condition_short)),
           p_srrm3, width = 8, height = 7)
    cat(sprintf("  - Saved: scatter_srrm3_vs_%s.pdf\n", condition_short))

    plots$srrm3 <- p_srrm3
  } else {
    cat(sprintf("  - Skipped Srrm3 scatter plot (insufficient data: n=%d)\n", stats$srrm3$n))
    plots$srrm3 <- NULL
  }

  # Scatter plot: Srrm4 KD vs this study
  if (stats$srrm4$n >= 3 && !is.na(stats$srrm4$pearson$estimate)) {
    srrm4_data <- stats$srrm4$data

    p_srrm4 <- ggplot(srrm4_data, aes(x = srrm4_dPSI, y = rmats_dPSI)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      geom_smooth(method = "lm", se = TRUE, color = "#B2182B", fill = "#FDDBC7", alpha = 0.3) +
      geom_point(aes(color = match_type), size = 3, alpha = 0.8) +
      geom_text(aes(label = gene), hjust = -0.1, vjust = -0.5, size = 2.5) +
      scale_color_manual(values = c("coordinate" = "#2166AC", "gene_name" = "#B2182B")) +
      labs(
        title = sprintf("Srrm4 Knockdown (Paper) vs %s (This Study)", condition_desc),
        subtitle = sprintf("Pearson r = %.3f, p = %.4f | Direction concordance = %.1f%%",
                          stats$srrm4$pearson$estimate,
                          stats$srrm4$pearson$p.value,
                          stats$srrm4$concordance),
        x = "Paper Srrm4 KD dPSI",
        y = sprintf("This Study %s dPSI", condition_short),
        color = "Match Type"
      ) +
      coord_fixed(ratio = 1) +
      xlim(c(-1, 0.5)) +
      ylim(c(-1, 0.5)) +
      theme_publication

    ggsave(file.path(output_dir, sprintf("scatter_srrm4_vs_%s.pdf", condition_short)),
           p_srrm4, width = 8, height = 7)
    cat(sprintf("  - Saved: scatter_srrm4_vs_%s.pdf\n", condition_short))

    plots$srrm4 <- p_srrm4
  } else {
    cat(sprintf("  - Skipped Srrm4 scatter plot (insufficient data: n=%d)\n", stats$srrm4$n))
    plots$srrm4 <- NULL
  }

  # Combined panel plot
  if (!is.null(plots$srrm3) && !is.null(plots$srrm4)) {
    pdf(file.path(output_dir, "scatter_combined.pdf"), width = 14, height = 7)
    grid.arrange(plots$srrm3, plots$srrm4, ncol = 2,
                 top = sprintf("Comparison with Gonatopoulos-Pournatzis et al. 2018 - %s", comparison))
    dev.off()
    cat("  - Saved: scatter_combined.pdf\n")
  }

  # Heatmap of dPSI values
  matched_for_stats <- comparison_df %>%
    filter(match_type != "unmatched", !is.na(rmats_dPSI))

  heatmap_data <- matched_for_stats %>%
    filter(!is.na(srrm3_dPSI) | !is.na(srrm4_dPSI)) %>%
    select(gene, srrm3_dPSI, srrm4_dPSI, rmats_dPSI) %>%
    rename(
      `Paper Srrm3 KD` = srrm3_dPSI,
      `Paper Srrm4 KD` = srrm4_dPSI,
      !!sprintf("This Study %s", condition_short) := rmats_dPSI
    ) %>%
    pivot_longer(cols = -gene, names_to = "Condition", values_to = "dPSI") %>%
    mutate(Condition = factor(Condition, levels = c("Paper Srrm3 KD", "Paper Srrm4 KD",
                                                     sprintf("This Study %s", condition_short))))

  if (nrow(heatmap_data) > 0) {
    p_heatmap <- ggplot(heatmap_data, aes(x = Condition, y = gene, fill = dPSI)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", dPSI)), size = 2.5, color = "black") +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-1, 1),
        name = "dPSI"
      ) +
      labs(
        title = sprintf("dPSI Values Across Conditions - %s", comparison),
        subtitle = "Comparing paper knockdowns with this study",
        x = "",
        y = "Gene"
      ) +
      theme_publication +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(output_dir, "dpsi_heatmap.pdf"), p_heatmap,
           width = 6, height = max(4, nrow(matched_for_stats) * 0.3))
    cat("  - Saved: dpsi_heatmap.pdf\n")
  }

  # Summary bar plot
  match_df <- comparison_df %>%
    select(event, match_type) %>%
    distinct()

  match_summary <- match_df %>%
    count(match_type) %>%
    mutate(
      match_type = factor(match_type, levels = c("coordinate", "gene_name", "unmatched"))
    )

  p_summary <- ggplot(match_summary, aes(x = match_type, y = n, fill = match_type)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("coordinate" = "#2166AC", "gene_name" = "#F4A582", "unmatched" = "#D1D1D1")) +
    labs(
      title = sprintf("Event Matching Summary - %s", comparison),
      subtitle = sprintf("Total paper events: %d", nrow(match_df)),
      x = "Match Type",
      y = "Number of Events"
    ) +
    theme_publication +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("Coordinate\nMatch", "Gene Name\nMatch", "Unmatched"))

  ggsave(file.path(output_dir, "summary_bar.pdf"), p_summary,
         width = 5, height = 4)
  cat("  - Saved: summary_bar.pdf\n")

  plots
}

#' Save output files for a comparison
#' @param comparison_df Comparison data frame
#' @param stats Statistics from calculate_statistics()
#' @param multi_match_df Data frame with multi-match details
#' @param comparison Comparison name
#' @param output_dir Output directory
save_outputs <- function(comparison_df, stats, multi_match_df, comparison, output_dir) {
  # Full comparison table
  write.csv(comparison_df, file.path(output_dir, "paper_comparison_full.csv"), row.names = FALSE)
  cat("  - Saved: paper_comparison_full.csv\n")

  # Matched events only
  matched_events <- comparison_df %>%
    filter(match_type != "unmatched")
  write.csv(matched_events, file.path(output_dir, "matched_events.csv"), row.names = FALSE)
  cat("  - Saved: matched_events.csv\n")

  # Unmatched events
  unmatched_events <- comparison_df %>%
    filter(match_type == "unmatched")
  write.csv(unmatched_events, file.path(output_dir, "unmatched_events.csv"), row.names = FALSE)
  cat("  - Saved: unmatched_events.csv\n")

  # Statistics summary
  write.csv(stats$summary, file.path(output_dir, "statistics_summary.csv"), row.names = FALSE)
  cat("  - Saved: statistics_summary.csv\n")

  # Multi-match report
  report_file <- file.path(output_dir, "multi_match_events.txt")
  if (nrow(multi_match_df) > 0) {
    n_events <- length(unique(multi_match_df$paper_event))
    con <- file(report_file, "w")
    writeLines(sprintf("Multi-Match Event Report: %s", comparison), con)
    writeLines(sprintf("Generated: %s", Sys.time()), con)
    writeLines(paste0(rep("=", 80), collapse = ""), con)
    writeLines("", con)
    writeLines(sprintf("Paper events with multiple rMATS candidates: %d", n_events), con)
    writeLines(sprintf("Total rMATS candidate entries: %d", nrow(multi_match_df)), con)
    writeLines("Selection criterion: prefer significant (FDR < 0.05), then strongest |dPSI|", con)
    writeLines("", con)
    writeLines(paste0(rep("-", 80), collapse = ""), con)

    for (evt in unique(multi_match_df$paper_event)) {
      evt_rows <- multi_match_df %>% filter(paper_event == evt)
      gene <- evt_rows$paper_gene[1]
      n_cand <- evt_rows$n_candidates[1]
      mtype <- evt_rows$match_type[1]

      writeLines("", con)
      writeLines(sprintf("Paper event: %s  |  Gene: %s  |  Match type: %s  |  Candidates: %d",
                         evt, gene, mtype, n_cand), con)

      for (k in 1:nrow(evt_rows)) {
        r <- evt_rows[k, ]
        tag <- if (r$selected) {
          if (r$FDR < 0.05) " << SELECTED (significant, strongest |dPSI|)"
          else " << SELECTED (none significant, strongest |dPSI|)"
        } else ""
        writeLines(sprintf(
          "  rMATS ID=%d  %s:%d-%d (%s)  up:%d-%d  down:%d-%d  dPSI=%.4f  FDR=%.2e%s",
          r$rmats_ID, r$chr, r$exonStart_0base, r$exonEnd, r$strand,
          r$upstreamES, r$upstreamEE, r$downstreamES, r$downstreamEE,
          r$IncLevelDifference, r$FDR, tag
        ), con)
      }
    }

    writeLines("", con)
    writeLines(paste0(rep("=", 80), collapse = ""), con)
    close(con)
    cat(sprintf("  - Saved: multi_match_events.txt (%d events with %d total candidates)\n",
                n_events, nrow(multi_match_df)))
  } else {
    con <- file(report_file, "w")
    writeLines(sprintf("Multi-Match Event Report: %s", comparison), con)
    writeLines(sprintf("Generated: %s", Sys.time()), con)
    writeLines(paste0(rep("=", 80), collapse = ""), con)
    writeLines("", con)
    writeLines("No multi-match events found. All paper events matched to a single rMATS event.", con)
    close(con)
    cat("  - Saved: multi_match_events.txt (no multi-matches)\n")
  }

  # Also save multi-match data as CSV for programmatic access
  if (nrow(multi_match_df) > 0) {
    write.csv(multi_match_df, file.path(output_dir, "multi_match_events.csv"), row.names = FALSE)
    cat("  - Saved: multi_match_events.csv\n")
  }
}

#' Process a single comparison
#' @param comparison Comparison name
#' @param paper_mm39 Paper data with mm39 coordinates
#' @param condition_label Condition label list
#' @return List with results
process_comparison <- function(comparison, paper_mm39, condition_label) {
  cat("\n")
  cat("=============================================================================\n")
  cat(sprintf("Processing: %s (%s)\n", comparison, condition_label$description))
  cat("=============================================================================\n\n")

  # Create output directory
  output_dir <- file.path(base_output_dir, comparison)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load rMATS data
  cat("Loading rMATS data...\n")
  rmats <- load_rmats_data(comparison)

  # Match events
  cat("\nMatching events...\n")
  match_result <- match_events(paper_mm39, rmats, tolerance = coordinate_tolerance)
  comparison_df <- match_result$comparison_df
  multi_match_df <- match_result$multi_match_records

  # Print matching summary
  n_coord <- sum(comparison_df$match_type == "coordinate")
  n_gene <- sum(comparison_df$match_type == "gene_name")
  n_unmatched <- sum(comparison_df$match_type == "unmatched")
  n_total <- nrow(comparison_df)

  cat(sprintf("  - Coordinate matches: %d\n", n_coord))
  cat(sprintf("  - Gene name matches: %d\n", n_gene))
  cat(sprintf("  - Unmatched: %d\n", n_unmatched))
  cat(sprintf("  - Total matched: %d / %d (%.1f%%)\n",
              n_coord + n_gene, n_total, 100 * (n_coord + n_gene) / n_total))

  # Report multi-match events
  if (nrow(multi_match_df) > 0) {
    n_multi <- length(unique(multi_match_df$paper_event))
    cat(sprintf("  - Events with multiple rMATS candidates: %d (prefer significant, then strongest |dPSI|)\n",
                n_multi))
  }

  # Calculate statistics
  cat("\nCalculating statistics...\n")
  stats <- calculate_statistics(comparison_df, condition_label$short)

  # Print correlation summary
  cat(sprintf("\n  Srrm3 KD vs %s:\n", condition_label$short))
  if (!is.na(stats$srrm3$pearson$estimate)) {
    cat(sprintf("    - N events: %d\n", stats$srrm3$n))
    cat(sprintf("    - Pearson r: %.3f (p = %.4f)\n",
                stats$srrm3$pearson$estimate, stats$srrm3$pearson$p.value))
    cat(sprintf("    - Direction concordance: %.1f%%\n", stats$srrm3$concordance))
  } else {
    cat("    - Insufficient data for correlation\n")
  }

  cat(sprintf("\n  Srrm4 KD vs %s:\n", condition_label$short))
  if (!is.na(stats$srrm4$pearson$estimate)) {
    cat(sprintf("    - N events: %d\n", stats$srrm4$n))
    cat(sprintf("    - Pearson r: %.3f (p = %.4f)\n",
                stats$srrm4$pearson$estimate, stats$srrm4$pearson$p.value))
    cat(sprintf("    - Direction concordance: %.1f%%\n", stats$srrm4$concordance))
  } else {
    cat("    - Insufficient data for correlation\n")
  }

  # Generate plots
  cat("\nGenerating visualizations...\n")
  plots <- generate_plots(comparison_df, stats, comparison, condition_label, output_dir)

  # Save outputs
  cat("\nSaving output files...\n")
  save_outputs(comparison_df, stats, multi_match_df, comparison, output_dir)

  list(
    comparison = comparison,
    comparison_df = comparison_df,
    stats = stats,
    plots = plots,
    rmats = rmats,
    multi_match_df = multi_match_df,
    n_matched = n_coord + n_gene,
    n_total = n_total
  )
}

# =============================================================================
# Main Analysis
# =============================================================================

cat("=============================================================================\n")
cat("SRRM3 Multi-Comparison Paper Validation Analysis\n")
cat("Gonatopoulos-Pournatzis et al. 2018 Comparison\n")
cat("=============================================================================\n\n")

# Step 1: Load paper data (done once)
cat("Step 1: Loading and preparing paper data...\n")
cat("-----------------------------------------------------------------------------\n")
paper <- load_paper_data()

# Step 2: LiftOver coordinates (done once)
cat("\nStep 2: Converting coordinates...\n")
cat("-----------------------------------------------------------------------------\n")
paper_mm39 <- liftover_coordinates(paper$paper_data)

# Step 3: Process each comparison
cat("\nStep 3: Processing comparisons...\n")
cat("-----------------------------------------------------------------------------\n")

all_results <- list()
for (comparison in comparisons) {
  condition_label <- condition_labels[[comparison]]
  all_results[[comparison]] <- process_comparison(comparison, paper_mm39, condition_label)
}

# =============================================================================
# Cross-Comparison Summary
# =============================================================================

cat("\n\n")
cat("=============================================================================\n")
cat("CROSS-COMPARISON SUMMARY\n")
cat("=============================================================================\n\n")

# Combine statistics from all comparisons
all_stats <- bind_rows(lapply(all_results, function(x) {
  x$stats$summary %>%
    mutate(study_comparison = x$comparison)
}))

# Print summary table
cat("Srrm3 KD (Paper) vs This Study:\n")
cat("-----------------------------------------------------------------------------\n")
srrm3_summary <- all_stats %>%
  filter(grepl("Srrm3", comparison)) %>%
  select(study_comparison, n_events, pearson_r, pearson_p, direction_concordance_pct)

for (i in 1:nrow(srrm3_summary)) {
  row <- srrm3_summary[i, ]
  cat(sprintf("  %s:\n", row$study_comparison))
  cat(sprintf("    N = %d, r = %.3f (p = %.4f), concordance = %.1f%%\n",
              row$n_events, row$pearson_r, row$pearson_p, row$direction_concordance_pct))
}

cat("\nSrrm4 KD (Paper) vs This Study:\n")
cat("-----------------------------------------------------------------------------\n")
srrm4_summary <- all_stats %>%
  filter(grepl("Srrm4", comparison)) %>%
  select(study_comparison, n_events, pearson_r, pearson_p, direction_concordance_pct)

for (i in 1:nrow(srrm4_summary)) {
  row <- srrm4_summary[i, ]
  cat(sprintf("  %s:\n", row$study_comparison))
  cat(sprintf("    N = %d, r = %.3f (p = %.4f), concordance = %.1f%%\n",
              row$n_events, row$pearson_r, row$pearson_p, row$direction_concordance_pct))
}

# Save cross-comparison summary
write.csv(all_stats, file.path(base_output_dir, "cross_comparison_statistics.csv"), row.names = FALSE)
cat(sprintf("\nSaved cross-comparison summary: %s/cross_comparison_statistics.csv\n",
            base_output_dir))

# =============================================================================
# Final Summary
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")

cat("Output directories created:\n")
for (comparison in comparisons) {
  output_dir <- file.path(base_output_dir, comparison)
  cat(sprintf("  - %s/\n", output_dir))
  cat("      scatter_srrm3_vs_*.pdf, scatter_srrm4_vs_*.pdf\n")
  cat("      scatter_combined.pdf, dpsi_heatmap.pdf, summary_bar.pdf\n")
  cat("      paper_comparison_full.csv, matched_events.csv\n")
  cat("      unmatched_events.csv, statistics_summary.csv\n")
  cat("      multi_match_events.txt, multi_match_events.csv\n")
}

cat(sprintf("\nCross-comparison summary: %s/cross_comparison_statistics.csv\n", base_output_dir))

cat("\n=============================================================================\n")
