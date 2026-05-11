#!/usr/bin/env Rscript
# =============================================================================
# 12_paper_comparison_double_kd_microexons.R
# =============================================================================
# Compare SRRM3 RNA-seq splicing data with Gonatopoulos-Pournatzis 2018
# DOUBLE KNOCKDOWN (Srrm3+Srrm4) data — MICROEXONS ONLY (≤30bp)
#
# This is the microexon-filtered variant of 12_paper_comparison_double_kd.R.
# SRRM3 is a neural-specific microexon splicing regulator, so filtering to
# exons ≤30bp should enrich for events directly regulated by SRRM3 and
# improve the signal-to-noise ratio of the correlation analysis.
#
# Comparisons analyzed:
#   - Neg_vs_Parental (partial SRRM3 knockdown)
#   - Pos_vs_Parental (SRRM3 overexpression)
#   - KO_vs_Parental  (SRRM3 knockout)
#
# Input:
#   - data/srrm3_and_srrm4/AS-change.csv (dPSI values)
#   - data/srrm3_and_srrm4/event-information.csv (event coordinates in mm9)
#   - results/05_splicing/{comparison}/SE.MATS.JC.txt (rMATS output in mm39)
#   - data/mm9ToMm39.over.chain (liftOver chain file)
#
# Output:
#   - results/12_paper_comparison_double_kd_microexons/{comparison}/
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

# Input files
as_change_file <- file.path(base_dir, "data/srrm3_and_srrm4/AS-change.csv")
event_info_file <- file.path(base_dir, "data/srrm3_and_srrm4/event-information.csv")
chain_file <- file.path(base_dir, "data/mm9ToMm39.over.chain")

# Base output directory
base_output_dir <- file.path(base_dir, "results/12_paper_comparison_double_kd_microexons")
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# Microexon size threshold
max_exon_length <- 30  # bp — SRRM3-dependent microexons

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
size_tolerance <- 5  # bp tolerance for exon size matching

# =============================================================================
# Helper Functions
# =============================================================================

#' Load double knockdown paper data (AS-change + event-information)
#' @return Data frame with merged event info and dPSI values
load_paper_data_double_kd <- function() {
  cat("Loading double knockdown paper data...\n")

  # Load dPSI values
  as_change <- read.csv(as_change_file, stringsAsFactors = FALSE)
  cat(sprintf("  - Loaded %d events from AS-change.csv\n", nrow(as_change)))
  cat(sprintf("  - Columns: %s\n", paste(colnames(as_change), collapse = ", ")))

  # Load event coordinates (mm9, broad regions)
  event_info <- read.csv(event_info_file, stringsAsFactors = FALSE)
  cat(sprintf("  - Loaded %d events from event-information.csv\n", nrow(event_info)))

  # Merge by VastDB.ID and gene
  paper_data <- merge(event_info, as_change, by = c("VastDB.ID", "gene"))
  cat(sprintf("  - Merged data: %d events\n", nrow(paper_data)))

  # Rename and extract dPSI columns
  # Convert from percentage (paper scale) to 0-1 scale (rMATS scale)
  paper_data <- paper_data %>%
    mutate(
      srrm3_srrm4_dPSI = Srrm3.Srrm4.dPSI / 100,
      srrm3_dPSI = Srrm3.dPSI / 100,
      srrm4_dPSI = Srrm4.dPSI / 100,
      srrm3_srrm4_signif = Srrm3.Srrm4.signif,
      srrm3_signif = Srrm3.signif,
      srrm4_signif = Srrm4.signif
    )

  # Filter to events with at least one non-NA dPSI value
  n_before <- nrow(paper_data)
  paper_data <- paper_data %>%
    filter(!is.na(srrm3_srrm4_dPSI) | !is.na(srrm3_dPSI) | !is.na(srrm4_dPSI))
  n_after <- nrow(paper_data)

  cat(sprintf("  - Events with at least one dPSI value: %d / %d\n", n_after, n_before))

  # *** MICROEXON FILTER ***
  n_pre_filter <- nrow(paper_data)
  paper_data <- paper_data %>%
    filter(length <= max_exon_length)
  cat(sprintf("  - Microexon filter (length <= %dbp): %d / %d events retained\n",
              max_exon_length, nrow(paper_data), n_pre_filter))

  cat(sprintf("  - Srrm3+Srrm4 KD: %d microexons with dPSI\n", sum(!is.na(paper_data$srrm3_srrm4_dPSI))))
  cat(sprintf("  - Srrm3 KD: %d microexons with dPSI\n", sum(!is.na(paper_data$srrm3_dPSI))))
  cat(sprintf("  - Srrm4 KD: %d microexons with dPSI\n", sum(!is.na(paper_data$srrm4_dPSI))))
  cat(sprintf("  - Exon size range: %d-%dbp\n", min(paper_data$length), max(paper_data$length)))

  paper_data
}

#' LiftOver broad event regions from mm9 to mm39
#' The coordinates in the double KD data are NOT exon boundaries — they define
#' a broad region containing the cassette exon (typically 3-400x the exon length).
#' We use these only for disambiguation, not for primary matching.
#' @param paper_data Data frame with mm9 coordinates
#' @return Data frame with VastDB.ID, chrom_mm39, broad_start_mm39, broad_end_mm39
liftover_broad_regions <- function(paper_data) {
  cat("LiftOver mm9 -> mm39 (broad event regions)...\n")

  # Load chain file
  chain <- import.chain(chain_file)
  cat("  - Loaded chain file\n")

  # Create GRanges object from paper data (mm9)
  paper_gr_mm9 <- GRanges(
    seqnames = paper_data$chrom,
    ranges = IRanges(start = paper_data$start, end = paper_data$end),
    strand = paper_data$strand,
    VastDB.ID = paper_data$VastDB.ID
  )

  # Perform liftOver
  paper_gr_mm39_list <- liftOver(paper_gr_mm9, chain)

  # Count successful conversions
  n_converted <- sum(lengths(paper_gr_mm39_list) > 0)
  cat(sprintf("  - Successfully converted: %d / %d events (%.1f%%)\n",
              n_converted, length(paper_gr_mm9), 100 * n_converted / length(paper_gr_mm9)))

  # Extract successfully converted coordinates
  successful_idx <- which(lengths(paper_gr_mm39_list) > 0)
  paper_gr_mm39 <- unlist(paper_gr_mm39_list[successful_idx])

  # Track failed conversions
  n_failed <- length(paper_gr_mm9) - n_converted
  if (n_failed > 0) {
    cat(sprintf("  - Failed liftOver: %d events\n", n_failed))
  }

  # Create data frame with mm39 broad region coordinates
  lifted_df <- data.frame(
    VastDB.ID = paper_gr_mm39$VastDB.ID,
    chrom_mm39 = as.character(seqnames(paper_gr_mm39)),
    broad_start_mm39 = start(paper_gr_mm39),
    broad_end_mm39 = end(paper_gr_mm39),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  - Lifted regions: %d\n", nrow(lifted_df)))

  lifted_df
}

#' Load rMATS data for a specific comparison
#' @param comparison Comparison name (e.g., "KO_vs_Parental")
#' @return Data frame with rMATS results
load_rmats_data <- function(comparison) {
  rmats_file <- file.path(base_dir, "results/05_splicing", comparison, "SE.MATS.JC.txt")

  if (!file.exists(rmats_file)) {
    stop(sprintf("rMATS file not found: %s", rmats_file))
  }

  rmats <- read.delim(rmats_file, stringsAsFactors = FALSE, check.names = TRUE)

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

#' Build gene-indexed lookup for rMATS events (for fast matching)
#' @param rmats rMATS data frame
#' @return Named list: gene_upper -> subset of rmats rows
build_rmats_gene_index <- function(rmats) {
  rmats$gene_upper <- toupper(rmats$geneSymbol)
  split(rmats, rmats$gene_upper)
}

#' Find match by gene name + exon size, with optional containment disambiguation
#' Uses base R indexing (no dplyr) for performance in hot loop.
#' @param rmats_by_gene Pre-indexed rMATS lookup (from build_rmats_gene_index)
#' @param gene_name Gene name from paper
#' @param paper_length Exon length from paper
#' @param lifted_row Row from liftover results (NULL if liftover failed for this event)
#' @param size_tol Size tolerance in bp
#' @return List with best row index, n_all candidates, and used_containment flag; or NULL
find_gene_size_match <- function(rmats_by_gene, gene_name, paper_length, lifted_row = NULL,
                                 size_tol = 5) {
  gene_key <- toupper(gene_name)

  # Step 1: Lookup by gene name
  candidates <- rmats_by_gene[[gene_key]]
  if (is.null(candidates) || nrow(candidates) == 0) return(NULL)

  # Step 2: Filter by exon size tolerance (base R)
  size_mask <- abs(candidates$exon_size - paper_length) <= size_tol
  if (!any(size_mask)) return(NULL)

  matches <- candidates[size_mask, , drop = FALSE]
  n_all <- nrow(matches)
  used_containment <- FALSE

  # Step 3: If multiple candidates and liftOver available, check containment
  if (nrow(matches) > 1 && !is.null(lifted_row)) {
    cont_mask <- matches$chr == lifted_row$chrom_mm39 &
      matches$exon_start_1based >= lifted_row$broad_start_mm39 &
      matches$exonEnd <= lifted_row$broad_end_mm39
    if (any(cont_mask) && sum(cont_mask) < nrow(matches)) {
      matches <- matches[cont_mask, , drop = FALSE]
      used_containment <- TRUE
    }
  }

  # Step 4: Tiebreaker — prefer FDR < 0.05, then strongest |dPSI|
  if (nrow(matches) > 1) {
    sig_mask <- matches$FDR < 0.05
    if (any(sig_mask, na.rm = TRUE)) {
      matches <- matches[which(sig_mask), , drop = FALSE]
    }
    best_idx <- which.max(abs(matches$IncLevelDifference))
    matches <- matches[best_idx, , drop = FALSE]
  }

  return(list(best = matches, n_all = n_all, used_containment = used_containment))
}

#' Match paper events to rMATS events using gene name + exon size
#' Memory-efficient: uses pre-allocated vectors instead of list-of-lists.
#' @param paper_data Paper data frame with event info and dPSI values
#' @param lifted_df LiftOver results (VastDB.ID, chrom_mm39, broad_start_mm39, broad_end_mm39)
#' @param rmats rMATS data frame
#' @param size_tol Size tolerance in bp
#' @return List with 'comparison_df' and 'multi_match_records'
match_events_double_kd <- function(paper_data, lifted_df, rmats, size_tol = 5) {
  # Pre-index rMATS by gene for fast lookup
  rmats_by_gene <- build_rmats_gene_index(rmats)
  cat(sprintf("  - Built gene index: %d unique genes in rMATS\n", length(rmats_by_gene)))

  # Create lifted lookup by VastDB.ID
  lifted_lookup <- split(lifted_df, lifted_df$VastDB.ID)

  n <- nrow(paper_data)

  # Pre-allocate result vectors
  res_match_type <- character(n)
  res_n_candidates <- integer(n)
  res_rmats_id <- integer(n)
  res_rmats_gene <- character(n)
  res_rmats_chr <- character(n)
  res_rmats_start <- integer(n)
  res_rmats_end <- integer(n)
  res_rmats_exon_size <- integer(n)
  res_rmats_dPSI <- numeric(n)
  res_rmats_FDR <- numeric(n)
  res_rmats_PValue <- numeric(n)

  # Initialize with NA
  res_match_type[] <- "unmatched"
  res_n_candidates[] <- 0L
  res_rmats_id[] <- NA_integer_
  res_rmats_gene[] <- NA_character_
  res_rmats_chr[] <- NA_character_
  res_rmats_start[] <- NA_integer_
  res_rmats_end[] <- NA_integer_
  res_rmats_exon_size[] <- NA_integer_
  res_rmats_dPSI[] <- NA_real_
  res_rmats_FDR[] <- NA_real_
  res_rmats_PValue[] <- NA_real_

  # Collect multi-match events (capped for memory)
  max_multi_events <- 500
  multi_match_list <- list()
  n_multi_recorded <- 0

  n_gene_size <- 0
  n_gene_size_contained <- 0
  n_unmatched <- 0
  n_multi_total <- 0

  for (i in seq_len(n)) {
    vastdb_id <- paper_data$VastDB.ID[i]
    gene_name <- paper_data$gene[i]
    paper_length <- paper_data$length[i]

    # Get lifted coordinates if available
    lifted_row <- lifted_lookup[[vastdb_id]]
    if (!is.null(lifted_row) && nrow(lifted_row) > 0) {
      lifted_row <- lifted_row[1, ]
    } else {
      lifted_row <- NULL
    }

    # Find match
    result <- find_gene_size_match(
      rmats_by_gene = rmats_by_gene,
      gene_name = gene_name,
      paper_length = paper_length,
      lifted_row = lifted_row,
      size_tol = size_tol
    )

    if (!is.null(result)) {
      best <- result$best
      match_type <- if (result$used_containment) "gene_size_contained" else "gene_size"

      if (result$used_containment) n_gene_size_contained <- n_gene_size_contained + 1
      else n_gene_size <- n_gene_size + 1

      res_match_type[i] <- match_type
      res_n_candidates[i] <- result$n_all
      res_rmats_id[i] <- best$rmats_ID[1]
      res_rmats_gene[i] <- best$geneSymbol[1]
      res_rmats_chr[i] <- best$chr[1]
      res_rmats_start[i] <- best$exon_start_1based[1]
      res_rmats_end[i] <- best$exonEnd[1]
      res_rmats_exon_size[i] <- best$exon_size[1]
      res_rmats_dPSI[i] <- best$IncLevelDifference[1]
      res_rmats_FDR[i] <- best$FDR[1]
      res_rmats_PValue[i] <- best$PValue[1]

      # Record multi-match events (capped)
      if (result$n_all > 1) {
        n_multi_total <- n_multi_total + 1
        if (n_multi_recorded < max_multi_events) {
          # Re-fetch candidates to record all of them
          gene_key <- toupper(gene_name)
          candidates <- rmats_by_gene[[gene_key]]
          all_matches <- candidates[abs(candidates$exon_size - paper_length) <= size_tol, , drop = FALSE]
          for (j in seq_len(nrow(all_matches))) {
            row <- all_matches[j, ]
            multi_match_list[[length(multi_match_list) + 1]] <- data.frame(
              VastDB.ID = vastdb_id,
              paper_gene = gene_name,
              paper_length = paper_length,
              match_type = match_type,
              n_candidates = result$n_all,
              selected = (row$rmats_ID == best$rmats_ID[1]),
              rmats_ID = row$rmats_ID,
              geneSymbol = row$geneSymbol,
              chr = row$chr,
              exonStart_0base = row$exonStart_0base,
              exonEnd = row$exonEnd,
              exon_size = row$exon_size,
              strand = row$strand,
              IncLevelDifference = row$IncLevelDifference,
              FDR = row$FDR,
              stringsAsFactors = FALSE
            )
          }
          n_multi_recorded <- n_multi_recorded + 1
        }
      }
    } else {
      n_unmatched <- n_unmatched + 1
    }

    # Progress reporting (every 5000 events)
    if (i %% 5000 == 0) {
      cat(sprintf("  - Processed %d / %d events...\n", i, n))
    }
  }

  cat(sprintf("  - Gene+size matches: %d\n", n_gene_size))
  cat(sprintf("  - Gene+size+contained matches: %d\n", n_gene_size_contained))
  cat(sprintf("  - Unmatched: %d\n", n_unmatched))
  if (n_multi_total > max_multi_events) {
    cat(sprintf("  - Multi-match events: %d total (%d recorded in detail)\n",
                n_multi_total, n_multi_recorded))
  } else if (n_multi_total > 0) {
    cat(sprintf("  - Multi-match events: %d\n", n_multi_total))
  }

  # Build match result data frame from pre-allocated vectors
  match_df <- data.frame(
    VastDB.ID = paper_data$VastDB.ID,
    match_type = res_match_type,
    n_rmats_candidates = res_n_candidates,
    rmats_id = res_rmats_id,
    rmats_gene = res_rmats_gene,
    rmats_chr = res_rmats_chr,
    rmats_start = res_rmats_start,
    rmats_end = res_rmats_end,
    rmats_exon_size = res_rmats_exon_size,
    rmats_dPSI = res_rmats_dPSI,
    rmats_FDR = res_rmats_FDR,
    rmats_PValue = res_rmats_PValue,
    stringsAsFactors = FALSE
  )

  # Merge with paper data
  comparison_df <- merge(paper_data, match_df, by = "VastDB.ID", all.x = TRUE)

  # Combine multi-match records
  if (length(multi_match_list) > 0) {
    multi_match_df <- bind_rows(multi_match_list)
  } else {
    multi_match_df <- data.frame()
  }

  list(comparison_df = comparison_df, multi_match_records = multi_match_df)
}

#' Calculate correlation statistics for 3 paper dPSI columns vs rMATS dPSI
#' @param comparison_df Comparison data frame with matched events
#' @param condition_short Short condition name (e.g., "KO")
#' @return List with statistics for each paper KD condition
calculate_statistics <- function(comparison_df, condition_short) {
  matched_for_stats <- comparison_df %>%
    filter(match_type != "unmatched", !is.na(rmats_dPSI))

  results <- list()

  # Define the 3 paper KD conditions
  paper_conditions <- list(
    srrm3_srrm4 = list(col = "srrm3_srrm4_dPSI", label = "Srrm3+Srrm4 DKD"),
    srrm3 = list(col = "srrm3_dPSI", label = "Srrm3 KD"),
    srrm4 = list(col = "srrm4_dPSI", label = "Srrm4 KD")
  )

  summary_rows <- list()

  for (cond_name in names(paper_conditions)) {
    cond <- paper_conditions[[cond_name]]
    col_name <- cond$col

    cond_stats <- matched_for_stats %>%
      filter(!is.na(.data[[col_name]]))

    if (nrow(cond_stats) >= 3) {
      pearson_test <- cor.test(cond_stats[[col_name]], cond_stats$rmats_dPSI,
                               method = "pearson")
      spearman_test <- cor.test(cond_stats[[col_name]], cond_stats$rmats_dPSI,
                                method = "spearman")
      concordance <- mean(sign(cond_stats[[col_name]]) == sign(cond_stats$rmats_dPSI),
                          na.rm = TRUE) * 100

      results[[cond_name]] <- list(
        data = cond_stats,
        n = nrow(cond_stats),
        pearson = pearson_test,
        spearman = spearman_test,
        concordance = concordance,
        label = cond$label
      )

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        comparison = sprintf("%s_vs_SRRM3_%s", gsub(" ", "_", cond$label), condition_short),
        paper_condition = cond$label,
        n_events = nrow(cond_stats),
        pearson_r = as.numeric(pearson_test$estimate),
        pearson_p = pearson_test$p.value,
        spearman_rho = as.numeric(spearman_test$estimate),
        spearman_p = spearman_test$p.value,
        direction_concordance_pct = concordance,
        stringsAsFactors = FALSE
      )
    } else {
      results[[cond_name]] <- list(
        data = cond_stats,
        n = nrow(cond_stats),
        pearson = list(estimate = NA, p.value = NA),
        spearman = list(estimate = NA, p.value = NA),
        concordance = NA,
        label = cond$label
      )

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        comparison = sprintf("%s_vs_SRRM3_%s", gsub(" ", "_", cond$label), condition_short),
        paper_condition = cond$label,
        n_events = nrow(cond_stats),
        pearson_r = NA,
        pearson_p = NA,
        spearman_rho = NA,
        spearman_p = NA,
        direction_concordance_pct = NA,
        stringsAsFactors = FALSE
      )
    }
  }

  results$summary <- bind_rows(summary_rows)
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
  plot_colors <- list(
    srrm3_srrm4 = list(line = "#7570B3", fill = "#BCBDDC", file_tag = "srrm3_srrm4"),
    srrm3 = list(line = "steelblue", fill = "lightblue", file_tag = "srrm3"),
    srrm4 = list(line = "#B2182B", fill = "#FDDBC7", file_tag = "srrm4")
  )

  for (cond_name in c("srrm3_srrm4", "srrm3", "srrm4")) {
    cond_stats <- stats[[cond_name]]
    colors <- plot_colors[[cond_name]]
    col_name <- paste0(cond_name, "_dPSI")

    if (cond_stats$n >= 3 && !is.na(cond_stats$pearson$estimate)) {
      plot_data <- cond_stats$data

      # For large datasets, only label the top genes by |dPSI|
      n_label <- min(30, nrow(plot_data))
      plot_data <- plot_data %>%
        mutate(
          label_gene = ifelse(
            rank(-abs(rmats_dPSI)) <= n_label | rank(-abs(.data[[col_name]])) <= n_label,
            gene, ""
          )
        )

      p <- ggplot(plot_data, aes(x = .data[[col_name]], y = rmats_dPSI)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
        geom_smooth(method = "lm", se = TRUE, color = colors$line, fill = colors$fill, alpha = 0.3) +
        geom_point(aes(color = match_type), size = 2, alpha = 0.6) +
        scale_color_manual(values = c("gene_size" = "#2166AC",
                                       "gene_size_contained" = "#F4A582"),
                           labels = c("gene_size" = "Gene + Exon Size",
                                      "gene_size_contained" = "Gene + Size + Region")) +
        labs(
          title = sprintf("%s (Paper) vs %s (This Study) — Microexons", cond_stats$label, condition_desc),
          subtitle = sprintf("N = %d | Pearson r = %.3f (p = %.2e) | Concordance = %.1f%%",
                            cond_stats$n,
                            as.numeric(cond_stats$pearson$estimate),
                            cond_stats$pearson$p.value,
                            cond_stats$concordance),
          x = sprintf("Paper %s dPSI", cond_stats$label),
          y = sprintf("This Study %s dPSI", condition_short),
          color = "Match Type"
        ) +
        theme_publication

      # Add gene labels only if not too many points
      if (nrow(plot_data) <= 200) {
        p <- p + geom_text(aes(label = label_gene), hjust = -0.1, vjust = -0.5, size = 2,
                           check_overlap = TRUE)
      }

      ggsave(file.path(output_dir, sprintf("scatter_%s_vs_%s.pdf", colors$file_tag, condition_short)),
             p, width = 8, height = 7)
      cat(sprintf("  - Saved: scatter_%s_vs_%s.pdf\n", colors$file_tag, condition_short))

      plots[[cond_name]] <- p
    } else {
      cat(sprintf("  - Skipped %s scatter plot (insufficient data: n=%d)\n",
                  cond_stats$label, cond_stats$n))
      plots[[cond_name]] <- NULL
    }
  }

  # Combined 3-panel scatter plot
  valid_plots <- plots[!sapply(plots, is.null)]
  if (length(valid_plots) >= 2) {
    pdf(file.path(output_dir, "scatter_combined.pdf"),
        width = 7 * length(valid_plots), height = 7)
    do.call(grid.arrange, c(valid_plots, ncol = length(valid_plots),
                            top = sprintf("Double KD Microexon Comparison - %s", comparison)))
    dev.off()
    cat("  - Saved: scatter_combined.pdf\n")
  }

  # Heatmap of dPSI values (top matched events)
  # Require ALL 3 paper dPSI columns to be non-NA (no grey/blank cells)
  # and require meaningful signal from at least one source
  matched_for_heatmap <- comparison_df %>%
    filter(match_type != "unmatched", !is.na(rmats_dPSI)) %>%
    filter(!is.na(srrm3_srrm4_dPSI) & !is.na(srrm3_dPSI) & !is.na(srrm4_dPSI)) %>%
    mutate(
      max_paper_dPSI = pmax(abs(srrm3_srrm4_dPSI), abs(srrm3_dPSI), abs(srrm4_dPSI)),
      combined_signal = max_paper_dPSI + abs(rmats_dPSI)
    ) %>%
    filter(max_paper_dPSI > 0.01 | abs(rmats_dPSI) > 0.05)

  # Select top events ranked by combined signal from both datasets
  max_heatmap_genes <- 50
  if (nrow(matched_for_heatmap) > max_heatmap_genes) {
    matched_for_heatmap <- matched_for_heatmap %>%
      arrange(desc(combined_signal)) %>%
      slice(1:max_heatmap_genes)
  }

  if (nrow(matched_for_heatmap) > 0) {
    heatmap_data <- matched_for_heatmap %>%
      select(gene, VastDB.ID, srrm3_srrm4_dPSI, srrm3_dPSI, srrm4_dPSI, rmats_dPSI) %>%
      mutate(gene_label = paste0(gene, " (", VastDB.ID, ")")) %>%
      select(-VastDB.ID) %>%
      rename(
        `Paper Srrm3+Srrm4 DKD` = srrm3_srrm4_dPSI,
        `Paper Srrm3 KD` = srrm3_dPSI,
        `Paper Srrm4 KD` = srrm4_dPSI,
        !!sprintf("This Study %s", condition_short) := rmats_dPSI
      ) %>%
      pivot_longer(cols = c("Paper Srrm3+Srrm4 DKD", "Paper Srrm3 KD", "Paper Srrm4 KD",
                            sprintf("This Study %s", condition_short)),
                   names_to = "Condition", values_to = "dPSI") %>%
      mutate(Condition = factor(Condition, levels = c("Paper Srrm3+Srrm4 DKD",
                                                       "Paper Srrm3 KD", "Paper Srrm4 KD",
                                                       sprintf("This Study %s", condition_short))))

    p_heatmap <- ggplot(heatmap_data, aes(x = Condition, y = gene_label, fill = dPSI)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-1, 1),
        name = "dPSI"
      ) +
      labs(
        title = sprintf("Microexon dPSI - Top %d Events - %s", nrow(matched_for_heatmap), comparison),
        subtitle = sprintf("Microexons (≤%dbp) — ranked by |rMATS dPSI|", max_exon_length),
        x = "",
        y = "Gene (VastDB ID)"
      ) +
      theme_publication +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7))

    # Add numeric dPSI values — scale text size to number of genes
    label_size <- if (nrow(matched_for_heatmap) <= 20) 2.5 else if (nrow(matched_for_heatmap) <= 35) 2 else 1.5
    p_heatmap <- p_heatmap +
      geom_text(aes(label = ifelse(!is.na(dPSI), sprintf("%.2f", dPSI), "")),
                size = label_size, color = "black")

    ggsave(file.path(output_dir, "dpsi_heatmap.pdf"), p_heatmap,
           width = 7, height = max(5, nrow(matched_for_heatmap) * 0.25 + 2))
    cat("  - Saved: dpsi_heatmap.pdf\n")
  }

  # Summary bar plot
  match_df <- comparison_df %>%
    select(VastDB.ID, match_type) %>%
    distinct()

  match_summary <- match_df %>%
    count(match_type) %>%
    mutate(
      match_type = factor(match_type,
                          levels = c("gene_size", "gene_size_contained", "unmatched"))
    )

  p_summary <- ggplot(match_summary, aes(x = match_type, y = n, fill = match_type)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("gene_size" = "#2166AC",
                                  "gene_size_contained" = "#F4A582",
                                  "unmatched" = "#D1D1D1")) +
    labs(
      title = sprintf("Microexon Matching Summary - %s", comparison),
      subtitle = sprintf("Microexons (≤%dbp): %d events | Gene+size matching (±%dbp tolerance)",
                         max_exon_length, nrow(match_df), size_tolerance),
      x = "Match Type",
      y = "Number of Events"
    ) +
    theme_publication +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("Gene+Size\nMatch", "Gene+Size\n+Contained", "Unmatched"))

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
  cat(sprintf("  - Saved: matched_events.csv (%d events)\n", nrow(matched_events)))

  # Unmatched events
  unmatched_events <- comparison_df %>%
    filter(match_type == "unmatched")
  write.csv(unmatched_events, file.path(output_dir, "unmatched_events.csv"), row.names = FALSE)
  cat(sprintf("  - Saved: unmatched_events.csv (%d events)\n", nrow(unmatched_events)))

  # Statistics summary
  write.csv(stats$summary, file.path(output_dir, "statistics_summary.csv"), row.names = FALSE)
  cat("  - Saved: statistics_summary.csv\n")

  # Multi-match report
  report_file <- file.path(output_dir, "multi_match_events.txt")
  if (nrow(multi_match_df) > 0) {
    n_events <- length(unique(multi_match_df$VastDB.ID))
    con <- file(report_file, "w")
    writeLines(sprintf("Multi-Match Event Report: %s (Double KD - Microexons)", comparison), con)
    writeLines(sprintf("Generated: %s", Sys.time()), con)
    writeLines(paste0(rep("=", 80), collapse = ""), con)
    writeLines("", con)
    writeLines(sprintf("Paper events with multiple rMATS candidates: %d", n_events), con)
    writeLines(sprintf("Total rMATS candidate entries: %d", nrow(multi_match_df)), con)
    writeLines("Selection criterion: prefer significant (FDR < 0.05), then strongest |dPSI|", con)
    writeLines("", con)
    writeLines(paste0(rep("-", 80), collapse = ""), con)

    for (evt in unique(multi_match_df$VastDB.ID)) {
      evt_rows <- multi_match_df %>% filter(VastDB.ID == evt)
      gene <- evt_rows$paper_gene[1]
      n_cand <- evt_rows$n_candidates[1]
      mtype <- evt_rows$match_type[1]
      paper_len <- evt_rows$paper_length[1]

      writeLines("", con)
      writeLines(sprintf("VastDB: %s  |  Gene: %s  |  Paper length: %dbp  |  Match: %s  |  Candidates: %d",
                         evt, gene, paper_len, mtype, n_cand), con)

      for (k in 1:nrow(evt_rows)) {
        r <- evt_rows[k, ]
        tag <- if (r$selected) {
          if (r$FDR < 0.05) " << SELECTED (significant, strongest |dPSI|)"
          else " << SELECTED (none significant, strongest |dPSI|)"
        } else ""
        writeLines(sprintf(
          "  rMATS ID=%d  %s:%d-%d (%s)  exon=%dbp  dPSI=%.4f  FDR=%.2e%s",
          r$rmats_ID, r$chr, r$exonStart_0base, r$exonEnd, r$strand,
          r$exon_size, r$IncLevelDifference, r$FDR, tag
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
    writeLines(sprintf("Multi-Match Event Report: %s (Double KD - Microexons)", comparison), con)
    writeLines(sprintf("Generated: %s", Sys.time()), con)
    writeLines(paste0(rep("=", 80), collapse = ""), con)
    writeLines("", con)
    writeLines("No multi-match events found.", con)
    close(con)
    cat("  - Saved: multi_match_events.txt (no multi-matches)\n")
  }

  # Also save multi-match data as CSV
  if (nrow(multi_match_df) > 0) {
    write.csv(multi_match_df, file.path(output_dir, "multi_match_events.csv"), row.names = FALSE)
    cat("  - Saved: multi_match_events.csv\n")
  }
}

#' Process a single comparison
#' @param comparison Comparison name
#' @param paper_data Paper data with dPSI values
#' @param lifted_df LiftOver results
#' @param condition_label Condition label list
#' @return List with results
process_comparison <- function(comparison, paper_data, lifted_df, condition_label) {
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
  cat("\nMatching events (gene name + exon size)...\n")
  match_result <- match_events_double_kd(paper_data, lifted_df, rmats, size_tol = size_tolerance)
  comparison_df <- match_result$comparison_df
  multi_match_df <- match_result$multi_match_records

  # Print matching summary
  n_gene_size <- sum(comparison_df$match_type == "gene_size")
  n_contained <- sum(comparison_df$match_type == "gene_size_contained")
  n_unmatched <- sum(comparison_df$match_type == "unmatched")
  n_total <- nrow(comparison_df)

  cat(sprintf("\n  Matching summary:\n"))
  cat(sprintf("  - Gene+size matches: %d\n", n_gene_size))
  cat(sprintf("  - Gene+size+contained matches: %d\n", n_contained))
  cat(sprintf("  - Unmatched: %d\n", n_unmatched))
  cat(sprintf("  - Total matched: %d / %d (%.1f%%)\n",
              n_gene_size + n_contained, n_total,
              100 * (n_gene_size + n_contained) / n_total))

  if (nrow(multi_match_df) > 0) {
    n_multi <- length(unique(multi_match_df$VastDB.ID))
    cat(sprintf("  - Events with multiple rMATS candidates: %d\n", n_multi))
  }

  # Calculate statistics
  cat("\nCalculating statistics...\n")
  stats <- calculate_statistics(comparison_df, condition_label$short)

  # Print correlation summary
  for (cond_name in c("srrm3_srrm4", "srrm3", "srrm4")) {
    cond <- stats[[cond_name]]
    cat(sprintf("\n  %s vs %s:\n", cond$label, condition_label$short))
    if (!is.na(cond$pearson$estimate)) {
      cat(sprintf("    - N events: %d\n", cond$n))
      cat(sprintf("    - Pearson r: %.3f (p = %.2e)\n",
                  as.numeric(cond$pearson$estimate), cond$pearson$p.value))
      cat(sprintf("    - Spearman rho: %.3f (p = %.2e)\n",
                  as.numeric(cond$spearman$estimate), cond$spearman$p.value))
      cat(sprintf("    - Direction concordance: %.1f%%\n", cond$concordance))
    } else {
      cat(sprintf("    - Insufficient data for correlation (n=%d)\n", cond$n))
    }
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
    n_matched = n_gene_size + n_contained,
    n_total = n_total
  )
}

# =============================================================================
# Main Analysis
# =============================================================================

cat("=============================================================================\n")
cat("SRRM3 Paper Validation - DOUBLE KNOCKDOWN - MICROEXONS ONLY\n")
cat("Gonatopoulos-Pournatzis et al. 2018 - Srrm3+Srrm4 Data (≤30bp)\n")
cat("=============================================================================\n\n")

# Step 1: Load paper data (done once)
cat("Step 1: Loading and preparing double knockdown paper data...\n")
cat("-----------------------------------------------------------------------------\n")
paper_data <- load_paper_data_double_kd()

# Step 2: LiftOver broad event regions (done once)
cat("\nStep 2: Converting broad event regions (mm9 -> mm39)...\n")
cat("-----------------------------------------------------------------------------\n")
lifted_df <- liftover_broad_regions(paper_data)

# Step 3: Process each comparison
cat("\nStep 3: Processing comparisons...\n")
cat("-----------------------------------------------------------------------------\n")

all_results <- list()
for (comparison in comparisons) {
  condition_label <- condition_labels[[comparison]]
  all_results[[comparison]] <- process_comparison(comparison, paper_data, lifted_df, condition_label)
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
for (cond_label in c("Srrm3+Srrm4 DKD", "Srrm3 KD", "Srrm4 KD")) {
  cat(sprintf("%s (Paper) vs This Study:\n", cond_label))
  cat("-----------------------------------------------------------------------------\n")
  cond_summary <- all_stats %>%
    filter(paper_condition == cond_label) %>%
    select(study_comparison, n_events, pearson_r, pearson_p, spearman_rho, direction_concordance_pct)

  for (i in 1:nrow(cond_summary)) {
    row <- cond_summary[i, ]
    if (!is.na(row$pearson_r)) {
      cat(sprintf("  %s:\n", row$study_comparison))
      cat(sprintf("    N = %d, r = %.3f (p = %.2e), rho = %.3f, concordance = %.1f%%\n",
                  row$n_events, row$pearson_r, row$pearson_p,
                  row$spearman_rho, row$direction_concordance_pct))
    } else {
      cat(sprintf("  %s: insufficient data (N = %d)\n", row$study_comparison, row$n_events))
    }
  }
  cat("\n")
}

# Save cross-comparison summary
write.csv(all_stats, file.path(base_output_dir, "cross_comparison_statistics.csv"), row.names = FALSE)
cat(sprintf("Saved cross-comparison summary: %s/cross_comparison_statistics.csv\n",
            base_output_dir))
cat(sprintf("  -> %d rows (3 comparisons x 3 paper KD conditions)\n", nrow(all_stats)))

# =============================================================================
# Final Summary
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")

cat("Output directories created:\n")
for (comparison in comparisons) {
  result <- all_results[[comparison]]
  output_dir <- file.path(base_output_dir, comparison)
  cat(sprintf("  - %s/  (matched: %d / %d events)\n",
              output_dir, result$n_matched, result$n_total))
  cat("      scatter_srrm3_srrm4_vs_*.pdf, scatter_srrm3_vs_*.pdf, scatter_srrm4_vs_*.pdf\n")
  cat("      scatter_combined.pdf, dpsi_heatmap.pdf, summary_bar.pdf\n")
  cat("      paper_comparison_full.csv, matched_events.csv\n")
  cat("      unmatched_events.csv, statistics_summary.csv\n")
  cat("      multi_match_events.txt, multi_match_events.csv\n")
}

cat(sprintf("\nCross-comparison summary: %s/cross_comparison_statistics.csv\n", base_output_dir))

cat("\n=============================================================================\n")
