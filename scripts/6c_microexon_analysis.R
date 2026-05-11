#!/usr/bin/env Rscript
# =============================================================================
# Microexon Analysis
# Project: 90-1239779069
# Analysis: Exon size stratification (0-30bp, 30-50bp, >50bp)
# Focus on skipped exons (SE) events
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(VennDiagram)
    library(grid)
    library(RColorBrewer)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/07_microexon_analysis")

# Thresholds (matching paper standards)
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1

# Exon size bins
MICROEXON_MAX <- 30      # 0-30bp: microexons
SMALL_EXON_MAX <- 50     # 30-50bp: small exons
                          # >50bp: regular exons

cat("============================================\n")
cat("Microexon Analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Size bins: Microexon (0-30bp), Small (30-50bp), Regular (>50bp)\n")
cat("Thresholds: FDR <", FDR_THRESHOLD, ", |dPSI| >", DPSI_THRESHOLD, "\n")
cat("============================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "size_distribution"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "microexon_specific"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "common_events"), showWarnings = FALSE)

# Define comparisons
parental_comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")
all_comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental",
                     "Pos_vs_Neg", "KO_vs_Neg", "KO_vs_Pos")

# Color palettes
size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")
comparison_colors <- c("Neg_vs_Parental" = "#377EB8",
                       "Pos_vs_Parental" = "#4DAF4A",
                       "KO_vs_Parental" = "#984EA3")

# =============================================================================
# Function: Load SE events and calculate exon size
# =============================================================================
load_se_with_size <- function(comparison, splicing_dir) {
    file_path <- file.path(splicing_dir, comparison, "SE.MATS.JC.txt")

    if (!file.exists(file_path)) {
        warning(paste("File not found:", file_path))
        return(NULL)
    }

    df <- read.delim(file_path, stringsAsFactors = FALSE)

    # Calculate exon length (exonEnd - exonStart_0base)
    df$exon_length <- df$exonEnd - df$exonStart_0base

    # Classify by exon size
    df$size_category <- case_when(
        df$exon_length <= MICROEXON_MAX ~ "Microexon (0-30bp)",
        df$exon_length <= SMALL_EXON_MAX ~ "Small (30-50bp)",
        TRUE ~ "Regular (>50bp)"
    )
    df$size_category <- factor(df$size_category,
                               levels = c("Microexon (0-30bp)", "Small (30-50bp)", "Regular (>50bp)"))

    # Mark significant events
    df$significant <- df$FDR < FDR_THRESHOLD & abs(df$IncLevelDifference) > DPSI_THRESHOLD

    df$comparison <- comparison
    
    # Calculate average coverage (IJC + SJC)
    sum_counts <- function(x) {
        if (is.null(x) || is.na(x)) return(0)
        sum(as.numeric(unlist(strsplit(as.character(x), ","))), na.rm = TRUE)
    }
    
    if(all(c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2") %in% colnames(df))) {
        df$avg_reads <- (
            sapply(df$IJC_SAMPLE_1, sum_counts) + 
            sapply(df$SJC_SAMPLE_1, sum_counts) + 
            sapply(df$IJC_SAMPLE_2, sum_counts) + 
            sapply(df$SJC_SAMPLE_2, sum_counts)
        ) / 6 
        
        # Apply filter: Average reads >= 10
        df <- df[df$avg_reads >= 10, ]
    }

    return(df)
}

# =============================================================================
# Function: Classify exons by size
# =============================================================================
classify_exons_by_size <- function(splicing_dir, comparisons) {
    cat("Classifying exons by size...\n")

    all_data <- list()

    for (comp in comparisons) {
        df <- load_se_with_size(comp, splicing_dir)
        if (!is.null(df)) {
            all_data[[comp]] <- df
        }
    }

    combined <- do.call(rbind, all_data)

    # Summary by size category
    size_summary <- combined %>%
        group_by(comparison, size_category) %>%
        summarize(
            total = n(),
            significant = sum(significant),
            sig_inclusion = sum(significant & IncLevelDifference > 0),
            sig_exclusion = sum(significant & IncLevelDifference < 0),
            mean_exon_length = mean(exon_length),
            .groups = "drop"
        )

    cat("\nExon size classification summary:\n")
    print(size_summary %>% filter(comparison %in% parental_comparisons))

    return(list(data = combined, summary = size_summary))
}

# =============================================================================
# Function: Plot exon size distribution
# =============================================================================
plot_size_distribution <- function(data, output_dir) {
    cat("Generating exon size distribution plots...\n")

    # 1. Histogram of all exon sizes
    p1 <- ggplot(data, aes(x = exon_length)) +
        geom_histogram(bins = 100, fill = "#3498db", alpha = 0.7) +
        geom_vline(xintercept = c(MICROEXON_MAX, SMALL_EXON_MAX),
                   linetype = "dashed", color = "red") +
        scale_x_log10() +
        facet_wrap(~comparison) +
        theme_bw() +
        labs(
            x = "Exon Length (bp, log scale)",
            y = "Count",
            title = "Distribution of Skipped Exon Sizes",
            subtitle = paste0("Vertical lines at ", MICROEXON_MAX, "bp and ", SMALL_EXON_MAX, "bp")
        )

    ggsave(file.path(output_dir, "size_distribution", "exon_size_histogram.pdf"),
           p1, width = 14, height = 10)

    # 2. Stacked bar by size category - all events
    size_counts <- data %>%
        group_by(comparison, size_category) %>%
        summarize(count = n(), .groups = "drop")

    p2 <- ggplot(size_counts, aes(x = comparison, y = count, fill = size_category)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = size_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Number of Events",
            title = "All Skipped Exon Events by Size Category",
            fill = "Exon Size"
        )

    ggsave(file.path(output_dir, "size_distribution", "all_events_by_size.pdf"),
           p2, width = 12, height = 8)

    # 3. Stacked bar by size category - significant only
    sig_counts <- data %>%
        filter(significant) %>%
        group_by(comparison, size_category) %>%
        summarize(count = n(), .groups = "drop")

    p3 <- ggplot(sig_counts, aes(x = comparison, y = count, fill = size_category)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = size_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Number of Significant Events",
            title = paste0("Significant Differential Splicing Events by Exon Size"),
            subtitle = paste0("FDR < ", FDR_THRESHOLD, ", |dPSI| > ", DPSI_THRESHOLD),
            fill = "Exon Size"
        )

    ggsave(file.path(output_dir, "size_distribution", "significant_by_size.pdf"),
           p3, width = 12, height = 8)

    # 4. Percentage breakdown
    pct_data <- data %>%
        filter(significant) %>%
        group_by(comparison, size_category) %>%
        summarize(count = n(), .groups = "drop") %>%
        group_by(comparison) %>%
        mutate(percentage = count / sum(count) * 100)

    p4 <- ggplot(pct_data, aes(x = comparison, y = percentage, fill = size_category)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = size_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Percentage",
            title = "Percentage Distribution of Significant Events by Exon Size",
            fill = "Exon Size"
        )

    ggsave(file.path(output_dir, "size_distribution", "percentage_by_size.pdf"),
           p4, width = 12, height = 8)
}

# =============================================================================
# Function: Analyze microexon splicing specifically
# =============================================================================
analyze_microexon_splicing <- function(data, output_dir) {
    cat("Analyzing microexon splicing patterns...\n")

    # Filter for microexons only
    microexons <- data %>%
        filter(size_category == "Microexon (0-30bp)")

    if (nrow(microexons) == 0) {
        warning("No microexons found")
        return(NULL)
    }

    cat(paste0("  Total microexon events: ", nrow(microexons), "\n"))
    cat(paste0("  Significant microexon events: ", sum(microexons$significant), "\n"))

    # 1. dPSI distribution for microexons vs regular exons
    sig_data <- data %>%
        filter(significant)

    p1 <- ggplot(sig_data, aes(x = IncLevelDifference, fill = size_category)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~comparison) +
        scale_fill_manual(values = size_colors) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        theme_bw() +
        labs(
            x = "dPSI (IncLevelDifference)",
            y = "Density",
            title = "dPSI Distribution: Microexons vs Regular Exons",
            fill = "Exon Size"
        )

    ggsave(file.path(output_dir, "microexon_specific", "dpsi_by_size.pdf"),
           p1, width = 14, height = 10)

    # 2. Directionality for microexons
    direction_data <- sig_data %>%
        mutate(direction = ifelse(IncLevelDifference > 0, "Inclusion", "Exclusion")) %>%
        group_by(comparison, size_category, direction) %>%
        summarize(count = n(), .groups = "drop")

    p2 <- ggplot(direction_data, aes(x = size_category, y = count, fill = direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = c("Inclusion" = "#e74c3c", "Exclusion" = "#3498db")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Exon Size Category",
            y = "Number of Events",
            title = "Splicing Directionality by Exon Size",
            fill = "Direction"
        )

    ggsave(file.path(output_dir, "microexon_specific", "directionality_by_size.pdf"),
           p2, width = 14, height = 10)

    # 3. Box plot of dPSI by size category
    p3 <- ggplot(sig_data, aes(x = size_category, y = IncLevelDifference, fill = size_category)) +
        geom_boxplot(alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = size_colors) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        ) +
        labs(
            x = "Exon Size Category",
            y = "dPSI (IncLevelDifference)",
            title = "dPSI Distribution by Exon Size Category"
        )

    ggsave(file.path(output_dir, "microexon_specific", "dpsi_boxplot_by_size.pdf"),
           p3, width = 14, height = 10)

    # 4. List of significant microexon genes
    sig_microexons <- microexons %>%
        filter(significant) %>%
        select(comparison, GeneID, geneSymbol, chr, exonStart_0base, exonEnd,
               exon_length, IncLevelDifference, FDR) %>%
        arrange(comparison, FDR)

    write.csv(sig_microexons,
              file.path(output_dir, "microexon_specific", "significant_microexons.csv"),
              row.names = FALSE)

    cat(paste0("  Saved ", nrow(sig_microexons), " significant microexon events\n"))

    return(sig_microexons)
}

# =============================================================================
# Function: Compare microexon vs regular exon patterns
# =============================================================================
compare_microexon_vs_regular <- function(data, output_dir) {
    cat("Comparing microexon vs regular exon patterns...\n")

    sig_data <- data %>%
        filter(significant)

    # Statistical comparison
    comparison_stats <- sig_data %>%
        group_by(comparison, size_category) %>%
        summarize(
            n_events = n(),
            mean_dPSI = mean(IncLevelDifference),
            median_dPSI = median(IncLevelDifference),
            sd_dPSI = sd(IncLevelDifference),
            pct_inclusion = sum(IncLevelDifference > 0) / n() * 100,
            pct_exclusion = sum(IncLevelDifference < 0) / n() * 100,
            mean_FDR = mean(FDR),
            .groups = "drop"
        )

    write.csv(comparison_stats,
              file.path(output_dir, "microexon_specific", "size_comparison_stats.csv"),
              row.names = FALSE)

    # Radar/Spider chart style comparison (using bar plot as alternative)
    metrics_long <- comparison_stats %>%
        filter(comparison %in% parental_comparisons) %>%
        select(comparison, size_category, mean_dPSI, pct_inclusion) %>%
        pivot_longer(cols = c(mean_dPSI, pct_inclusion),
                    names_to = "metric", values_to = "value")

    p <- ggplot(metrics_long, aes(x = size_category, y = value, fill = comparison)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~metric, scales = "free_y") +
        scale_fill_manual(values = comparison_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Exon Size Category",
            y = "Value",
            title = "Splicing Metrics Comparison by Exon Size",
            fill = "Comparison"
        )

    ggsave(file.path(output_dir, "microexon_specific", "metrics_comparison.pdf"),
           p, width = 12, height = 8)

    return(comparison_stats)
}

# =============================================================================
# Function: Common microexon events between comparisons
# =============================================================================
find_common_microexon_events <- function(data, parental_comparisons, output_dir) {
    cat("Finding common microexon events between comparisons...\n")

    # Get significant microexon events for each comparison
    sig_microexons <- list()

    for (comp in parental_comparisons) {
        events <- data %>%
            filter(comparison == comp,
                   size_category == "Microexon (0-30bp)",
                   significant) %>%
            mutate(event_id = paste(GeneID, chr, exonStart_0base, exonEnd, sep = "_"))

        sig_microexons[[comp]] <- unique(events$event_id)
        cat(paste0("  ", comp, ": ", length(sig_microexons[[comp]]),
                  " significant microexon events\n"))
    }

    # Create Venn diagram for microexons
    if (length(parental_comparisons) >= 2 && all(sapply(sig_microexons, length) > 0)) {
        venn_list <- sig_microexons
        names(venn_list) <- gsub("_vs_Parental", "", names(venn_list))

        pdf(file.path(output_dir, "common_events", "microexon_venn.pdf"),
            width = 10, height = 10)

        if (length(venn_list) == 3) {
            venn_plot <- draw.triple.venn(
                area1 = length(venn_list[[1]]),
                area2 = length(venn_list[[2]]),
                area3 = length(venn_list[[3]]),
                n12 = length(intersect(venn_list[[1]], venn_list[[2]])),
                n23 = length(intersect(venn_list[[2]], venn_list[[3]])),
                n13 = length(intersect(venn_list[[1]], venn_list[[3]])),
                n123 = length(Reduce(intersect, venn_list)),
                category = names(venn_list),
                fill = c("#377EB8", "#4DAF4A", "#984EA3"),
                alpha = 0.5,
                cat.pos = c(-20, 20, 180),
                main = "Common Microexon Events"
            )
            grid.draw(venn_plot)
        }

        dev.off()

        # Save common microexon events
        common_microexons <- Reduce(intersect, sig_microexons)

        if (length(common_microexons) > 0) {
            # Get details - remove duplicates by taking first occurrence
            common_details <- data %>%
                filter(size_category == "Microexon (0-30bp)",
                       significant) %>%
                mutate(event_id = paste(GeneID, chr, exonStart_0base, exonEnd, sep = "_")) %>%
                filter(event_id %in% common_microexons) %>%
                select(comparison, event_id, GeneID, geneSymbol, exon_length,
                       IncLevelDifference, FDR) %>%
                distinct(event_id, comparison, .keep_all = TRUE)

            # Pivot to wide format with values_fn to handle any remaining duplicates
            common_wide <- common_details %>%
                pivot_wider(
                    id_cols = c(event_id, GeneID, geneSymbol, exon_length),
                    names_from = comparison,
                    values_from = c(IncLevelDifference, FDR),
                    values_fn = first
                )

            write.csv(common_wide,
                     file.path(output_dir, "common_events", "common_microexons_details.csv"),
                     row.names = FALSE)

            cat(paste0("  Common microexon events: ", length(common_microexons), "\n"))
        }
    }

    # Also do the same for small exons (30-50bp)
    sig_small <- list()

    for (comp in parental_comparisons) {
        events <- data %>%
            filter(comparison == comp,
                   size_category == "Small (30-50bp)",
                   significant) %>%
            mutate(event_id = paste(GeneID, chr, exonStart_0base, exonEnd, sep = "_"))

        sig_small[[comp]] <- unique(events$event_id)
    }

    if (all(sapply(sig_small, length) > 0)) {
        common_small <- Reduce(intersect, sig_small)

        if (length(common_small) > 0) {
            common_details <- data %>%
                filter(size_category == "Small (30-50bp)",
                       significant) %>%
                mutate(event_id = paste(GeneID, chr, exonStart_0base, exonEnd, sep = "_")) %>%
                filter(event_id %in% common_small) %>%
                select(comparison, event_id, GeneID, geneSymbol, exon_length,
                       IncLevelDifference, FDR) %>%
                distinct(event_id, comparison, .keep_all = TRUE)

            common_wide <- common_details %>%
                pivot_wider(
                    id_cols = c(event_id, GeneID, geneSymbol, exon_length),
                    names_from = comparison,
                    values_from = c(IncLevelDifference, FDR),
                    values_fn = first
                )

            write.csv(common_wide,
                     file.path(output_dir, "common_events", "common_small_exons_details.csv"),
                     row.names = FALSE)

            cat(paste0("  Common small exon events: ", length(common_small), "\n"))
        }
    }

    return(sig_microexons)
}

# =============================================================================
# Main Execution
# =============================================================================

# 1. Classify exons by size
cat("\n=== Step 1: Exon Size Classification ===\n")
result <- classify_exons_by_size(SPLICING_DIR, all_comparisons)
all_data <- result$data
size_summary <- result$summary

# Save summary
write.csv(size_summary, file.path(OUTPUT_DIR, "exon_size_summary.csv"), row.names = FALSE)

# 2. Plot size distributions
cat("\n=== Step 2: Size Distribution Plots ===\n")
plot_size_distribution(all_data, OUTPUT_DIR)

# 3. Microexon-specific analysis
cat("\n=== Step 3: Microexon-Specific Analysis ===\n")
sig_microexons <- analyze_microexon_splicing(all_data, OUTPUT_DIR)

# 4. Compare microexon vs regular patterns
cat("\n=== Step 4: Microexon vs Regular Comparison ===\n")
comparison_stats <- compare_microexon_vs_regular(all_data, OUTPUT_DIR)

# 5. Common microexon events (focus on parental comparisons)
cat("\n=== Step 5: Common Microexon Events ===\n")
common_events <- find_common_microexon_events(all_data, parental_comparisons, OUTPUT_DIR)

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n============================================\n")
cat("Microexon Analysis Summary\n")
cat("============================================\n")

# Print summary for parental comparisons
for (comp in parental_comparisons) {
    comp_data <- all_data %>% filter(comparison == comp)
    sig_data <- comp_data %>% filter(significant)

    cat(paste0("\n", comp, ":\n"))
    cat(paste0("  Total SE events: ", nrow(comp_data), "\n"))
    cat(paste0("  Significant SE events: ", nrow(sig_data), "\n"))

    for (size_cat in levels(all_data$size_category)) {
        n_total <- sum(comp_data$size_category == size_cat)
        n_sig <- sum(comp_data$size_category == size_cat & comp_data$significant)
        cat(paste0("    ", size_cat, ": ", n_sig, "/", n_total,
                  " (", round(n_sig/n_total*100, 1), "%)\n"))
    }
}

cat("\n============================================\n")
cat("Microexon Analysis Complete!\n")
cat("============================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - size_distribution/: Exon size histograms and bar charts\n")
cat("  - microexon_specific/: Microexon-focused analyses\n")
cat("  - common_events/: Venn diagrams for microexon overlaps\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
