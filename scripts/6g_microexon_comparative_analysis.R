#!/usr/bin/env Rscript
# =============================================================================
# Microexon Comparative Analysis
# Project: 90-1239779069
# Analysis: Common/unique microexons with directionality across SRRM3 perturbations
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(VennDiagram)
    library(UpSetR)
    library(grid)
    library(RColorBrewer)
    library(pheatmap)
})

# Try to load ggalluvial for Sankey-like plots (optional)
sankey_available <- FALSE # Forced to FALSE to ensure pipeline completion

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/09_microexon_extended/comparative_analysis")

# Thresholds
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1

# Exon size bins
MICROEXON_MAX <- 30
SMALL_EXON_MAX <- 50

cat("============================================\n")
cat("Microexon Comparative Analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define comparisons
parental_comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")

# Color palette
comparison_colors <- c("Neg_vs_Parental" = "#377EB8",
                       "Pos_vs_Parental" = "#4DAF4A",
                       "KO_vs_Parental" = "#984EA3")
direction_colors <- c("Inclusion" = "#e74c3c", "Exclusion" = "#3498db")

# =============================================================================
# Function: Load SE events with size and direction
# =============================================================================
load_se_data <- function(comparison, splicing_dir) {
    file_path <- file.path(splicing_dir, comparison, "SE.MATS.JC.txt")

    if (!file.exists(file_path)) {
        warning(paste("File not found:", file_path))
        return(NULL)
    }

    df <- read.delim(file_path, stringsAsFactors = FALSE)

    # Calculate exon length
    df$exon_length <- df$exonEnd - df$exonStart_0base

    # Classify by size
    df$size_category <- case_when(
        df$exon_length <= MICROEXON_MAX ~ "Microexon (0-30bp)",
        df$exon_length <= SMALL_EXON_MAX ~ "Small (30-50bp)",
        TRUE ~ "Regular (>50bp)"
    )
    df$size_category <- factor(df$size_category,
                               levels = c("Microexon (0-30bp)", "Small (30-50bp)", "Regular (>50bp)"))

    # Mark significant and direction
    df$significant <- df$FDR < FDR_THRESHOLD & abs(df$IncLevelDifference) > DPSI_THRESHOLD
    df$dPSI <- df$IncLevelDifference
    df$direction <- ifelse(df$IncLevelDifference > 0, "Inclusion", "Exclusion")
    df$comparison <- comparison

    # Create unique event ID
    # Create unique event ID
    df$event_id <- paste(df$GeneID, df$chr, df$exonStart_0base, df$exonEnd, sep = "_")
    
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
# Function: Build dPSI comparison matrix
# =============================================================================
build_dPSI_comparison_matrix <- function(all_data, comparisons) {
    cat("Building dPSI comparison matrix...\n")

    # Get all unique events
    all_events <- unique(do.call(rbind, lapply(all_data, function(x) {
        x %>% select(event_id, GeneID, geneSymbol, chr, exonStart_0base, exonEnd, exon_length, size_category)
    })))

    # Remove duplicates by event_id
    all_events <- all_events[!duplicated(all_events$event_id), ]

    # Build matrix
    for (comp in comparisons) {
        comp_data <- all_data[[comp]] %>%
            select(event_id, dPSI, FDR, significant, direction) %>%
            rename_with(~paste0(gsub("_vs_Parental", "", comp), "_", .), -event_id)

        all_events <- merge(all_events, comp_data, by = "event_id", all.x = TRUE)
    }

    cat(paste0("  Total unique events: ", nrow(all_events), "\n"))

    return(all_events)
}

# =============================================================================
# Function: Classify direction patterns
# =============================================================================
classify_direction_patterns <- function(comparison_matrix, comparisons) {
    cat("Classifying direction patterns...\n")

    # Get direction columns
    dir_cols <- paste0(gsub("_vs_Parental", "", comparisons), "_direction")

    # Only consider events significant in at least one comparison
    sig_cols <- paste0(gsub("_vs_Parental", "", comparisons), "_significant")
    sig_any <- rowSums(comparison_matrix[, sig_cols], na.rm = TRUE) > 0

    sig_matrix <- comparison_matrix[sig_any, ]

    # Create pattern string
    sig_matrix$pattern <- apply(sig_matrix[, dir_cols], 1, function(row) {
        pattern <- ifelse(is.na(row), "NS", substr(row, 1, 3))  # Inc, Exc, or NS (not sig)
        paste(pattern, collapse = "/")
    })

    # Count patterns
    pattern_counts <- sig_matrix %>%
        group_by(pattern) %>%
        summarize(
            count = n(),
            n_microexons = sum(size_category == "Microexon (0-30bp)", na.rm = TRUE),
            n_small = sum(size_category == "Small (30-50bp)", na.rm = TRUE),
            n_regular = sum(size_category == "Regular (>50bp)", na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(count))

    cat("\nTop direction patterns:\n")
    print(head(pattern_counts, 10))

    return(list(matrix = sig_matrix, patterns = pattern_counts))
}

# =============================================================================
# Function: Find common microexons
# =============================================================================
find_common_microexons <- function(all_data, comparisons, size_filter = NULL) {
    cat("Finding common microexon events...\n")

    # Get significant events for each comparison
    sig_events <- list()
    for (comp in comparisons) {
        data <- all_data[[comp]]
        if (!is.null(size_filter)) {
            data <- data %>% filter(size_category == size_filter)
        }
        sig_events[[comp]] <- unique(data$event_id[data$significant])
        cat(paste0("  ", comp, ": ", length(sig_events[[comp]]), " significant events\n"))
    }

    # Find intersections
    common_all <- Reduce(intersect, sig_events)
    common_neg_ko <- intersect(sig_events[["Neg_vs_Parental"]], sig_events[["KO_vs_Parental"]])
    common_neg_pos <- intersect(sig_events[["Neg_vs_Parental"]], sig_events[["Pos_vs_Parental"]])
    common_pos_ko <- intersect(sig_events[["Pos_vs_Parental"]], sig_events[["KO_vs_Parental"]])

    cat(paste0("\n  Common to all: ", length(common_all), "\n"))
    cat(paste0("  Common Neg & KO: ", length(common_neg_ko), "\n"))
    cat(paste0("  Common Neg & Pos: ", length(common_neg_pos), "\n"))
    cat(paste0("  Common Pos & KO: ", length(common_pos_ko), "\n"))

    return(list(
        all_three = common_all,
        neg_ko = common_neg_ko,
        neg_pos = common_neg_pos,
        pos_ko = common_pos_ko,
        by_comparison = sig_events
    ))
}

# =============================================================================
# Function: Find unique microexons
# =============================================================================
find_unique_microexons <- function(common_results, all_data, comparisons) {
    cat("Finding comparison-specific microexon events...\n")

    sig_events <- common_results$by_comparison
    unique_events <- list()

    for (comp in comparisons) {
        others <- setdiff(comparisons, comp)
        other_events <- unique(unlist(sig_events[others]))
        unique_events[[comp]] <- setdiff(sig_events[[comp]], other_events)
        cat(paste0("  Unique to ", comp, ": ", length(unique_events[[comp]]), "\n"))
    }

    return(unique_events)
}

# =============================================================================
# Function: Analyze concordance
# =============================================================================
analyze_concordance <- function(comparison_matrix, comparisons) {
    cat("Analyzing directional concordance...\n")

    # Get significant columns
    sig_cols <- paste0(gsub("_vs_Parental", "", comparisons), "_significant")
    dir_cols <- paste0(gsub("_vs_Parental", "", comparisons), "_direction")

    results <- list()

    # Pairwise comparisons
    pairs <- list(
        c("Neg_vs_Parental", "KO_vs_Parental"),  # Both loss-of-function: expect same direction
        c("Neg_vs_Parental", "Pos_vs_Parental"), # Opposite: expect different direction
        c("KO_vs_Parental", "Pos_vs_Parental")   # Opposite: expect different direction
    )

    for (pair in pairs) {
        comp1 <- gsub("_vs_Parental", "", pair[1])
        comp2 <- gsub("_vs_Parental", "", pair[2])

        # Events significant in both
        both_sig <- comparison_matrix[[paste0(comp1, "_significant")]] &
                    comparison_matrix[[paste0(comp2, "_significant")]]
        both_sig[is.na(both_sig)] <- FALSE

        if (sum(both_sig) > 0) {
            subset <- comparison_matrix[both_sig, ]

            same_direction <- sum(subset[[paste0(comp1, "_direction")]] ==
                                 subset[[paste0(comp2, "_direction")]], na.rm = TRUE)
            opposite_direction <- sum(subset[[paste0(comp1, "_direction")]] !=
                                     subset[[paste0(comp2, "_direction")]], na.rm = TRUE)
            total <- same_direction + opposite_direction

            results[[paste(pair, collapse = "_vs_")]] <- data.frame(
                comparison1 = pair[1],
                comparison2 = pair[2],
                total_common = total,
                same_direction = same_direction,
                opposite_direction = opposite_direction,
                pct_concordant = round(same_direction / total * 100, 1)
            )

            cat(paste0("\n  ", pair[1], " vs ", pair[2], ":\n"))
            cat(paste0("    Common significant: ", total, "\n"))
            cat(paste0("    Same direction: ", same_direction, " (", round(same_direction/total*100, 1), "%)\n"))
            cat(paste0("    Opposite direction: ", opposite_direction, " (", round(opposite_direction/total*100, 1), "%)\n"))
        }
    }

    return(do.call(rbind, results))
}

# =============================================================================
# Function: Plot UpSet diagram
# =============================================================================
plot_upset_microexons <- function(all_data, comparisons, output_dir, size_filter = NULL) {
    cat("Generating UpSet plot...\n")

    # Create binary matrix for UpSetR
    sig_events <- list()
    for (comp in comparisons) {
        data <- all_data[[comp]]
        if (!is.null(size_filter)) {
            data <- data %>% filter(size_category == size_filter)
        }
        sig_events[[gsub("_vs_Parental", "", comp)]] <- unique(data$event_id[data$significant])
    }

    # Get all events
    all_event_ids <- unique(unlist(sig_events))

    if (length(all_event_ids) == 0) {
        warning("No significant events for UpSet plot")
        return(NULL)
    }

    # Create binary matrix
    upset_data <- data.frame(event_id = all_event_ids)
    for (comp_name in names(sig_events)) {
        upset_data[[comp_name]] <- as.integer(upset_data$event_id %in% sig_events[[comp_name]])
    }

    # Generate plot
    suffix <- ifelse(is.null(size_filter), "all_events", gsub("[^a-zA-Z]", "", size_filter))

    pdf(file.path(output_dir, paste0("upset_plot_", suffix, ".pdf")), width = 10, height = 7)
    print(
        upset(
            upset_data,
            sets = names(sig_events),
            order.by = "freq",
            main.bar.color = "#3498db",
            sets.bar.color = c("#377EB8", "#4DAF4A", "#984EA3"),
            text.scale = 1.3,
            mb.ratio = c(0.6, 0.4)
        )
    )
    dev.off()
}

# =============================================================================
# Function: Plot direction flow (alluvial/Sankey-like)
# =============================================================================
plot_direction_flow <- function(comparison_matrix, output_dir) {
    cat("Generating direction flow plot...\n")

    # Filter to events significant in at least 2 comparisons
    sig_cols <- c("Neg_significant", "Pos_significant", "KO_significant")
    dir_cols <- c("Neg_direction", "Pos_direction", "KO_direction")

    matrix_filtered <- comparison_matrix
    matrix_filtered$n_sig <- rowSums(matrix_filtered[, sig_cols], na.rm = TRUE)
    matrix_filtered <- matrix_filtered %>% filter(n_sig >= 2)

    if (nrow(matrix_filtered) == 0) {
        warning("No events significant in 2+ comparisons for direction flow")
        return(NULL)
    }

    # Create long format for visualization
    flow_data <- matrix_filtered %>%
        select(event_id, size_category, Neg_direction, Pos_direction, KO_direction) %>%
        mutate(
            Neg = ifelse(is.na(Neg_direction), "NS", Neg_direction),
            Pos = ifelse(is.na(Pos_direction), "NS", Pos_direction),
            KO = ifelse(is.na(KO_direction), "NS", KO_direction)
        )

    # Count flow patterns
    flow_counts <- flow_data %>%
        group_by(Neg, Pos, KO, size_category) %>%
        summarize(count = n(), .groups = "drop")

    if (sankey_available) {
        # Limit to top 500 events to prevent rendering issues and illegible plots
        if(nrow(flow_data) > 500) {
            cat(paste0("  (Subsampling 500 events from ", nrow(flow_data), " for Sankey plot)\n"))
            set.seed(123)
            flow_data <- flow_data %>% sample_n(500)
        }

        # Use ggalluvial for Sankey-like plot
        flow_long <- flow_data %>%
            select(event_id, size_category, Neg, Pos, KO) %>%
            pivot_longer(cols = c(Neg, Pos, KO), names_to = "comparison", values_to = "direction") %>%
            mutate(
                comparison = factor(comparison, levels = c("Neg", "Pos", "KO")),
                direction = factor(direction, levels = c("Inclusion", "Exclusion", "NS")),
                event_id = as.factor(event_id)
            )

        # Plot using event_id as alluvium
        p <- ggplot(flow_long, aes(x = comparison, stratum = direction, alluvium = event_id,
                                    fill = direction, label = direction)) +
            geom_flow(alpha = 0.6) +
            geom_stratum(alpha = 0.8) +
            geom_text(stat = "stratum", size = 3) +
            scale_fill_manual(values = c("Inclusion" = "#e74c3c", "Exclusion" = "#3498db", "NS" = "#bdc3c7")) +
            theme_minimal() +
            labs(
                x = "Comparison (vs Parental)",
                y = "Number of Events",
                title = "Direction Flow Across SRRM3 Perturbations",
                subtitle = "Events significant in at least 2 comparisons",
                fill = "Direction"
            )

        ggsave(file.path(output_dir, "direction_sankey.pdf"), p, width = 10, height = 8)
    } else {
        # Fallback: use heatmap-style visualization
        cat("  (ggalluvial not available, using alternative visualization)\n")

        pattern_summary <- flow_counts %>%
            mutate(pattern = paste(Neg, Pos, KO, sep = "/")) %>%
            group_by(pattern) %>%
            summarize(total = sum(count), .groups = "drop") %>%
            arrange(desc(total))

        p <- ggplot(head(pattern_summary, 15), aes(x = reorder(pattern, total), y = total)) +
            geom_bar(stat = "identity", fill = "#3498db") +
            coord_flip() +
            theme_bw() +
            labs(
                x = "Direction Pattern (Neg/Pos/KO)",
                y = "Number of Events",
                title = "Direction Patterns Across Comparisons",
                subtitle = "Inc = Inclusion, Exc = Exclusion, NS = Not Significant"
            )

        ggsave(file.path(output_dir, "direction_patterns.pdf"), p, width = 10, height = 8)
    }
}

# =============================================================================
# Function: Plot dPSI heatmap for common events
# =============================================================================
plot_dPSI_heatmap <- function(comparison_matrix, common_events, output_dir) {
    cat("Generating dPSI heatmap...\n")

    if (length(common_events) < 5) {
        warning("Too few common events for heatmap")
        return(NULL)
    }

    # Filter to common events
    common_data <- comparison_matrix %>%
        filter(event_id %in% common_events)

    # Create matrix for heatmap
    dpsi_matrix <- common_data %>%
        select(geneSymbol, Neg_dPSI, Pos_dPSI, KO_dPSI)

    # Use gene symbol as row names (or event_id if duplicates)
    row_labels <- common_data$geneSymbol
    if (any(duplicated(row_labels))) {
        row_labels <- paste(common_data$geneSymbol, seq_along(row_labels), sep = "_")
    }
    rownames(dpsi_matrix) <- row_labels
    dpsi_matrix <- as.matrix(dpsi_matrix[, c("Neg_dPSI", "Pos_dPSI", "KO_dPSI")])
    colnames(dpsi_matrix) <- c("Neg", "Pos", "KO")

    # Remove rows with NA
    dpsi_matrix <- dpsi_matrix[complete.cases(dpsi_matrix), ]

    if (nrow(dpsi_matrix) < 3) {
        warning("Too few complete cases for heatmap")
        return(NULL)
    }

    # Limit to top 50 by variance if too many
    if (nrow(dpsi_matrix) > 50) {
        row_vars <- apply(dpsi_matrix, 1, var)
        dpsi_matrix <- dpsi_matrix[order(row_vars, decreasing = TRUE)[1:50], ]
    }

    # Create heatmap
    pdf(file.path(output_dir, "dPSI_heatmap_common.pdf"), width = 8, height = 12)
    pheatmap(
        dpsi_matrix,
        main = "dPSI Values for Common Significant Events",
        color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(50),
        breaks = seq(-0.5, 0.5, length.out = 51),
        cluster_cols = FALSE,
        fontsize_row = 8,
        border_color = NA
    )
    dev.off()
}

# =============================================================================
# Main Execution
# =============================================================================

# 1. Load all data
cat("\n=== Step 1: Loading SE Data ===\n")
all_data <- list()
for (comp in parental_comparisons) {
    cat(paste0("Loading ", comp, "...\n"))
    data <- load_se_data(comp, SPLICING_DIR)
    if (!is.null(data)) {
        all_data[[comp]] <- data
    }
}

# 2. Build comparison matrix
cat("\n=== Step 2: Building Comparison Matrix ===\n")
comparison_matrix <- build_dPSI_comparison_matrix(all_data, parental_comparisons)
write.csv(comparison_matrix, file.path(OUTPUT_DIR, "dPSI_comparison_matrix.csv"), row.names = FALSE)

# 3. Classify direction patterns
cat("\n=== Step 3: Classifying Direction Patterns ===\n")
pattern_results <- classify_direction_patterns(comparison_matrix, parental_comparisons)
write.csv(pattern_results$patterns, file.path(OUTPUT_DIR, "directionality_patterns.csv"), row.names = FALSE)
write.csv(pattern_results$matrix, file.path(OUTPUT_DIR, "events_with_patterns.csv"), row.names = FALSE)

# 4. Find common and unique events
cat("\n=== Step 4: Finding Common/Unique Events ===\n")

# All events
common_all <- find_common_microexons(all_data, parental_comparisons)

# Microexons only
common_micro <- find_common_microexons(all_data, parental_comparisons, "Microexon (0-30bp)")

# Find unique events
unique_all <- find_unique_microexons(common_all, all_data, parental_comparisons)
unique_micro <- find_unique_microexons(common_micro, all_data, parental_comparisons)

# Save common events with details
common_details <- comparison_matrix %>%
    filter(event_id %in% common_all$all_three) %>%
    arrange(geneSymbol)
write.csv(common_details, file.path(OUTPUT_DIR, "common_microexons_all.csv"), row.names = FALSE)

# Save unique events
for (comp in parental_comparisons) {
    if (length(unique_all[[comp]]) > 0) {
        unique_details <- all_data[[comp]] %>%
            filter(event_id %in% unique_all[[comp]], significant) %>%
            select(event_id, GeneID, geneSymbol, chr, exonStart_0base, exonEnd,
                   exon_length, size_category, dPSI, FDR, direction)
        write.csv(unique_details, file.path(OUTPUT_DIR, paste0("unique_", comp, ".csv")), row.names = FALSE)
    }
}

# 5. Analyze concordance
cat("\n=== Step 5: Analyzing Concordance ===\n")
concordance_results <- analyze_concordance(comparison_matrix, parental_comparisons)
if (!is.null(concordance_results)) {
    write.csv(concordance_results, file.path(OUTPUT_DIR, "concordance_analysis.csv"), row.names = FALSE)
}

# 6. Generate plots
cat("\n=== Step 6: Generating Plots ===\n")

# UpSet plots
plot_upset_microexons(all_data, parental_comparisons, OUTPUT_DIR)
plot_upset_microexons(all_data, parental_comparisons, OUTPUT_DIR, "Microexon (0-30bp)")

# Direction flow
plot_direction_flow(comparison_matrix, OUTPUT_DIR)

# dPSI heatmap
plot_dPSI_heatmap(comparison_matrix, common_all$all_three, OUTPUT_DIR)

# 7. Venn diagram
cat("Generating Venn diagram...\n")
venn_list <- common_micro$by_comparison
names(venn_list) <- gsub("_vs_Parental", "", names(venn_list))

if (all(sapply(venn_list, length) > 0)) {
    pdf(file.path(OUTPUT_DIR, "microexon_venn.pdf"), width = 10, height = 10)
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
        cat.cex = 1.5,
        cex = 1.3,
        main = "Significant Microexons (0-30bp)"
    )
    grid.draw(venn_plot)
    dev.off()
}

# 8. Comparison bar plot
cat("Generating comparison summary plot...\n")

summary_data <- data.frame()
for (comp in parental_comparisons) {
    for (size in c("Microexon (0-30bp)", "Small (30-50bp)", "Regular (>50bp)")) {
        for (dir in c("Inclusion", "Exclusion")) {
            count <- sum(all_data[[comp]]$size_category == size &
                        all_data[[comp]]$significant &
                        all_data[[comp]]$direction == dir)
            summary_data <- rbind(summary_data, data.frame(
                comparison = comp,
                size_category = size,
                direction = dir,
                count = count
            ))
        }
    }
}

p_summary <- ggplot(summary_data, aes(x = comparison, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~size_category) +
    scale_fill_manual(values = direction_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        x = "Comparison",
        y = "Number of Significant Events",
        title = "Splicing Direction by Comparison and Exon Size",
        fill = "Direction"
    )

ggsave(file.path(OUTPUT_DIR, "comparison_summary.pdf"), p_summary, width = 14, height = 8)

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n============================================\n")
cat("Comparative Analysis Summary\n")
cat("============================================\n")

cat("\nSignificant events per comparison:\n")
for (comp in parental_comparisons) {
    total_sig <- sum(all_data[[comp]]$significant)
    micro_sig <- sum(all_data[[comp]]$significant & all_data[[comp]]$size_category == "Microexon (0-30bp)")
    inc <- sum(all_data[[comp]]$significant & all_data[[comp]]$direction == "Inclusion")
    exc <- sum(all_data[[comp]]$significant & all_data[[comp]]$direction == "Exclusion")

    cat(paste0("\n", comp, ":\n"))
    cat(paste0("  Total significant: ", total_sig, "\n"))
    cat(paste0("  Microexons: ", micro_sig, "\n"))
    cat(paste0("  Inclusion: ", inc, " (", round(inc/total_sig*100, 1), "%)\n"))
    cat(paste0("  Exclusion: ", exc, " (", round(exc/total_sig*100, 1), "%)\n"))
}

cat("\nCommon events:\n")
cat(paste0("  All three comparisons: ", length(common_all$all_three), "\n"))
cat(paste0("  Neg & KO (both loss-of-function): ", length(common_all$neg_ko), "\n"))

if (!is.null(concordance_results)) {
    cat("\nConcordance (Neg vs KO - both loss-of-function):\n")
    neg_ko_row <- concordance_results %>% filter(comparison1 == "Neg_vs_Parental", comparison2 == "KO_vs_Parental")
    if (nrow(neg_ko_row) > 0) {
        cat(paste0("  Same direction: ", neg_ko_row$same_direction, " (",
                  neg_ko_row$pct_concordant, "%)\n"))
    }
}

cat("\n============================================\n")
cat("Comparative Analysis Complete!\n")
cat("============================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - dPSI_comparison_matrix.csv: Wide-format dPSI across comparisons\n")
cat("  - common_microexons_all.csv: Events significant in all comparisons\n")
cat("  - unique_*.csv: Comparison-specific events\n")
cat("  - directionality_patterns.csv: Direction pattern counts\n")
cat("  - concordance_analysis.csv: Pairwise concordance\n")
cat("  - upset_plot_*.pdf: UpSet diagrams\n")
cat("  - direction_*.pdf: Direction flow visualization\n")
cat("  - microexon_venn.pdf: Venn diagram\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
