#!/usr/bin/env Rscript
# =============================================================================
# Splicing Downstream Analysis
# Project: 90-1239779069
# Analysis: PCA on PSI, event distributions, directionality, cross-plots
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
    library(patchwork)
    library(VennDiagram)
    library(grid)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/06_splicing_downstream")

# Thresholds (matching paper standards)
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1

cat("============================================\n")
cat("Splicing Downstream Analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Thresholds: FDR <", FDR_THRESHOLD, ", |dPSI| >", DPSI_THRESHOLD, "\n")
cat("============================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "PCA"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "distributions"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "crossplots"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "common_events"), showWarnings = FALSE)

# Sample metadata
sample_info <- data.frame(
    sample_id = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15),
    sample_name = c("Parental_1", "Parental_2", "Parental_3",
                    "Neg_1", "Neg_2", "Neg_3",
                    "Pos_1", "Pos_2", "Pos_3",
                    "KO_1", "KO_2", "KO_3"),
    group = factor(c(rep("Parental", 3), rep("Neg", 3),
                     rep("Pos", 3), rep("KO", 3)),
                   levels = c("Parental", "Neg", "Pos", "KO"))
)

# Define comparisons (focus on vs Parental as requested)
comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental",
                 "Pos_vs_Neg", "KO_vs_Neg", "KO_vs_Pos")
parental_comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")
event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

# Color palettes
group_colors <- c(Parental = "#E41A1C", Neg = "#377EB8", Pos = "#4DAF4A", KO = "#984EA3")
event_colors <- c(SE = "#1f77b4", A5SS = "#ff7f0e", A3SS = "#2ca02c",
                  MXE = "#d62728", RI = "#9467bd")

# =============================================================================
# Function: Load rMATS results
# =============================================================================
load_rmats_results <- function(comparison, event_type, splicing_dir) {
    # Use JC (junction counts) file
    file_path <- file.path(splicing_dir, comparison, paste0(event_type, ".MATS.JC.txt"))
    if (!file.exists(file_path)) {
        warning(paste("File not found:", file_path))
        return(NULL)
    }
    df <- read.delim(file_path, stringsAsFactors = FALSE)
    df$event_type <- event_type
    df$comparison <- comparison
    # Calculate average coverage (IJC + SJC) across all samples
    # IJC/SJC columns are comma-separated strings of counts per sample
    # We need to parse them to filter
    
    # Function to sum counts from comma-separated string
    sum_counts <- function(x) {
        if (is.null(x) || is.na(x)) return(0)
        sum(as.numeric(unlist(strsplit(as.character(x), ","))), na.rm = TRUE)
    }
    
    # Calculate total reads per event (sum of all samples in both groups)
    # Note: simple approach to avoid complex parsing for every row if unnecessary
    # But for accurate filtering we need to parse.
    
    # Add read count filter if columns exist
    if(all(c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2") %in% colnames(df))) {
        # Calculate row-wise average read coverage
        # (Total IJC + Total SJC) / (N_samples_1 + N_samples_2)
        # Assuming standard 3+3 design = 6 samples. If unknown, we estimate.
        
        # Vectorized calculation for efficiency
        df$avg_reads <- (
            sapply(df$IJC_SAMPLE_1, sum_counts) + 
            sapply(df$SJC_SAMPLE_1, sum_counts) + 
            sapply(df$IJC_SAMPLE_2, sum_counts) + 
            sapply(df$SJC_SAMPLE_2, sum_counts)
        ) / 6 # Normalizing by approximate total samples (3+3)
        
        # Apply filter: Average reads > 10
        # This roughly means ~60 reads total across 6 samples
        filtered_df <- df[df$avg_reads >= 10, ]
        
        n_orig <- nrow(df)
        n_filt <- nrow(filtered_df)
        if(n_orig > n_filt) {
            # silently filter, or log if needed
            # cat(paste0("  Filtered ", n_orig - n_filt, " low-count events\n"))
        }
        df <- filtered_df
    }
    
    return(df)
}

# =============================================================================
# Helpers for coordinate-based event matching (rMATS IDs are NOT consistent
# across separate runs — only genomic coordinates are reliable)
# =============================================================================
make_event_key <- function(df) {
    paste(df$chr, df$strand, df$exonStart_0base, df$exonEnd,
          df$upstreamEE, df$downstreamES, sep = ":")
}

# Vectorized parsing of comma-separated value columns into a numeric matrix
parse_csv_column <- function(x, expected_n = 3) {
    splits <- strsplit(as.character(x), ",")
    lens   <- vapply(splits, length, integer(1))
    mat    <- matrix(NA_real_, nrow = length(x), ncol = expected_n)
    valid  <- lens == expected_n
    if (any(valid)) {
        mat[valid, ] <- do.call(rbind, lapply(splits[valid], as.numeric))
    }
    mat
}

# =============================================================================
# Function: Extract PSI matrix for PCA
# Uses coordinate-based keys (not rMATS IDs) to match events across comparisons
# =============================================================================
extract_psi_matrix <- function(splicing_dir, event_types) {
    cat("Extracting PSI matrix for PCA (coordinate-based matching)...\n")

    # --- Load all event types from the three vs-Parental comparisons ---
    # Neg_vs_Parental: IncLevel1 = Neg (3 samples), IncLevel2 = Parental (3 samples)
    # Pos_vs_Parental: IncLevel1 = Pos (3 samples), IncLevel2 = Parental (3 samples)
    # KO_vs_Parental:  IncLevel1 = KO  (3 samples), IncLevel2 = Parental (3 samples)

    neg_keys_all <- character(); neg_par_all <- list(); neg_grp_all <- list()
    pos_keys_all <- character(); pos_grp_all <- list()
    ko_keys_all  <- character(); ko_grp_all  <- list()

    for (event in event_types) {
        # --- Neg_vs_Parental (Parental + Neg PSI) ---
        neg_file <- file.path(splicing_dir, "Neg_vs_Parental", paste0(event, ".MATS.JC.txt"))
        if (file.exists(neg_file)) {
            df <- read.delim(neg_file, stringsAsFactors = FALSE)
            keys <- paste0(event, ":", make_event_key(df))
            inc1 <- parse_csv_column(df$IncLevel1, 3)  # Neg
            inc2 <- parse_csv_column(df$IncLevel2, 3)  # Parental
            valid <- complete.cases(inc1) & complete.cases(inc2)
            neg_keys_all <- c(neg_keys_all, keys[valid])
            neg_par_all  <- c(neg_par_all, split(inc2[valid, , drop = FALSE], seq_len(sum(valid))))
            neg_grp_all  <- c(neg_grp_all, split(inc1[valid, , drop = FALSE], seq_len(sum(valid))))
        }

        # --- Pos_vs_Parental (Pos PSI only) ---
        pos_file <- file.path(splicing_dir, "Pos_vs_Parental", paste0(event, ".MATS.JC.txt"))
        if (file.exists(pos_file)) {
            df <- read.delim(pos_file, stringsAsFactors = FALSE)
            keys <- paste0(event, ":", make_event_key(df))
            inc1 <- parse_csv_column(df$IncLevel1, 3)  # Pos
            valid <- complete.cases(inc1)
            pos_keys_all <- c(pos_keys_all, keys[valid])
            pos_grp_all  <- c(pos_grp_all, split(inc1[valid, , drop = FALSE], seq_len(sum(valid))))
        }

        # --- KO_vs_Parental (KO PSI only) ---
        ko_file <- file.path(splicing_dir, "KO_vs_Parental", paste0(event, ".MATS.JC.txt"))
        if (file.exists(ko_file)) {
            df <- read.delim(ko_file, stringsAsFactors = FALSE)
            keys <- paste0(event, ":", make_event_key(df))
            inc1 <- parse_csv_column(df$IncLevel1, 3)  # KO
            valid <- complete.cases(inc1)
            ko_keys_all <- c(ko_keys_all, keys[valid])
            ko_grp_all  <- c(ko_grp_all, split(inc1[valid, , drop = FALSE], seq_len(sum(valid))))
        }
    }

    cat(sprintf("  Loaded events — Neg_vs_Par: %d, Pos_vs_Par: %d, KO_vs_Par: %d\n",
                length(neg_keys_all), length(pos_keys_all), length(ko_keys_all)))

    # --- Find events common to all three comparisons by coordinate key ---
    common_keys <- Reduce(intersect, list(neg_keys_all, pos_keys_all, ko_keys_all))
    cat(sprintf("  Common events (coordinate-matched): %d\n", length(common_keys)))

    if (length(common_keys) < 100) {
        warning("Less than 100 common events found")
        return(NULL)
    }

    # --- Build named lookups ---
    names(neg_par_all) <- neg_keys_all
    names(neg_grp_all) <- neg_keys_all
    names(pos_grp_all) <- pos_keys_all
    names(ko_grp_all)  <- ko_keys_all

    # --- Assemble 12-sample PSI matrix ---
    psi_matrix <- matrix(NA, nrow = length(common_keys), ncol = 12)
    rownames(psi_matrix) <- common_keys
    colnames(psi_matrix) <- c("Parental_1", "Parental_2", "Parental_3",
                               "Neg_1", "Neg_2", "Neg_3",
                               "Pos_1", "Pos_2", "Pos_3",
                               "KO_1", "KO_2", "KO_3")

    for (i in seq_along(common_keys)) {
        k <- common_keys[i]
        psi_matrix[i, 1:3]   <- neg_par_all[[k]]   # Parental (from Neg_vs_Par IncLevel2)
        psi_matrix[i, 4:6]   <- neg_grp_all[[k]]   # Neg      (from Neg_vs_Par IncLevel1)
        psi_matrix[i, 7:9]   <- pos_grp_all[[k]]   # Pos      (from Pos_vs_Par IncLevel1)
        psi_matrix[i, 10:12] <- ko_grp_all[[k]]    # KO       (from KO_vs_Par  IncLevel1)
    }

    # Remove rows with any NA
    psi_matrix <- psi_matrix[complete.cases(psi_matrix), ]

    cat("PSI matrix dimensions:", nrow(psi_matrix), "events x", ncol(psi_matrix), "samples\n")
    return(psi_matrix)
}

# =============================================================================
# Function: Plot Splicing PCA
# =============================================================================
plot_splicing_pca <- function(psi_matrix, sample_info, output_dir) {
    cat("Generating PCA plot...\n")

    # Filter events with variance > 0
    vars <- apply(psi_matrix, 1, var, na.rm = TRUE)
    psi_filtered <- psi_matrix[vars > 0.001, ]

    if (nrow(psi_filtered) < 50) {
        warning("Not enough variable events for PCA")
        return(NULL)
    }

    # Perform PCA
    pca_result <- prcomp(t(psi_filtered), center = TRUE, scale. = TRUE)

    # Extract variance explained
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

    # Create plot data
    pca_data <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        sample = rownames(pca_result$x),
        group = sample_info$group[match(rownames(pca_result$x), sample_info$sample_name)]
    )

    # PCA plot
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
        geom_point(size = 4) +
        geom_text_repel(size = 3, max.overlaps = 20) +
        scale_color_manual(values = group_colors) +
        theme_bw() +
        theme(legend.position = "right") +
        labs(
            x = paste0("PC1 (", var_explained[1], "% variance)"),
            y = paste0("PC2 (", var_explained[2], "% variance)"),
            title = "PCA of Splicing Events (PSI values)",
            subtitle = paste0("Based on ", nrow(psi_filtered), " variable events")
        )

    ggsave(file.path(output_dir, "PCA", "splicing_PCA.pdf"), p, width = 10, height = 8)
    ggsave(file.path(output_dir, "PCA", "splicing_PCA.png"), p, width = 10, height = 8, dpi = 300)

    # Also save PC3 vs PC4 if available
    if (ncol(pca_result$x) >= 4) {
        pca_data$PC3 <- pca_result$x[, 3]
        pca_data$PC4 <- pca_result$x[, 4]

        p2 <- ggplot(pca_data, aes(x = PC3, y = PC4, color = group, label = sample)) +
            geom_point(size = 4) +
            geom_text_repel(size = 3, max.overlaps = 20) +
            scale_color_manual(values = group_colors) +
            theme_bw() +
            labs(
                x = paste0("PC3 (", var_explained[3], "% variance)"),
                y = paste0("PC4 (", var_explained[4], "% variance)"),
                title = "PCA of Splicing Events (PC3 vs PC4)"
            )

        ggsave(file.path(output_dir, "PCA", "splicing_PCA_PC3_PC4.pdf"), p2, width = 10, height = 8)
    }

    # Save PCA results
    saveRDS(pca_result, file.path(output_dir, "PCA", "pca_result.rds"))
    write.csv(pca_data, file.path(output_dir, "PCA", "pca_coordinates.csv"), row.names = FALSE)

    return(pca_result)
}

# =============================================================================
# Function: Plot Event Distribution
# =============================================================================
plot_event_distribution <- function(splicing_dir, comparisons, event_types, output_dir) {
    cat("Generating event distribution plots...\n")

    # Collect all events - only keep common columns needed for analysis
    all_events <- list()

    for (comp in comparisons) {
        for (event in event_types) {
            df <- load_rmats_results(comp, event, splicing_dir)
            if (!is.null(df) && nrow(df) > 0) {
                # Select only the common columns needed for analysis
                df_subset <- df %>%
                    select(ID, GeneID, geneSymbol, FDR, IncLevelDifference, event_type, comparison)
                all_events[[paste(comp, event, sep = "_")]] <- df_subset
            }
        }
    }

    # Combine all events (now they all have the same columns)
    combined <- do.call(rbind, all_events)

    # 1. Total events by type and comparison
    summary_df <- combined %>%
        group_by(comparison, event_type) %>%
        summarize(
            total = n(),
            significant = sum(FDR < FDR_THRESHOLD & abs(IncLevelDifference) > DPSI_THRESHOLD, na.rm = TRUE),
            .groups = "drop"
        )

    # Stacked bar plot - total events
    p1 <- ggplot(summary_df, aes(x = comparison, y = total, fill = event_type)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = event_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Number of Events",
            title = "Distribution of Splicing Events by Type",
            fill = "Event Type"
        )

    ggsave(file.path(output_dir, "distributions", "event_distribution_total.pdf"),
           p1, width = 12, height = 8)

    # 2. Significant events only
    p2 <- ggplot(summary_df, aes(x = comparison, y = significant, fill = event_type)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = event_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Number of Significant Events",
            title = paste0("Significant Splicing Events (FDR<", FDR_THRESHOLD,
                          ", |dPSI|>", DPSI_THRESHOLD, ")"),
            fill = "Event Type"
        )

    ggsave(file.path(output_dir, "distributions", "event_distribution_significant.pdf"),
           p2, width = 12, height = 8)

    # 3. Pie chart for each comparison
    for (comp in comparisons) {
        comp_data <- summary_df %>%
            filter(comparison == comp) %>%
            mutate(percentage = significant / sum(significant) * 100)

        if (sum(comp_data$significant) > 0) {
            p <- ggplot(comp_data, aes(x = "", y = significant, fill = event_type)) +
                geom_bar(stat = "identity", width = 1) +
                coord_polar("y", start = 0) +
                scale_fill_manual(values = event_colors) +
                theme_void() +
                labs(
                    title = paste0("Event Type Distribution: ", comp),
                    subtitle = paste0("Total significant events: ", sum(comp_data$significant)),
                    fill = "Event Type"
                ) +
                geom_text(aes(label = ifelse(significant > 0,
                                             paste0(event_type, "\n", significant), "")),
                         position = position_stack(vjust = 0.5),
                         size = 3)

            ggsave(file.path(output_dir, "distributions",
                            paste0("pie_chart_", comp, ".pdf")), p, width = 8, height = 8)
        }
    }

    # Save summary table
    write.csv(summary_df, file.path(output_dir, "distributions", "event_summary.csv"),
              row.names = FALSE)

    return(summary_df)
}

# =============================================================================
# Function: Plot dPSI Distribution (Directionality Analysis)
# =============================================================================
plot_dpsi_distribution <- function(splicing_dir, comparisons, event_types, output_dir) {
    cat("Generating dPSI distribution plots (directionality analysis)...\n")

    # Collect all significant events - only keep common columns
    all_sig <- list()

    for (comp in comparisons) {
        for (event in event_types) {
            df <- load_rmats_results(comp, event, splicing_dir)
            if (!is.null(df) && nrow(df) > 0) {
                sig <- df %>%
                    filter(FDR < FDR_THRESHOLD, abs(IncLevelDifference) > DPSI_THRESHOLD) %>%
                    select(ID, GeneID, geneSymbol, FDR, IncLevelDifference, event_type, comparison)
                if (nrow(sig) > 0) {
                    all_sig[[paste(comp, event, sep = "_")]] <- sig
                }
            }
        }
    }

    if (length(all_sig) == 0) {
        warning("No significant events found")
        return(NULL)
    }

    combined_sig <- do.call(rbind, all_sig)

    # 1. dPSI histogram for all comparisons
    p1 <- ggplot(combined_sig, aes(x = IncLevelDifference, fill = comparison)) +
        geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
        facet_wrap(~comparison, scales = "free_y") +
        theme_bw() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        labs(
            x = "dPSI (IncLevelDifference)",
            y = "Count",
            title = "Distribution of dPSI for Significant Splicing Events"
        ) +
        theme(legend.position = "none")

    ggsave(file.path(output_dir, "distributions", "dPSI_histogram_all.pdf"),
           p1, width = 14, height = 10)

    # 2. Density plot by event type
    p2 <- ggplot(combined_sig, aes(x = IncLevelDifference, color = event_type)) +
        geom_density(linewidth = 1) +
        facet_wrap(~comparison) +
        scale_color_manual(values = event_colors) +
        theme_bw() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        labs(
            x = "dPSI (IncLevelDifference)",
            y = "Density",
            title = "dPSI Distribution by Event Type",
            color = "Event Type"
        )

    ggsave(file.path(output_dir, "distributions", "dPSI_density_by_type.pdf"),
           p2, width = 14, height = 10)

    # 3. Directionality summary
    direction_summary <- combined_sig %>%
        mutate(direction = ifelse(IncLevelDifference > 0, "Inclusion", "Exclusion")) %>%
        group_by(comparison, event_type, direction) %>%
        summarize(count = n(), .groups = "drop")

    p3 <- ggplot(direction_summary, aes(x = event_type, y = count, fill = direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = c("Inclusion" = "#e74c3c", "Exclusion" = "#3498db")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Event Type",
            y = "Number of Events",
            title = "Splicing Directionality: Inclusion vs Exclusion",
            fill = "Direction"
        )

    ggsave(file.path(output_dir, "distributions", "directionality_barplot.pdf"),
           p3, width = 14, height = 10)

    # Save directionality summary
    direction_wide <- direction_summary %>%
        pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
        mutate(
            total = Inclusion + Exclusion,
            inclusion_ratio = Inclusion / total
        )

    write.csv(direction_wide, file.path(output_dir, "distributions", "directionality_summary.csv"),
              row.names = FALSE)

    return(direction_summary)
}

# =============================================================================
# Function: Cross-Comparison Plots
# =============================================================================
plot_comparison_crossplot <- function(splicing_dir, parental_comparisons, output_dir) {
    cat("Generating cross-comparison scatter plots...\n")

    # Focus on SE events for cross-plots
    event_type <- "SE"

    # Load data for parental comparisons
    dpsi_data <- list()
    for (comp in parental_comparisons) {
        df <- load_rmats_results(comp, event_type, splicing_dir)
        if (!is.null(df) && nrow(df) > 0) {
            # Create unique event ID
            df$event_id <- paste(df$GeneID, df$chr, df$exonStart_0base, df$exonEnd, sep = "_")
            dpsi_data[[comp]] <- df %>%
                select(event_id, GeneID, geneSymbol, IncLevelDifference, FDR)
        }
    }

    if (length(dpsi_data) < 2) {
        warning("Not enough data for cross-plots")
        return(NULL)
    }

    # Create pairwise comparisons
    comp_pairs <- combn(names(dpsi_data), 2, simplify = FALSE)

    for (pair in comp_pairs) {
        comp1 <- pair[1]
        comp2 <- pair[2]

        # Merge data
        merged <- merge(
            dpsi_data[[comp1]] %>% select(event_id, geneSymbol,
                                           dPSI_1 = IncLevelDifference, FDR_1 = FDR),
            dpsi_data[[comp2]] %>% select(event_id,
                                           dPSI_2 = IncLevelDifference, FDR_2 = FDR),
            by = "event_id"
        )

        if (nrow(merged) < 10) next

        # Classify events
        merged <- merged %>%
            mutate(
                sig_1 = FDR_1 < FDR_THRESHOLD & abs(dPSI_1) > DPSI_THRESHOLD,
                sig_2 = FDR_2 < FDR_THRESHOLD & abs(dPSI_2) > DPSI_THRESHOLD,
                category = case_when(
                    sig_1 & sig_2 ~ "Both",
                    sig_1 & !sig_2 ~ comp1,
                    !sig_1 & sig_2 ~ comp2,
                    TRUE ~ "Neither"
                )
            )

        # Calculate correlation (on all events for reference)
        cor_all <- cor(merged$dPSI_1, merged$dPSI_2, use = "complete.obs")

        # Filter to show only significant events (remove "Neither")
        merged_sig <- merged %>% filter(category != "Neither")

        # Calculate correlation on significant events only
        cor_sig <- if(nrow(merged_sig) > 2) {
            cor(merged_sig$dPSI_1, merged_sig$dPSI_2, use = "complete.obs")
        } else { NA }

        # Color palette (no "Neither" needed)
        cat_colors <- c("Both" = "#e74c3c")
        cat_colors[comp1] <- "#3498db"
        cat_colors[comp2] <- "#2ecc71"

        # Cross-plot with only significant events
        p <- ggplot(merged_sig, aes(x = dPSI_1, y = dPSI_2, color = category)) +
            geom_point(alpha = 0.7, size = 2) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
            geom_hline(yintercept = c(-DPSI_THRESHOLD, DPSI_THRESHOLD),
                      linetype = "dotted", color = "grey70") +
            geom_vline(xintercept = c(-DPSI_THRESHOLD, DPSI_THRESHOLD),
                      linetype = "dotted", color = "grey70") +
            scale_color_manual(values = cat_colors) +
            theme_bw() +
            labs(
                x = paste0("dPSI: ", comp1),
                y = paste0("dPSI: ", comp2),
                title = paste0("Cross-comparison: ", comp1, " vs ", comp2),
                subtitle = paste0("Correlation (sig): ", round(cor_sig, 3),
                                 " | Significant events: ", nrow(merged_sig),
                                 " / ", nrow(merged), " total"),
                color = "Significant in"
            ) +
            coord_fixed()

        # Add labels for highly significant in both
        top_both <- merged_sig %>%
            filter(category == "Both") %>%
            arrange(desc(abs(dPSI_1) + abs(dPSI_2))) %>%
            head(10)

        if (nrow(top_both) > 0) {
            p <- p + geom_text_repel(
                data = top_both,
                aes(label = geneSymbol),
                size = 2.5,
                max.overlaps = 15
            )
        }

        # Only save if we have significant events to plot
        if (nrow(merged_sig) > 0) {
            filename <- paste0("crossplot_", gsub("_vs_Parental", "", comp1),
                              "_vs_", gsub("_vs_Parental", "", comp2), ".pdf")
            ggsave(file.path(output_dir, "crossplots", filename), p, width = 10, height = 10)
        }

        # Save merged data
        csv_filename <- paste0("crossplot_data_", gsub("_vs_Parental", "", comp1),
                               "_vs_", gsub("_vs_Parental", "", comp2), ".csv")
        write.csv(merged, file.path(output_dir, "crossplots", csv_filename), row.names = FALSE)
    }
}

# =============================================================================
# Function: Find Common Events
# =============================================================================
find_common_events <- function(splicing_dir, parental_comparisons, event_types, output_dir) {
    cat("Finding common events between comparisons...\n")

    # Collect significant events for each comparison
    sig_events <- list()

    for (comp in parental_comparisons) {
        comp_events <- character()

        for (event in event_types) {
            df <- load_rmats_results(comp, event, splicing_dir)
            if (!is.null(df) && nrow(df) > 0) {
                sig <- df %>%
                    filter(FDR < FDR_THRESHOLD, abs(IncLevelDifference) > DPSI_THRESHOLD)

                if (nrow(sig) > 0) {
                    # Create unique event ID
                    event_ids <- paste(event, sig$GeneID, sig$chr,
                                      sig$exonStart_0base, sig$exonEnd, sep = "_")
                    comp_events <- c(comp_events, event_ids)
                }
            }
        }

        sig_events[[comp]] <- unique(comp_events)
        cat(paste0("  ", comp, ": ", length(sig_events[[comp]]), " significant events\n"))
    }

    # Calculate overlaps
    if (length(sig_events) >= 2) {
        # Pairwise overlaps
        for (i in 1:(length(sig_events) - 1)) {
            for (j in (i + 1):length(sig_events)) {
                comp1 <- names(sig_events)[i]
                comp2 <- names(sig_events)[j]
                overlap <- length(intersect(sig_events[[comp1]], sig_events[[comp2]]))
                cat(paste0("  Overlap ", comp1, " & ", comp2, ": ", overlap, " events\n"))
            }
        }

        # All three overlaps
        if (length(sig_events) >= 3) {
            common_all <- Reduce(intersect, sig_events[1:3])
            cat(paste0("  Common to all three: ", length(common_all), " events\n"))
        }
    }

    # Create Venn diagram if we have 2-3 comparisons
    if (length(parental_comparisons) >= 2 && length(parental_comparisons) <= 3) {
        venn_list <- sig_events[parental_comparisons]
        names(venn_list) <- gsub("_vs_Parental", "", names(venn_list))

        # Generate Venn diagram
        pdf(file.path(output_dir, "common_events", "venn_diagram.pdf"), width = 10, height = 10)

        if (length(venn_list) == 2) {
            venn_plot <- draw.pairwise.venn(
                area1 = length(venn_list[[1]]),
                area2 = length(venn_list[[2]]),
                cross.area = length(intersect(venn_list[[1]], venn_list[[2]])),
                category = names(venn_list),
                fill = c("#3498db", "#2ecc71"),
                alpha = 0.5,
                cat.pos = c(-20, 20)
            )
        } else {
            venn_plot <- draw.triple.venn(
                area1 = length(venn_list[[1]]),
                area2 = length(venn_list[[2]]),
                area3 = length(venn_list[[3]]),
                n12 = length(intersect(venn_list[[1]], venn_list[[2]])),
                n23 = length(intersect(venn_list[[2]], venn_list[[3]])),
                n13 = length(intersect(venn_list[[1]], venn_list[[3]])),
                n123 = length(Reduce(intersect, venn_list)),
                category = names(venn_list),
                fill = c("#3498db", "#2ecc71", "#e74c3c"),
                alpha = 0.5,
                cat.pos = c(-20, 20, 180)
            )
        }

        grid.draw(venn_plot)
        dev.off()
    }

    # Save common events
    if (length(sig_events) >= 3) {
        common_all <- Reduce(intersect, sig_events[parental_comparisons])

        if (length(common_all) > 0) {
            # Get details for common events
            common_details <- list()

            for (comp in parental_comparisons) {
                for (event in event_types) {
                    df <- load_rmats_results(comp, event, splicing_dir)
                    if (!is.null(df) && nrow(df) > 0) {
                        df$event_id <- paste(event, df$GeneID, df$chr,
                                            df$exonStart_0base, df$exonEnd, sep = "_")
                        common_df <- df %>%
                            filter(event_id %in% common_all) %>%
                            select(event_id, GeneID, geneSymbol, event_type,
                                   IncLevelDifference, FDR) %>%
                            rename(!!paste0("dPSI_", gsub("_vs_Parental", "", comp)) := IncLevelDifference,
                                   !!paste0("FDR_", gsub("_vs_Parental", "", comp)) := FDR)

                        if (nrow(common_df) > 0) {
                            common_details[[paste(comp, event, sep = "_")]] <- common_df
                        }
                    }
                }
            }

            if (length(common_details) > 0) {
                # Merge all details
                merged_common <- Reduce(function(x, y) {
                    merge(x, y, by = c("event_id", "GeneID", "geneSymbol", "event_type"),
                          all = TRUE)
                }, common_details)

                write.csv(merged_common,
                         file.path(output_dir, "common_events", "common_events_details.csv"),
                         row.names = FALSE)
            }
        }
    }

    # Save overlap statistics
    overlap_stats <- data.frame(
        comparison = names(sig_events),
        n_significant = sapply(sig_events, length)
    )
    write.csv(overlap_stats, file.path(output_dir, "common_events", "overlap_statistics.csv"),
              row.names = FALSE)

    return(sig_events)
}

# =============================================================================
# Main Execution
# =============================================================================

# 1. Extract PSI matrix and run PCA
cat("\n=== Step 1: PCA on Splicing Data ===\n")
psi_matrix <- extract_psi_matrix(SPLICING_DIR, event_types)
if (!is.null(psi_matrix)) {
    pca_result <- plot_splicing_pca(psi_matrix, sample_info, OUTPUT_DIR)

    # Save PSI matrix
    write.csv(psi_matrix, file.path(OUTPUT_DIR, "PCA", "psi_matrix.csv"))
}

# 2. Event distribution analysis
cat("\n=== Step 2: Event Distribution Analysis ===\n")
event_summary <- plot_event_distribution(SPLICING_DIR, comparisons, event_types, OUTPUT_DIR)

# 3. dPSI distribution and directionality
cat("\n=== Step 3: dPSI Distribution and Directionality ===\n")
direction_summary <- plot_dpsi_distribution(SPLICING_DIR, comparisons, event_types, OUTPUT_DIR)

# 4. Cross-comparison plots
cat("\n=== Step 4: Cross-Comparison Plots ===\n")
plot_comparison_crossplot(SPLICING_DIR, parental_comparisons, OUTPUT_DIR)

# 5. Common events analysis
cat("\n=== Step 5: Common Events Analysis ===\n")
common_events <- find_common_events(SPLICING_DIR, parental_comparisons, event_types, OUTPUT_DIR)

# =============================================================================
# Summary
# =============================================================================

cat("\n============================================\n")
cat("Splicing Downstream Analysis Complete!\n")
cat("============================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - PCA/: PCA plots and coordinates\n")
cat("  - distributions/: Event type distributions, dPSI histograms, directionality\n")
cat("  - crossplots/: Comparison scatter plots\n")
cat("  - common_events/: Venn diagrams, overlap statistics\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
