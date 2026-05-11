#!/usr/bin/env Rscript
# =============================================================================
# Microexon Expression Correlation Analysis
# Project: 90-1239779069
# Analysis: Correlate dPSI (splicing changes) with log2FC (expression changes)
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
    library(tibble)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
DESEQ_DIR <- file.path(BASE_DIR, "results/04_deseq2")
OUTPUT_DIR <- file.path(BASE_DIR, "results/09_microexon_extended/expression_correlation")

# Thresholds
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1

# Exon size bins
MICROEXON_MAX <- 30
SMALL_EXON_MAX <- 50

cat("============================================\n")
cat("Microexon Expression Correlation Analysis\n")
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
size_colors <- c("Microexon (0-30bp)" = "#e74c3c",
                 "Small (30-50bp)" = "#f39c12",
                 "Regular (>50bp)" = "#3498db")

# =============================================================================
# Function: Load SE events with size classification
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

    # Mark significant
    df$sig_splicing <- df$FDR < FDR_THRESHOLD & abs(df$IncLevelDifference) > DPSI_THRESHOLD
    df$dPSI <- df$IncLevelDifference
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
# Function: Load DESeq2 results
# =============================================================================
load_deseq_data <- function(comparison, deseq_dir) {
    file_path <- file.path(deseq_dir, paste0(comparison, "_all.csv"))

    if (!file.exists(file_path)) {
        warning(paste("File not found:", file_path))
        return(NULL)
    }

    df <- read.csv(file_path, stringsAsFactors = FALSE)

    # Clean gene ID column name
    if ("gene_id" %in% names(df)) {
        df$GeneID <- df$gene_id
    } else if ("X" %in% names(df)) {
        df$GeneID <- df$X
    }

    # Mark significant DE genes
    df$sig_expression <- !is.na(df$padj) & df$padj < PADJ_THRESHOLD & abs(df$log2FoldChange) > LFC_THRESHOLD

    return(df)
}

# =============================================================================
# Function: Merge splicing and expression data
# =============================================================================
merge_splicing_expression <- function(splicing_df, deseq_df) {
    cat("Merging splicing and expression data...\n")

    # Clean gene IDs (remove quotes and version numbers for matching)
    splicing_df$clean_gene_id <- gsub("\"", "", splicing_df$GeneID)
    splicing_df$clean_gene_id <- gsub("\\..*", "", splicing_df$clean_gene_id)

    deseq_df$clean_gene_id <- gsub("\"", "", deseq_df$GeneID)
    deseq_df$clean_gene_id <- gsub("\\..*", "", deseq_df$clean_gene_id)

    # Merge by gene ID
    merged <- merge(
        splicing_df %>% select(GeneID, geneSymbol, chr, exonStart_0base, exonEnd,
                               exon_length, size_category, dPSI, FDR, sig_splicing,
                               comparison, clean_gene_id),
        deseq_df %>% select(clean_gene_id, log2FoldChange, padj, sig_expression),
        by = "clean_gene_id",
        all.x = TRUE
    )

    cat(paste0("  Matched ", sum(!is.na(merged$log2FoldChange)), " events with expression data\n"))
    cat(paste0("  Total events: ", nrow(merged), "\n"))

    return(merged)
}

# =============================================================================
# Function: Calculate correlations by size category
# =============================================================================
calculate_correlation_by_size <- function(merged_df) {
    cat("Calculating correlations by exon size...\n")

    correlation_results <- merged_df %>%
        filter(!is.na(log2FoldChange), !is.na(dPSI)) %>%
        group_by(comparison, size_category) %>%
        summarize(
            n_events = n(),
            n_sig_both = sum(sig_splicing & sig_expression, na.rm = TRUE),
            n_sig_splicing_only = sum(sig_splicing & !sig_expression, na.rm = TRUE),
            n_sig_expression_only = sum(!sig_splicing & sig_expression, na.rm = TRUE),
            n_neither_sig = sum(!sig_splicing & !sig_expression, na.rm = TRUE),
            pearson_r = cor(dPSI, log2FoldChange, use = "complete.obs", method = "pearson"),
            spearman_rho = cor(dPSI, log2FoldChange, use = "complete.obs", method = "spearman"),
            .groups = "drop"
        )

    # Add p-values for correlations
    for (i in seq_len(nrow(correlation_results))) {
        comp <- correlation_results$comparison[i]
        size <- correlation_results$size_category[i]

        subset_data <- merged_df %>%
            filter(comparison == comp, size_category == size,
                   !is.na(log2FoldChange), !is.na(dPSI))

        if (nrow(subset_data) >= 3) {
            pearson_test <- cor.test(subset_data$dPSI, subset_data$log2FoldChange, method = "pearson")
            spearman_test <- cor.test(subset_data$dPSI, subset_data$log2FoldChange, method = "spearman")

            correlation_results$pearson_p[i] <- pearson_test$p.value
            correlation_results$spearman_p[i] <- spearman_test$p.value
        } else {
            correlation_results$pearson_p[i] <- NA
            correlation_results$spearman_p[i] <- NA
        }
    }

    cat("\nCorrelation summary:\n")
    print(correlation_results)

    return(correlation_results)
}

# =============================================================================
# Function: Plot splicing-expression scatter
# =============================================================================
plot_splicing_expression_scatter <- function(merged_df, comparison_name, output_dir) {
    cat(paste0("Generating scatter plots for ", comparison_name, "...\n"))

    plot_data <- merged_df %>%
        filter(!is.na(log2FoldChange), !is.na(dPSI))

    if (nrow(plot_data) == 0) {
        warning(paste("No data for scatter plot:", comparison_name))
        return(NULL)
    }

    # Add significance category
    plot_data <- plot_data %>%
        mutate(sig_category = case_when(
            sig_splicing & sig_expression ~ "Both significant",
            sig_splicing & !sig_expression ~ "Splicing only",
            !sig_splicing & sig_expression ~ "Expression only",
            TRUE ~ "Neither"
        ))

    # 1. Basic scatter by size category
    p1 <- ggplot(plot_data, aes(x = log2FoldChange, y = dPSI)) +
        geom_point(aes(color = size_category), alpha = 0.5, size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
        scale_color_manual(values = size_colors) +
        theme_bw() +
        labs(
            x = "log2 Fold Change (Expression)",
            y = "dPSI (Splicing)",
            title = paste0("Expression vs Splicing: ", comparison_name),
            color = "Exon Size"
        )

    ggsave(file.path(output_dir, paste0(comparison_name, "_scatter_by_size.pdf")),
           p1, width = 10, height = 8)

    # 2. Faceted by size category with regression lines
    p2 <- ggplot(plot_data, aes(x = log2FoldChange, y = dPSI)) +
        geom_point(alpha = 0.4, size = 1, color = "#3498db") +
        geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
        facet_wrap(~size_category, scales = "free") +
        theme_bw() +
        labs(
            x = "log2 Fold Change (Expression)",
            y = "dPSI (Splicing)",
            title = paste0("Expression-Splicing Correlation by Exon Size: ", comparison_name)
        )

    ggsave(file.path(output_dir, paste0(comparison_name, "_scatter_faceted.pdf")),
           p2, width = 14, height = 5)

    # 3. Significant events highlighted
    p3 <- ggplot(plot_data, aes(x = log2FoldChange, y = dPSI)) +
        geom_point(aes(color = sig_category), alpha = 0.6, size = 1.5) +
        geom_hline(yintercept = c(-DPSI_THRESHOLD, DPSI_THRESHOLD),
                   linetype = "dashed", color = "gray60") +
        geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
                   linetype = "dashed", color = "gray60") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
        scale_color_manual(values = c("Both significant" = "#e74c3c",
                                      "Splicing only" = "#3498db",
                                      "Expression only" = "#2ecc71",
                                      "Neither" = "#bdc3c7")) +
        theme_bw() +
        labs(
            x = "log2 Fold Change (Expression)",
            y = "dPSI (Splicing)",
            title = paste0("Significant Events: ", comparison_name),
            subtitle = paste0("Thresholds: |dPSI| > ", DPSI_THRESHOLD, ", |log2FC| > ", LFC_THRESHOLD),
            color = "Significance"
        )

    ggsave(file.path(output_dir, paste0(comparison_name, "_scatter_significance.pdf")),
           p3, width = 10, height = 8)

    # 4. Label top genes (significant in both with extreme values)
    top_genes <- plot_data %>%
        filter(sig_category == "Both significant") %>%
        mutate(extremity = abs(dPSI) + abs(log2FoldChange)) %>%
        arrange(desc(extremity)) %>%
        head(20)

    if (nrow(top_genes) > 0) {
        p4 <- ggplot(plot_data, aes(x = log2FoldChange, y = dPSI)) +
            geom_point(aes(color = sig_category), alpha = 0.4, size = 1) +
            geom_point(data = top_genes, color = "red", size = 2) +
            geom_text_repel(data = top_genes, aes(label = geneSymbol),
                           size = 3, max.overlaps = 15) +
            scale_color_manual(values = c("Both significant" = "#e74c3c",
                                          "Splicing only" = "#3498db",
                                          "Expression only" = "#2ecc71",
                                          "Neither" = "#bdc3c7")) +
            theme_bw() +
            labs(
                x = "log2 Fold Change (Expression)",
                y = "dPSI (Splicing)",
                title = paste0("Top Genes: ", comparison_name),
                color = "Significance"
            )

        ggsave(file.path(output_dir, paste0(comparison_name, "_scatter_labeled.pdf")),
               p4, width = 12, height = 10)
    }

    return(p1)
}

# =============================================================================
# Function: Generate correlation heatmap
# =============================================================================
generate_correlation_heatmap <- function(all_correlations, output_dir) {
    cat("Generating correlation heatmap...\n")

    # Reshape for heatmap
    cor_matrix <- all_correlations %>%
        select(comparison, size_category, spearman_rho) %>%
        pivot_wider(names_from = size_category, values_from = spearman_rho) %>%
        column_to_rownames("comparison")

    # Create heatmap
    if (nrow(cor_matrix) > 0 && ncol(cor_matrix) > 0) {
        pdf(file.path(output_dir, "correlation_heatmap.pdf"), width = 8, height = 6)
        pheatmap(
            as.matrix(cor_matrix),
            main = "Spearman Correlation (dPSI vs log2FC)",
            color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(50),
            display_numbers = TRUE,
            number_format = "%.3f",
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            fontsize = 12
        )
        dev.off()
    }

    # Also create bar plot version
    p <- ggplot(all_correlations, aes(x = comparison, y = spearman_rho, fill = size_category)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = size_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Spearman Correlation (rho)",
            title = "Expression-Splicing Correlation by Exon Size",
            fill = "Exon Size"
        )

    ggsave(file.path(output_dir, "correlation_barplot.pdf"), p, width = 10, height = 6)
}

# =============================================================================
# Function: Analyze quadrant distribution
# =============================================================================
analyze_quadrants <- function(merged_df, output_dir) {
    cat("Analyzing quadrant distribution...\n")

    quadrant_data <- merged_df %>%
        filter(!is.na(log2FoldChange), !is.na(dPSI), sig_splicing) %>%
        mutate(quadrant = case_when(
            dPSI > 0 & log2FoldChange > 0 ~ "Inc+Up",
            dPSI > 0 & log2FoldChange < 0 ~ "Inc+Down",
            dPSI < 0 & log2FoldChange > 0 ~ "Exc+Up",
            dPSI < 0 & log2FoldChange < 0 ~ "Exc+Down",
            TRUE ~ "Neutral"
        ))

    quadrant_summary <- quadrant_data %>%
        group_by(comparison, size_category, quadrant) %>%
        summarize(count = n(), .groups = "drop") %>%
        pivot_wider(names_from = quadrant, values_from = count, values_fill = 0)

    write.csv(quadrant_summary, file.path(output_dir, "quadrant_distribution.csv"), row.names = FALSE)

    # Plot quadrant distribution
    quadrant_plot_data <- quadrant_data %>%
        group_by(comparison, size_category, quadrant) %>%
        summarize(count = n(), .groups = "drop")

    p <- ggplot(quadrant_plot_data, aes(x = quadrant, y = count, fill = size_category)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = size_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Quadrant (Splicing + Expression Direction)",
            y = "Number of Events",
            title = "Splicing-Expression Quadrant Analysis",
            subtitle = "Inc = Inclusion, Exc = Exclusion, Up = Upregulated, Down = Downregulated",
            fill = "Exon Size"
        )

    ggsave(file.path(output_dir, "quadrant_distribution.pdf"), p, width = 14, height = 10)

    return(quadrant_summary)
}

# =============================================================================
# Main Execution
# =============================================================================

all_merged_data <- list()
all_correlations <- list()

for (comp in parental_comparisons) {
    cat(paste0("\n=== Processing: ", comp, " ===\n"))

    # 1. Load data
    splicing_data <- load_se_data(comp, SPLICING_DIR)
    deseq_data <- load_deseq_data(comp, DESEQ_DIR)

    if (is.null(splicing_data) || is.null(deseq_data)) {
        warning(paste("Skipping", comp, "due to missing data"))
        next
    }

    # 2. Merge data
    merged_data <- merge_splicing_expression(splicing_data, deseq_data)
    merged_data$comparison <- comp
    all_merged_data[[comp]] <- merged_data

    # 3. Calculate correlations
    correlations <- calculate_correlation_by_size(merged_data)
    all_correlations[[comp]] <- correlations

    # 4. Generate scatter plots
    plot_splicing_expression_scatter(merged_data, comp, OUTPUT_DIR)

    # 5. Save correlation data
    write.csv(correlations, file.path(OUTPUT_DIR, paste0(comp, "_correlation.csv")), row.names = FALSE)

    # 6. Save merged data for significant events
    sig_merged <- merged_data %>%
        filter(sig_splicing | sig_expression) %>%
        select(GeneID, geneSymbol, chr, exonStart_0base, exonEnd, exon_length,
               size_category, dPSI, FDR, sig_splicing, log2FoldChange, padj, sig_expression)

    write.csv(sig_merged, file.path(OUTPUT_DIR, paste0(comp, "_merged_significant.csv")), row.names = FALSE)
}

# Combine all correlations
cat("\n=== Combined Analysis ===\n")
combined_correlations <- do.call(rbind, all_correlations)
write.csv(combined_correlations, file.path(OUTPUT_DIR, "correlation_summary.csv"), row.names = FALSE)

# Generate combined heatmap
generate_correlation_heatmap(combined_correlations, OUTPUT_DIR)

# Combine all merged data and analyze quadrants
combined_merged <- do.call(rbind, all_merged_data)
quadrant_results <- analyze_quadrants(combined_merged, OUTPUT_DIR)

# =============================================================================
# Combined scatter plot (all comparisons)
# =============================================================================
cat("Generating combined scatter plot...\n")

plot_data <- combined_merged %>%
    filter(!is.na(log2FoldChange), !is.na(dPSI), sig_splicing)

p_combined <- ggplot(plot_data, aes(x = log2FoldChange, y = dPSI)) +
    geom_point(aes(color = comparison), alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    facet_wrap(~size_category) +
    scale_color_manual(values = comparison_colors) +
    theme_bw() +
    labs(
        x = "log2 Fold Change (Expression)",
        y = "dPSI (Splicing)",
        title = "Expression vs Splicing: All Comparisons",
        subtitle = "Significant splicing events only",
        color = "Comparison"
    )

ggsave(file.path(OUTPUT_DIR, "combined_scatter_all_comparisons.pdf"),
       p_combined, width = 14, height = 6)

# =============================================================================
# Summary
# =============================================================================

cat("\n============================================\n")
cat("Expression Correlation Analysis Summary\n")
cat("============================================\n")

for (comp in parental_comparisons) {
    if (comp %in% names(all_merged_data)) {
        merged <- all_merged_data[[comp]]
        corr <- all_correlations[[comp]]

        cat(paste0("\n", comp, ":\n"))
        cat(paste0("  Total events merged: ", nrow(merged), "\n"))
        cat(paste0("  Events with expression data: ", sum(!is.na(merged$log2FoldChange)), "\n"))
        cat(paste0("  Significant in both: ", sum(merged$sig_splicing & merged$sig_expression, na.rm = TRUE), "\n"))

        # Print correlations for microexons
        micro_corr <- corr %>% filter(size_category == "Microexon (0-30bp)")
        if (nrow(micro_corr) > 0) {
            cat(paste0("  Microexon Spearman rho: ", round(micro_corr$spearman_rho[1], 3),
                      " (p=", format(micro_corr$spearman_p[1], digits = 3), ")\n"))
        }
    }
}

cat("\n============================================\n")
cat("Expression Correlation Analysis Complete!\n")
cat("============================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - {Comparison}_correlation.csv: Correlation statistics\n")
cat("  - {Comparison}_merged_significant.csv: Merged significant events\n")
cat("  - {Comparison}_scatter_*.pdf: Scatter plots\n")
cat("  - correlation_summary.csv: Combined correlation statistics\n")
cat("  - correlation_heatmap.pdf: Correlation heatmap\n")
cat("  - quadrant_distribution.csv/pdf: Quadrant analysis\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
