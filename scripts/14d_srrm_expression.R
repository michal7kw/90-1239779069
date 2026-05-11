#!/usr/bin/env Rscript
# =============================================================================
# 14d_srrm_expression.R
# Visualize Srrm3 and Srrm4 expression across conditions with statistical
# annotations from DESeq2.
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
DESEQ_DIR <- file.path(BASE_DIR, "results/04_deseq2")
OUTPUT_DIR <- file.path(BASE_DIR, "results/14_todo_analysis/task3_srrm_expression")

group_colors <- c(Parental = "#E41A1C", Neg = "#377EB8", Pos = "#4DAF4A", KO = "#984EA3")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load and reshape expression data
# =============================================================================
cat("========================================\n")
cat("14d: Srrm3 and Srrm4 Expression\n")
cat("========================================\n")

counts <- read.csv(file.path(DESEQ_DIR, "normalized_counts.csv"),
                   stringsAsFactors = FALSE)

srrm_genes <- c("Srrm3", "Srrm4")
srrm_data <- counts %>% filter(gene_symbol %in% srrm_genes)

if (nrow(srrm_data) == 0) {
    stop("Srrm3/Srrm4 not found in normalized_counts.csv")
}
cat(sprintf("  Found %d SRRM genes\n", nrow(srrm_data)))

sample_cols <- colnames(counts)[3:ncol(counts)]
srrm_long <- srrm_data %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "expression") %>%
    mutate(
        group = gsub("_[0-9]+$", "", sample),
        group = factor(group, levels = c("Parental", "Neg", "Pos", "KO"))
    )

write.csv(srrm_long, file.path(OUTPUT_DIR, "srrm_expression_data.csv"), row.names = FALSE)

# Group summary
srrm_summary <- srrm_long %>%
    group_by(gene_symbol, group) %>%
    summarize(
        mean = mean(expression),
        sd = sd(expression),
        sem = sd(expression) / sqrt(n()),
        n = n(),
        .groups = "drop"
    )
write.csv(srrm_summary, file.path(OUTPUT_DIR, "srrm_group_summary.csv"), row.names = FALSE)

cat("\n  Expression summary:\n")
print(as.data.frame(srrm_summary))

# =============================================================================
# Load DESeq2 significance
# =============================================================================
comparisons_vs_parental <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")
deseq_stats <- list()
for (comp in comparisons_vs_parental) {
    f <- file.path(DESEQ_DIR, paste0(comp, "_all.csv"))
    if (file.exists(f)) {
        df <- read.csv(f, stringsAsFactors = FALSE)
        for (g in srrm_genes) {
            row <- df %>% filter(gene_symbol == g)
            if (nrow(row) > 0) {
                deseq_stats[[paste(g, comp, sep = "_")]] <- data.frame(
                    gene_symbol = g,
                    comparison = comp,
                    log2FC = row$log2FoldChange[1],
                    padj = row$padj[1],
                    stringsAsFactors = FALSE
                )
            }
        }
    }
}
deseq_df <- bind_rows(deseq_stats)
deseq_df$sig_label <- ifelse(deseq_df$padj < 0.001, "***",
                        ifelse(deseq_df$padj < 0.01, "**",
                          ifelse(deseq_df$padj < 0.05, "*", "ns")))
write.csv(deseq_df, file.path(OUTPUT_DIR, "srrm_deseq2_stats.csv"), row.names = FALSE)

cat("\n  DESeq2 significance:\n")
print(deseq_df %>% select(gene_symbol, comparison, log2FC, padj, sig_label))

# =============================================================================
# Bar plot with individual points + SEM + significance brackets
# =============================================================================
sig_annot <- deseq_df %>%
    mutate(
        group = gsub("_vs_Parental", "", comparison),
        group = factor(group, levels = c("Neg", "Pos", "KO"))
    )

ymax_per_gene <- srrm_long %>%
    group_by(gene_symbol) %>%
    summarize(ymax = max(expression) * 1.05, .groups = "drop")

p_bar <- ggplot(srrm_summary, aes(x = group, y = mean, fill = group)) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    geom_jitter(data = srrm_long, aes(x = group, y = expression, fill = group),
                width = 0.15, size = 2, shape = 21, color = "black", alpha = 0.7) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = "Srrm3 and Srrm4 Normalized Expression",
         subtitle = "DESeq2 normalized counts, mean +/- SEM",
         x = "Condition", y = "Normalized Expression") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold.italic", size = 13))

# Add significance brackets
for (g in srrm_genes) {
    gene_sig <- sig_annot %>% filter(gene_symbol == g, sig_label != "ns")
    if (nrow(gene_sig) > 0) {
        ymax_g <- ymax_per_gene$ymax[ymax_per_gene$gene_symbol == g]
        for (i in seq_len(nrow(gene_sig))) {
            bracket_y <- ymax_g + (ymax_g * 0.08 * i)
            # sig_annot$group has levels Neg=1,Pos=2,KO=3 but plot x-axis
            # has Parental=1,Neg=2,Pos=3,KO=4, so offset by +1
            xend_pos <- as.numeric(gene_sig$group[i]) + 1
            p_bar <- p_bar +
                geom_segment(data = data.frame(
                    gene_symbol = g,
                    x = 1, xend = xend_pos,
                    y = bracket_y, yend = bracket_y
                ), aes(x = x, xend = xend, y = y, yend = yend),
                inherit.aes = FALSE, linewidth = 0.4) +
                geom_text(data = data.frame(
                    gene_symbol = g,
                    x = (1 + xend_pos) / 2,
                    y = bracket_y + ymax_g * 0.02,
                    label = gene_sig$sig_label[i]
                ), aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 4, vjust = 0)
        }
    }
}

ggsave(file.path(OUTPUT_DIR, "srrm_barplot.pdf"), p_bar, width = 9, height = 6)

# =============================================================================
# Box plot variant
# =============================================================================
p_box <- ggplot(srrm_long, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, shape = 21, color = "black", alpha = 0.8) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = "Srrm3 and Srrm4 Expression — Box Plot",
         x = "Condition", y = "Normalized Expression") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold.italic", size = 13))
ggsave(file.path(OUTPUT_DIR, "srrm_boxplot.pdf"), p_box, width = 9, height = 6)

# =============================================================================
# Fold-change comparison
# =============================================================================
fc_data <- srrm_summary %>%
    group_by(gene_symbol) %>%
    mutate(
        parental_mean = mean[group == "Parental"],
        fold_change = mean / parental_mean,
        fc_sem = sem / parental_mean
    ) %>%
    ungroup() %>%
    filter(group != "Parental")

p_fc <- ggplot(fc_data, aes(x = group, y = fold_change, fill = gene_symbol)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = fold_change - fc_sem, ymax = fold_change + fc_sem),
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("Srrm3" = "#E41A1C", "Srrm4" = "#377EB8"),
                      name = "Gene") +
    labs(title = "Relative Expression Change vs Parental",
         subtitle = "Fold change of mean normalized expression",
         x = "Condition", y = "Fold Change (vs Parental)") +
    theme_minimal(base_size = 13) +
    theme(legend.text = element_text(face = "italic"))
ggsave(file.path(OUTPUT_DIR, "srrm_foldchange_comparison.pdf"), p_fc, width = 7, height = 5)

cat("\n14d complete.\n")
