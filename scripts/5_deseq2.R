#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis
# Project: 90-1239779069
# Reference: Parental
# Comparisons: All pairwise (Neg vs Parental, Pos vs Parental, KO vs Parental,
#              Pos vs Neg, KO vs Neg, KO vs Pos)
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(tidyr)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/Azenta_projects/90-1239779069"
COUNT_DIR <- file.path(BASE_DIR, "results/02_aligned")
OUTPUT_DIR <- file.path(BASE_DIR, "results/04_deseq2")

cat("============================================\n")
cat("DESeq2 Differential Expression Analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

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
rownames(sample_info) <- sample_info$sample_name

cat("Sample information:\n")
print(sample_info)
cat("\n")

# Load count data from STAR ReadsPerGene.out.tab files
# STAR outputs 4 columns: gene_id, unstranded, forward, reverse
# We use column 2 (unstranded) or column 4 (reverse) depending on library prep
# For standard Illumina TruSeq, column 2 is usually appropriate

load_star_counts <- function(sample_id, count_dir) {
    file_path <- file.path(count_dir, sample_id,
                           paste0(sample_id, "_ReadsPerGene.out.tab"))
    if (!file.exists(file_path)) {
        stop(paste("Count file not found:", file_path))
    }
    counts <- read.table(file_path, header = FALSE, row.names = 1,
                         stringsAsFactors = FALSE)
    colnames(counts) <- c("unstranded", "forward", "reverse")
    # Remove first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
    counts <- counts[!grepl("^N_", rownames(counts)), ]
    return(counts$unstranded)  # Use unstranded counts
}

cat("Loading count data...\n")
count_list <- lapply(sample_info$sample_id, function(id) {
    load_star_counts(id, COUNT_DIR)
})

# Get gene names from first sample
first_file <- file.path(COUNT_DIR, sample_info$sample_id[1],
                        paste0(sample_info$sample_id[1], "_ReadsPerGene.out.tab"))
gene_info <- read.table(first_file, header = FALSE, row.names = 1,
                        stringsAsFactors = FALSE)
gene_names <- rownames(gene_info)[!grepl("^N_", rownames(gene_info))]

# Create count matrix
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- sample_info$sample_name

cat("Count matrix dimensions:", nrow(count_matrix), "genes x",
    ncol(count_matrix), "samples\n\n")

# Create gene ID to symbol mapping
cat("Mapping Ensembl IDs to gene symbols...\n")
# Remove version numbers from Ensembl IDs if present
gene_ids_clean <- gsub("\\..*", "", gene_names)
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = gene_ids_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
# Create mapping dataframe
gene_mapping <- data.frame(
    gene_id = gene_names,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE
)
# Use gene ID if symbol not found
gene_mapping$gene_symbol[is.na(gene_mapping$gene_symbol)] <-
    gene_mapping$gene_id[is.na(gene_mapping$gene_symbol)]
rownames(gene_mapping) <- gene_mapping$gene_id
cat("Gene symbols mapped:", sum(!is.na(gene_symbols)), "of", length(gene_names), "\n\n")

# Filter low-count genes (at least 10 counts in at least 3 samples)
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix_filtered <- count_matrix[keep, ]
cat("After filtering:", nrow(count_matrix_filtered), "genes retained\n\n")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                               colData = sample_info,
                               design = ~ group)

# Set Parental as reference level
dds$group <- relevel(dds$group, ref = "Parental")

# Run DESeq2
cat("Running DESeq2...\n")
dds <- DESeq(dds)
cat("DESeq2 analysis complete.\n\n")

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file.path(OUTPUT_DIR, "normalized_counts.csv"))

# VST transformation for visualization
vst_data <- vst(dds, blind = FALSE)

# PCA plot
cat("Generating PCA plot...\n")
pca_data <- plotPCA(vst_data, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.5, size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_bw() +
    theme(legend.position = "right") +
    scale_color_brewer(palette = "Set1") +
    ggtitle("PCA of samples")

ggsave(file.path(OUTPUT_DIR, "PCA_plot.pdf"), pca_plot, width = 10, height = 8)
ggsave(file.path(OUTPUT_DIR, "PCA_plot.png"), pca_plot, width = 10, height = 8, dpi = 300)

# Sample distance heatmap
cat("Generating sample distance heatmap...\n")
sample_dists <- dist(t(assay(vst_data)))
sample_dist_matrix <- as.matrix(sample_dists)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
annotation_col <- data.frame(Group = sample_info$group)
rownames(annotation_col) <- sample_info$sample_name

pdf(file.path(OUTPUT_DIR, "sample_distance_heatmap.pdf"), width = 10, height = 8)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         annotation_col = annotation_col,
         main = "Sample Distance Matrix")
dev.off()

# Define all pairwise comparisons
comparisons <- list(
    c("Neg", "Parental"),
    c("Pos", "Parental"),
    c("KO", "Parental"),
    c("Pos", "Neg"),
    c("KO", "Neg"),
    c("KO", "Pos")
)

# Function to extract and save results
extract_results <- function(dds, contrast, output_dir) {
    comparison_name <- paste0(contrast[1], "_vs_", contrast[2])
    cat("Processing:", comparison_name, "\n")

    # Get results
    res <- results(dds, contrast = c("group", contrast[1], contrast[2]))
    res <- res[order(res$padj), ]

    # Summary
    cat("  Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
    cat("  Up-regulated:", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE), "\n")
    cat("  Down-regulated:", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE), "\n\n")

    # Save full results
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                          "stat", "pvalue", "padj")]
    write.csv(res_df, file.path(output_dir, paste0(comparison_name, "_all.csv")),
              row.names = FALSE)

    # Save significant results
    sig_res <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
    write.csv(sig_res, file.path(output_dir, paste0(comparison_name, "_significant.csv")),
              row.names = FALSE)

    # MA plot
    pdf(file.path(output_dir, paste0(comparison_name, "_MA_plot.pdf")), width = 8, height = 6)
    plotMA(res, main = paste("MA Plot:", comparison_name), ylim = c(-5, 5))
    dev.off()

    # Volcano plot
    res_volcano <- as.data.frame(res)
    res_volcano$gene_id <- rownames(res_volcano)
    # Add gene symbols from mapping
    res_volcano$gene_symbol <- gene_mapping[res_volcano$gene_id, "gene_symbol"]
    res_volcano$significant <- ifelse(!is.na(res_volcano$padj) & res_volcano$padj < 0.05,
                                       ifelse(res_volcano$log2FoldChange > 0, "Up", "Down"),
                                       "NS")

    # Identify genes to label: top by significance and top by fold change
    sig_genes <- res_volcano[!is.na(res_volcano$padj) & res_volcano$padj < 0.05, ]

    # Top 10 most significant genes
    top_sig <- sig_genes[order(sig_genes$padj), ]
    top_sig <- head(top_sig, 10)

    # Top 10 by absolute log2FC (among significant)
    top_fc <- sig_genes[order(abs(sig_genes$log2FoldChange), decreasing = TRUE), ]
    top_fc <- head(top_fc, 10)

    # Combine and remove duplicates
    genes_to_label <- unique(c(rownames(top_sig), rownames(top_fc)))
    res_volcano$label <- ifelse(rownames(res_volcano) %in% genes_to_label,
                                 res_volcano$gene_symbol, "")

    volcano_plot <- ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(pvalue),
                                             color = significant)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_text_repel(aes(label = label),
                        size = 2.5,
                        max.overlaps = 20,
                        box.padding = 0.3,
                        point.padding = 0.2,
                        segment.color = "grey50",
                        segment.size = 0.2,
                        show.legend = FALSE) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
        theme_bw() +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
        xlab("log2 Fold Change") +
        ylab("-log10(p-value)") +
        ggtitle(paste("Volcano Plot:", comparison_name)) +
        theme(legend.position = "right")

    ggsave(file.path(output_dir, paste0(comparison_name, "_volcano.pdf")),
           volcano_plot, width = 10, height = 8)
    ggsave(file.path(output_dir, paste0(comparison_name, "_volcano.png")),
           volcano_plot, width = 10, height = 8, dpi = 300)

    return(res)
}

# Run all comparisons
cat("Running pairwise comparisons...\n")
cat("============================================\n")
all_results <- lapply(comparisons, function(contrast) {
    extract_results(dds, contrast, OUTPUT_DIR)
})
names(all_results) <- sapply(comparisons, function(x) paste0(x[1], "_vs_", x[2]))

# Heatmap of top DE genes across all comparisons
cat("Generating top DE genes heatmap...\n")

# Get union of top 50 DE genes from each comparison
top_genes <- unique(unlist(lapply(all_results, function(res) {
    res_sorted <- res[order(res$padj), ]
    head(rownames(res_sorted)[!is.na(res_sorted$padj)], 50)
})))

# Limit to top 100 for visualization
if (length(top_genes) > 100) {
    # Prioritize by minimum padj across comparisons
    min_padj <- sapply(top_genes, function(g) {
        min(sapply(all_results, function(res) {
            if (g %in% rownames(res)) res[g, "padj"] else 1
        }), na.rm = TRUE)
    })
    top_genes <- top_genes[order(min_padj)][1:100]
}

# Extract normalized counts for heatmap
heatmap_data <- assay(vst_data)[top_genes, ]
heatmap_data <- t(scale(t(heatmap_data)))  # Z-score normalization

# Heatmap annotation
annotation_col <- data.frame(Group = sample_info$group)
rownames(annotation_col) <- sample_info$sample_name
annotation_colors <- list(Group = c(Parental = "#E41A1C", Neg = "#377EB8",
                                     Pos = "#4DAF4A", KO = "#984EA3"))

pdf(file.path(OUTPUT_DIR, "top_DE_genes_heatmap.pdf"), width = 12, height = 14)
pheatmap(heatmap_data,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         fontsize_row = 6,
         clustering_method = "ward.D2",
         main = "Top Differentially Expressed Genes (Z-score)")
dev.off()

# Save DESeq2 object for downstream analysis
saveRDS(dds, file.path(OUTPUT_DIR, "dds_object.rds"))
saveRDS(vst_data, file.path(OUTPUT_DIR, "vst_object.rds"))

# Summary table
summary_df <- data.frame(
    Comparison = names(all_results),
    Total_Sig = sapply(all_results, function(x) sum(x$padj < 0.05, na.rm = TRUE)),
    Up = sapply(all_results, function(x) sum(x$padj < 0.05 & x$log2FoldChange > 0, na.rm = TRUE)),
    Down = sapply(all_results, function(x) sum(x$padj < 0.05 & x$log2FoldChange < 0, na.rm = TRUE))
)
write.csv(summary_df, file.path(OUTPUT_DIR, "DEG_summary.csv"), row.names = FALSE)

cat("\n============================================\n")
cat("DESeq2 Analysis Summary\n")
cat("============================================\n")
print(summary_df)

cat("\n============================================\n")
cat("DESeq2 analysis complete!\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
