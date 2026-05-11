#!/usr/bin/env Rscript
# =============================================================================
# 14e_srrm3_dosage_response.R
# Correlate per-sample SRRM3 expression with per-sample microexon PSI values
# across all 12 samples to test the dosage-response hypothesis.
#
# Rationale: Pos condition does not significantly upregulate SRRM3 mRNA
# (log2FC = -0.11, padj = 0.069), so condition labels alone are misleading.
# By treating each sample as a point on a continuous SRRM3 dosage axis, we
# directly test whether SRRM3 levels drive microexon inclusion.
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
})

# =============================================================================
# Configuration
# =============================================================================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
DESEQ_DIR   <- file.path(BASE_DIR, "results/04_deseq2")
OUTPUT_DIR   <- file.path(BASE_DIR, "results/14_todo_analysis/task4_dosage_response")

SAMPLE_NAMES <- c("Parental_1", "Parental_2", "Parental_3",
                   "Neg_1", "Neg_2", "Neg_3",
                   "Pos_1", "Pos_2", "Pos_3",
                   "KO_1", "KO_2", "KO_3")

GROUP_COLORS <- c(Parental = "#E41A1C", Neg = "#377EB8",
                  Pos = "#4DAF4A", KO = "#984EA3")

MIN_AVG_READS <- 10
MICROEXON_MAX <- 30
SMALL_EXON_MAX <- 50

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("14e: SRRM3 Dosage-Microexon PSI Correlation\n")
cat("========================================\n\n")

# =============================================================================
# Helper: vectorized parsing of comma-separated value columns
# Returns a numeric matrix (nrow = length(x), ncol = expected number of values)
# =============================================================================
parse_csv_column <- function(x, expected_n = 3) {
    splits <- strsplit(as.character(x), ",")
    lens <- vapply(splits, length, integer(1))
    mat <- matrix(NA_real_, nrow = length(x), ncol = expected_n)
    valid <- lens == expected_n
    if (any(valid)) {
        mat[valid, ] <- do.call(rbind, lapply(splits[valid], as.numeric))
    }
    mat
}

# =============================================================================
# Step 1: Build 12-sample PSI matrix for all SE events (vectorized)
# Approach: Combine IncLevel values from 3 vs-Parental comparisons.
# Same logic as extract_psi_matrix() in 6b_splicing_downstream.R
# =============================================================================
cat("Step 1: Building 12-sample PSI matrix...\n")

neg_file <- file.path(SPLICING_DIR, "Neg_vs_Parental/SE.MATS.JC.txt")
pos_file <- file.path(SPLICING_DIR, "Pos_vs_Parental/SE.MATS.JC.txt")
ko_file  <- file.path(SPLICING_DIR, "KO_vs_Parental/SE.MATS.JC.txt")

for (f in c(neg_file, pos_file, ko_file)) {
    if (!file.exists(f)) stop("Missing: ", f)
}

neg_df <- read.delim(neg_file, stringsAsFactors = FALSE)
pos_df <- read.delim(pos_file, stringsAsFactors = FALSE)
ko_df  <- read.delim(ko_file,  stringsAsFactors = FALSE)

cat(sprintf("  Neg_vs_Parental: %d SE events\n", nrow(neg_df)))
cat(sprintf("  Pos_vs_Parental: %d SE events\n", nrow(pos_df)))
cat(sprintf("  KO_vs_Parental:  %d SE events\n", nrow(ko_df)))

# --- Build coordinate-based event keys for matching across comparisons ---
# rMATS IDs are NOT consistent across separate runs, so we match by exon coordinates.
make_event_key <- function(df) {
    paste(df$chr, df$strand, df$exonStart_0base, df$exonEnd,
          df$upstreamEE, df$downstreamES, sep = ":")
}

neg_keys <- make_event_key(neg_df)
pos_keys <- make_event_key(pos_df)
ko_keys  <- make_event_key(ko_df)

cat(sprintf("  Unique event keys — Neg: %d, Pos: %d, KO: %d\n",
            length(unique(neg_keys)), length(unique(pos_keys)), length(unique(ko_keys))))

# --- Vectorized parse of Neg_vs_Parental ---
# IncLevel1 = Neg (3 samples), IncLevel2 = Parental (3 samples)
neg_inc1 <- parse_csv_column(neg_df$IncLevel1, 3)   # Neg PSI
neg_inc2 <- parse_csv_column(neg_df$IncLevel2, 3)   # Parental PSI

# Junction counts for coverage filter
neg_ijc1 <- parse_csv_column(neg_df$IJC_SAMPLE_1, 3)
neg_sjc1 <- parse_csv_column(neg_df$SJC_SAMPLE_1, 3)
neg_ijc2 <- parse_csv_column(neg_df$IJC_SAMPLE_2, 3)
neg_sjc2 <- parse_csv_column(neg_df$SJC_SAMPLE_2, 3)

# Valid rows: all 6 PSI values present
neg_valid <- complete.cases(neg_inc1) & complete.cases(neg_inc2)

# Build metadata for valid Neg_vs_Parental events
avg_reads_vec <- rowMeans(neg_ijc1 + neg_sjc1 + neg_ijc2 + neg_sjc2, na.rm = TRUE) / 2

meta_df <- data.frame(
    event_key   = neg_keys[neg_valid],
    gene_id     = neg_df$GeneID[neg_valid],
    gene_symbol = neg_df$geneSymbol[neg_valid],
    chr         = neg_df$chr[neg_valid],
    strand      = neg_df$strand[neg_valid],
    exon_start  = neg_df$exonStart_0base[neg_valid],
    exon_end    = neg_df$exonEnd[neg_valid],
    exon_length = neg_df$exonEnd[neg_valid] - neg_df$exonStart_0base[neg_valid],
    avg_reads   = avg_reads_vec[neg_valid],
    stringsAsFactors = FALSE
)

cat(sprintf("  Valid Neg_vs_Parental events: %d\n", sum(neg_valid)))

# --- Vectorized parse of Pos_vs_Parental (IncLevel1 = Pos) ---
pos_inc1  <- parse_csv_column(pos_df$IncLevel1, 3)
pos_valid <- complete.cases(pos_inc1)

# --- Vectorized parse of KO_vs_Parental (IncLevel1 = KO) ---
ko_inc1  <- parse_csv_column(ko_df$IncLevel1, 3)
ko_valid <- complete.cases(ko_inc1)

# --- Intersect by coordinate key to events present in all 3 comparisons ---
neg_valid_keys <- neg_keys[neg_valid]
pos_valid_keys <- pos_keys[pos_valid]
ko_valid_keys  <- ko_keys[ko_valid]
common_keys <- Reduce(intersect, list(neg_valid_keys, pos_valid_keys, ko_valid_keys))

cat(sprintf("  Events with complete 12-sample PSI: %d\n", length(common_keys)))

if (length(common_keys) < 50) {
    stop("Too few common events (", length(common_keys), "). Check rMATS output.")
}

# --- Assemble 12-sample PSI matrix using coordinate-based index lookup ---
neg_idx <- match(common_keys, neg_keys)
pos_idx <- match(common_keys, pos_keys)
ko_idx  <- match(common_keys, ko_keys)

psi_mat <- cbind(
    neg_inc2[neg_idx, ],   # Parental (cols 1-3)
    neg_inc1[neg_idx, ],   # Neg      (cols 4-6)
    pos_inc1[pos_idx, ],   # Pos      (cols 7-9)
    ko_inc1[ko_idx, ]      # KO       (cols 10-12)
)
rownames(psi_mat) <- common_keys
colnames(psi_mat) <- SAMPLE_NAMES

psi_mat <- psi_mat[complete.cases(psi_mat), ]

# --- Subset metadata to matched events and add size category ---
meta_df <- meta_df[match(rownames(psi_mat), meta_df$event_key), ]
meta_df$size_category <- ifelse(
    meta_df$exon_length <= MICROEXON_MAX, "microexon (<=30bp)",
    ifelse(meta_df$exon_length <= SMALL_EXON_MAX, "small (30-50bp)", "regular (>50bp)")
)
meta_df$size_category <- factor(meta_df$size_category,
    levels = c("microexon (<=30bp)", "small (30-50bp)", "regular (>50bp)"))

cat(sprintf("  After NA removal: %d events\n", nrow(psi_mat)))
cat("  Size distribution:\n")
print(table(meta_df$size_category))

# --- Apply junction read filter ---
kept <- meta_df$event_key[meta_df$avg_reads >= MIN_AVG_READS]
psi_mat <- psi_mat[kept, ]
meta_df <- meta_df %>% filter(event_key %in% kept)

cat(sprintf("\n  After avg_reads >= %d filter: %d events\n", MIN_AVG_READS, nrow(psi_mat)))
cat("  Size distribution after filter:\n")
print(table(meta_df$size_category))

micro_ids <- meta_df$event_key[meta_df$size_category == "microexon (<=30bp)"]
cat(sprintf("\n  Microexon events (<=30bp): %d\n", length(micro_ids)))

# --- Verification 1: expect >= 100 microexon events ---
if (length(micro_ids) < 100) {
    warning(sprintf("Only %d microexon events found (expected >= 100). Proceeding anyway.",
                    length(micro_ids)))
}

# =============================================================================
# Step 2: Extract per-sample SRRM3 / SRRM4 expression
# =============================================================================
cat("\nStep 2: Loading SRRM3/SRRM4 expression...\n")

counts <- read.csv(file.path(DESEQ_DIR, "normalized_counts.csv"),
                   stringsAsFactors = FALSE)

sample_cols <- colnames(counts)[3:ncol(counts)]
srrm3_row <- counts %>% filter(gene_symbol == "Srrm3")
srrm4_row <- counts %>% filter(gene_symbol == "Srrm4")

if (nrow(srrm3_row) == 0) stop("Srrm3 not found in normalized_counts.csv")
if (nrow(srrm4_row) == 0) stop("Srrm4 not found in normalized_counts.csv")

srrm3_expr <- as.numeric(srrm3_row[1, sample_cols])
srrm4_expr <- as.numeric(srrm4_row[1, sample_cols])
names(srrm3_expr) <- SAMPLE_NAMES
names(srrm4_expr) <- SAMPLE_NAMES

cat(sprintf("  SRRM3 range: %.1f - %.1f\n", min(srrm3_expr), max(srrm3_expr)))
cat(sprintf("  SRRM4 range: %.1f - %.1f\n", min(srrm4_expr), max(srrm4_expr)))

# =============================================================================
# Step 3: Compute per-event correlations (PSI ~ SRRM3 and PSI ~ SRRM4)
# =============================================================================
cat("\nStep 3: Computing per-event correlations...\n")

# Vectorized correlation: cor() on transposed PSI matrix vs expression vectors
# Then compute p-values from r and n using the t-distribution
n_samples <- ncol(psi_mat)
row_sds <- apply(psi_mat, 1, sd)
nonconst <- row_sds > 0  # skip zero-variance events

r_srrm3 <- rep(NA_real_, nrow(psi_mat))
r_srrm4 <- rep(NA_real_, nrow(psi_mat))
r_srrm3[nonconst] <- cor(t(psi_mat[nonconst, ]), srrm3_expr)
r_srrm4[nonconst] <- cor(t(psi_mat[nonconst, ]), srrm4_expr)

# p-value from Pearson r: t = r * sqrt((n-2)/(1-r^2)), df = n-2
pearson_p <- function(r, n) {
    t_stat <- r * sqrt((n - 2) / (1 - r^2))
    2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
}

cor_results <- data.frame(
    event_key = rownames(psi_mat),
    r_srrm3  = r_srrm3,
    p_srrm3  = pearson_p(r_srrm3, n_samples),
    r_srrm4  = r_srrm4,
    p_srrm4  = pearson_p(r_srrm4, n_samples),
    stringsAsFactors = FALSE
)

# Merge with metadata
all_cor <- cor_results %>%
    filter(!is.na(r_srrm3)) %>%
    left_join(meta_df, by = "event_key")

# Classify correlation strength
all_cor$srrm3_class <- ifelse(abs(all_cor$r_srrm3) > 0.7, "strong",
                        ifelse(abs(all_cor$r_srrm3) > 0.4, "moderate", "weak"))

cat(sprintf("  Events with computable correlations: %d\n", nrow(all_cor)))
cat("  SRRM3 correlation strength by size:\n")
print(table(all_cor$srrm3_class, all_cor$size_category))

# Microexon subset
micro_cor <- all_cor %>% filter(size_category == "microexon (<=30bp)")
cat(sprintf("\n  Microexon SRRM3 correlation summary:\n"))
cat(sprintf("    Mean r  = %.3f\n", mean(micro_cor$r_srrm3)))
cat(sprintf("    Median r = %.3f\n", median(micro_cor$r_srrm3)))
cat(sprintf("    Strong (|r|>0.7):   %d (%.1f%%)\n",
            sum(micro_cor$srrm3_class == "strong"),
            100 * mean(micro_cor$srrm3_class == "strong")))
cat(sprintf("    Moderate (0.4-0.7): %d (%.1f%%)\n",
            sum(micro_cor$srrm3_class == "moderate"),
            100 * mean(micro_cor$srrm3_class == "moderate")))
cat(sprintf("    Positive r: %d / %d (%.1f%%)\n",
            sum(micro_cor$r_srrm3 > 0), nrow(micro_cor),
            100 * mean(micro_cor$r_srrm3 > 0)))

# =============================================================================
# Step 4: Save output CSVs
# =============================================================================
cat("\nStep 4: Saving CSVs...\n")

# 4a. Microexon PSI matrix
micro_psi_out <- as.data.frame(psi_mat[micro_ids[micro_ids %in% rownames(psi_mat)], ])
micro_psi_out$event_key <- rownames(micro_psi_out)
micro_psi_out <- left_join(
    micro_psi_out %>% select(event_key, everything()),
    meta_df %>% select(event_key, gene_symbol, chr, strand,
                       exon_start, exon_end, exon_length),
    by = "event_key"
)
write.csv(micro_psi_out, file.path(OUTPUT_DIR, "microexon_psi_matrix.csv"),
          row.names = FALSE)
cat(sprintf("  microexon_psi_matrix.csv: %d events x 12 samples\n",
            nrow(micro_psi_out)))

# 4b. Dosage correlations (all events)
write.csv(
    all_cor %>% select(event_key, gene_symbol, chr, exon_start, exon_end,
                       exon_length, size_category, r_srrm3, p_srrm3,
                       r_srrm4, p_srrm4, srrm3_class),
    file.path(OUTPUT_DIR, "dosage_correlations.csv"),
    row.names = FALSE
)
cat(sprintf("  dosage_correlations.csv: %d events\n", nrow(all_cor)))

# 4c. Dosage summary by size (built alongside Plot 5 below)

# =============================================================================
# Step 5: Plots
# =============================================================================
cat("\nStep 5: Generating plots...\n")

# ---- Plot 1: Global scatter — mean microexon PSI vs SRRM3 expression ----
cat("  Plot 1: Global scatter...\n")

mean_micro_psi <- colMeans(psi_mat[micro_ids[micro_ids %in% rownames(psi_mat)], ,
                                    drop = FALSE], na.rm = TRUE)

scatter_df <- data.frame(
    sample   = SAMPLE_NAMES,
    group    = factor(gsub("_[0-9]+$", "", SAMPLE_NAMES),
                      levels = c("Parental", "Neg", "Pos", "KO")),
    srrm3    = srrm3_expr,
    mean_psi = mean_micro_psi
)

global_cor <- cor.test(scatter_df$srrm3, scatter_df$mean_psi)

p1 <- ggplot(scatter_df, aes(x = srrm3, y = mean_psi, color = group)) +
    geom_point(size = 4) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "grey40", linetype = "dashed") +
    scale_color_manual(values = GROUP_COLORS) +
    labs(
        title = "Mean Microexon PSI vs SRRM3 Expression",
        subtitle = sprintf(
            "Pearson r = %.3f, p = %.2e  (n = 12 samples, %d microexons)",
            global_cor$estimate, global_cor$p.value, length(micro_ids)),
        x = "SRRM3 Normalized Expression",
        y = "Mean Microexon PSI",
        color = "Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "plot1_global_scatter.pdf"), p1,
       width = 8, height = 6)

# ---- Plot 1b: Global scatter — mean microexon PSI vs SRRM4 expression ----
cat("  Plot 1b: Global scatter (SRRM4)...\n")

scatter_df$srrm4 <- srrm4_expr
global_cor_srrm4 <- cor.test(scatter_df$srrm4, scatter_df$mean_psi)

p1b <- ggplot(scatter_df, aes(x = srrm4, y = mean_psi, color = group)) +
    geom_point(size = 4) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "grey40", linetype = "dashed") +
    scale_color_manual(values = GROUP_COLORS) +
    labs(
        title = "Mean Microexon PSI vs SRRM4 Expression",
        subtitle = sprintf(
            "Pearson r = %.3f, p = %.2e  (n = 12 samples, %d microexons)",
            global_cor_srrm4$estimate, global_cor_srrm4$p.value, length(micro_ids)),
        x = "SRRM4 Normalized Expression",
        y = "Mean Microexon PSI",
        color = "Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "plot1b_global_scatter_srrm4.pdf"), p1b,
       width = 8, height = 6)

# ---- Plot 2: Per-event correlation histogram ----
cat("  Plot 2: Correlation histogram...\n")

p2 <- ggplot(micro_cor, aes(x = r_srrm3)) +
    geom_histogram(bins = 30, fill = "#E41A1C", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = median(micro_cor$r_srrm3), linetype = "solid",
               color = "darkblue", linewidth = 0.8) +
    annotate("text",
             x = median(micro_cor$r_srrm3), y = Inf, vjust = 2,
             label = sprintf("median r = %.2f", median(micro_cor$r_srrm3)),
             color = "darkblue", size = 4, hjust = -0.1) +
    labs(
        title = "Distribution of SRRM3-PSI Correlations Across Microexon Events",
        subtitle = sprintf("n = %d microexons (<=30bp), 12 samples per event",
                           nrow(micro_cor)),
        x = "Pearson r (PSI vs SRRM3 expression)",
        y = "Number of Events"
    ) +
    theme_minimal(base_size = 13)

ggsave(file.path(OUTPUT_DIR, "plot2_correlation_histogram.pdf"), p2,
       width = 8, height = 5)

# ---- Plot 2b: Per-event SRRM4 correlation histogram ----
cat("  Plot 2b: Correlation histogram (SRRM4)...\n")

p2b <- ggplot(micro_cor, aes(x = r_srrm4)) +
    geom_histogram(bins = 30, fill = "#377EB8", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = median(micro_cor$r_srrm4, na.rm = TRUE),
               linetype = "solid", color = "darkred", linewidth = 0.8) +
    annotate("text",
             x = median(micro_cor$r_srrm4, na.rm = TRUE), y = Inf, vjust = 2,
             label = sprintf("median r = %.2f",
                             median(micro_cor$r_srrm4, na.rm = TRUE)),
             color = "darkred", size = 4, hjust = -0.1) +
    labs(
        title = "Distribution of SRRM4-PSI Correlations Across Microexon Events",
        subtitle = sprintf("n = %d microexons (<=30bp), 12 samples per event",
                           nrow(micro_cor)),
        x = "Pearson r (PSI vs SRRM4 expression)",
        y = "Number of Events"
    ) +
    theme_minimal(base_size = 13)

ggsave(file.path(OUTPUT_DIR, "plot2b_correlation_histogram_srrm4.pdf"), p2b,
       width = 8, height = 5)

# ---- Plot 3: Top correlated events heatmap ----
cat("  Plot 3: Top correlated heatmap...\n")

top_n <- min(30, nrow(micro_cor))
top_events <- micro_cor %>% arrange(desc(r_srrm3)) %>% head(top_n)
top_psi <- psi_mat[top_events$event_key, , drop = FALSE]

# Order samples by SRRM3 expression (low → high)
sample_order <- order(srrm3_expr)
top_psi_ordered <- top_psi[, sample_order]

# Column annotation
col_annot <- data.frame(
    Group = factor(gsub("_[0-9]+$", "", SAMPLE_NAMES[sample_order]),
                   levels = c("Parental", "Neg", "Pos", "KO")),
    SRRM3 = srrm3_expr[sample_order]
)
rownames(col_annot) <- SAMPLE_NAMES[sample_order]

# Row labels: gene_symbol (r = x.xx)
row_labels <- sprintf("%s (r=%.2f)",
                      gsub('"', '', top_events$gene_symbol),
                      top_events$r_srrm3)
rownames(top_psi_ordered) <- row_labels

annot_colors <- list(Group = GROUP_COLORS)

pdf(file.path(OUTPUT_DIR, "plot3_top_correlated_heatmap.pdf"),
    width = 10, height = 10)
pheatmap(top_psi_ordered,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_col = col_annot,
         annotation_colors = annot_colors,
         color = colorRampPalette(c("#3C5488", "white", "#DC0000"))(100),
         main = sprintf("Top %d SRRM3-Correlated Microexons\n(samples ordered by SRRM3 expression)",
                        top_n),
         fontsize_row = 8,
         fontsize_col = 9)
dev.off()

# ---- Plot 4: SRRM3 vs SRRM4 correlation comparison ----
cat("  Plot 4: SRRM3 vs SRRM4 comparison...\n")

srrm34_cor <- cor.test(micro_cor$r_srrm3, micro_cor$r_srrm4)

p4 <- ggplot(micro_cor, aes(x = r_srrm3, y = r_srrm4)) +
    geom_point(alpha = 0.5, size = 1.5, color = "#3C5488") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", linewidth = 0.8) +
    labs(
        title = "SRRM3 vs SRRM4 Correlation with Microexon PSI",
        subtitle = sprintf("r(SRRM3) vs r(SRRM4): Pearson = %.3f, p = %.2e",
                           srrm34_cor$estimate, srrm34_cor$p.value),
        x = "Per-event Pearson r (PSI ~ SRRM3)",
        y = "Per-event Pearson r (PSI ~ SRRM4)"
    ) +
    theme_minimal(base_size = 13) +
    coord_equal()

ggsave(file.path(OUTPUT_DIR, "plot4_srrm3_vs_srrm4.pdf"), p4,
       width = 7, height = 7)

# ---- Plot 5: Dosage-response by exon size category ----
cat("  Plot 5: Dosage-response by size category...\n")

size_cats <- levels(meta_df$size_category)
size_summary_list <- list()

for (sc in size_cats) {
    sc_ids <- meta_df$event_key[meta_df$size_category == sc]
    sc_ids <- sc_ids[sc_ids %in% rownames(psi_mat)]
    if (length(sc_ids) < 2) next

    mean_psi_sc <- colMeans(psi_mat[sc_ids, , drop = FALSE], na.rm = TRUE)
    size_summary_list[[sc]] <- data.frame(
        sample        = SAMPLE_NAMES,
        group         = factor(gsub("_[0-9]+$", "", SAMPLE_NAMES),
                               levels = c("Parental", "Neg", "Pos", "KO")),
        srrm3         = srrm3_expr,
        mean_psi      = mean_psi_sc,
        size_category = sc,
        n_events      = length(sc_ids),
        stringsAsFactors = FALSE
    )
}
size_summary <- bind_rows(size_summary_list)
size_summary$size_category <- factor(size_summary$size_category, levels = size_cats)

# Save CSV (Step 4c)
write.csv(size_summary, file.path(OUTPUT_DIR, "dosage_summary_by_size.csv"),
          row.names = FALSE)
cat(sprintf("  dosage_summary_by_size.csv: %d rows\n", nrow(size_summary)))

# Compute per-category correlations for facet labels
size_cors <- size_summary %>%
    group_by(size_category) %>%
    summarize(r = cor(srrm3, mean_psi), n = first(n_events), .groups = "drop")

size_summary <- size_summary %>%
    left_join(size_cors %>% select(size_category, r, n), by = "size_category") %>%
    mutate(facet_label = sprintf("%s\n(r=%.2f, n=%d)", size_category, r, n))

# Preserve factor ordering for facets
facet_levels <- size_summary %>%
    distinct(size_category, facet_label) %>%
    arrange(size_category) %>%
    pull(facet_label)
size_summary$facet_label <- factor(size_summary$facet_label, levels = facet_levels)

p5 <- ggplot(size_summary, aes(x = srrm3, y = mean_psi, color = group)) +
    geom_point(size = 3) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "grey40", linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_color_manual(values = GROUP_COLORS) +
    labs(
        title = "SRRM3 Dosage-Response by Exon Size Category",
        subtitle = "Mean PSI vs SRRM3 normalized expression (12 samples)",
        x = "SRRM3 Normalized Expression",
        y = "Mean PSI",
        color = "Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(size = 10))

ggsave(file.path(OUTPUT_DIR, "plot5_dosage_by_size.pdf"), p5,
       width = 12, height = 5)

# ---- Plot 5b: SRRM4 Dosage-response by exon size category ----
cat("  Plot 5b: Dosage-response by size category (SRRM4)...\n")

size_summary$srrm4 <- srrm4_expr[match(size_summary$sample, SAMPLE_NAMES)]

size_cors_srrm4 <- size_summary %>%
    group_by(size_category) %>%
    summarize(r = cor(srrm4, mean_psi), n = first(n_events), .groups = "drop")

size_summary <- size_summary %>%
    left_join(size_cors_srrm4 %>%
                  select(size_category, r_s4 = r, n_s4 = n),
              by = "size_category") %>%
    mutate(facet_label_s4 = sprintf("%s\n(r=%.2f, n=%d)", size_category, r_s4, n_s4))

facet_levels_s4 <- size_summary %>%
    distinct(size_category, facet_label_s4) %>%
    arrange(size_category) %>%
    pull(facet_label_s4)
size_summary$facet_label_s4 <- factor(size_summary$facet_label_s4,
                                       levels = facet_levels_s4)

p5b <- ggplot(size_summary, aes(x = srrm4, y = mean_psi, color = group)) +
    geom_point(size = 3) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "grey40", linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~ facet_label_s4, scales = "free_y") +
    scale_color_manual(values = GROUP_COLORS) +
    labs(
        title = "SRRM4 Dosage-Response by Exon Size Category",
        subtitle = "Mean PSI vs SRRM4 normalized expression (12 samples)",
        x = "SRRM4 Normalized Expression",
        y = "Mean PSI",
        color = "Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(size = 10))

ggsave(file.path(OUTPUT_DIR, "plot5b_dosage_by_size_srrm4.pdf"), p5b,
       width = 12, height = 5)

# =============================================================================
# Verification & Summary
# =============================================================================
cat("\n========================================\n")
cat("Verification & Summary\n")
cat("========================================\n")

cat(sprintf("Total SE events (12-sample, avg_reads>=%d): %d\n",
            MIN_AVG_READS, nrow(psi_mat)))
cat(sprintf("  Microexons (<=30bp):  %d\n",
            sum(meta_df$size_category == "microexon (<=30bp)")))
cat(sprintf("  Small (30-50bp):      %d\n",
            sum(meta_df$size_category == "small (30-50bp)")))
cat(sprintf("  Regular (>50bp):      %d\n",
            sum(meta_df$size_category == "regular (>50bp)")))

cat(sprintf("\nGlobal SRRM3 ~ mean microexon PSI: r = %.3f, p = %.2e\n",
            global_cor$estimate, global_cor$p.value))
cat(sprintf("Microexon per-event SRRM3 correlation: median r = %.3f\n",
            median(micro_cor$r_srrm3)))
cat(sprintf("  Strong  (|r|>0.7):   %d events\n",
            sum(abs(micro_cor$r_srrm3) > 0.7)))
cat(sprintf("  Moderate (0.4-0.7):  %d events\n",
            sum(abs(micro_cor$r_srrm3) > 0.4 & abs(micro_cor$r_srrm3) <= 0.7)))
cat(sprintf("  Weak     (<0.4):     %d events\n",
            sum(abs(micro_cor$r_srrm3) <= 0.4)))

# Verification 3: majority of microexon events should have r > 0
pct_positive <- 100 * mean(micro_cor$r_srrm3 > 0)
cat(sprintf("\n  Positive correlations: %.1f%% of microexons", pct_positive))
if (pct_positive > 50) {
    cat(" [PASS]\n")
} else {
    cat(" [UNEXPECTED — expected majority positive]\n")
}

# Verification 4: microexons should have stronger correlations than regular exons
reg_cor <- all_cor %>% filter(size_category == "regular (>50bp)")
if (nrow(reg_cor) > 0) {
    cat(sprintf("  Median |r| — microexons: %.3f, regular: %.3f",
                median(abs(micro_cor$r_srrm3)),
                median(abs(reg_cor$r_srrm3))))
    if (median(abs(micro_cor$r_srrm3)) > median(abs(reg_cor$r_srrm3))) {
        cat(" [PASS — microexons more SRRM3-sensitive]\n")
    } else {
        cat(" [UNEXPECTED — regular exons equally/more sensitive]\n")
    }
}

# Size-dependent correlation summary
cat("\n  Per-size-category SRRM3 correlation (mean PSI ~ SRRM3):\n")
for (i in seq_len(nrow(size_cors))) {
    cat(sprintf("    %s: r = %.3f (n = %d events)\n",
                size_cors$size_category[i], size_cors$r[i], size_cors$n[i]))
}

# ---- SRRM4 summary statistics ----
cat(sprintf("\nGlobal SRRM4 ~ mean microexon PSI: r = %.3f, p = %.2e\n",
            global_cor_srrm4$estimate, global_cor_srrm4$p.value))
cat(sprintf("Microexon per-event SRRM4 correlation: median r = %.3f\n",
            median(micro_cor$r_srrm4, na.rm = TRUE)))
cat(sprintf("  Strong  (|r|>0.7):   %d events\n",
            sum(abs(micro_cor$r_srrm4) > 0.7, na.rm = TRUE)))
cat(sprintf("  Moderate (0.4-0.7):  %d events\n",
            sum(abs(micro_cor$r_srrm4) > 0.4 & abs(micro_cor$r_srrm4) <= 0.7,
                na.rm = TRUE)))
cat(sprintf("  Weak     (<0.4):     %d events\n",
            sum(abs(micro_cor$r_srrm4) <= 0.4, na.rm = TRUE)))
pct_positive_s4 <- 100 * mean(micro_cor$r_srrm4 > 0, na.rm = TRUE)
cat(sprintf("  Positive correlations: %.1f%% of microexons\n", pct_positive_s4))

cat("\n  Per-size-category SRRM4 correlation (mean PSI ~ SRRM4):\n")
for (i in seq_len(nrow(size_cors_srrm4))) {
    cat(sprintf("    %s: r = %.3f (n = %d events)\n",
                size_cors_srrm4$size_category[i],
                size_cors_srrm4$r[i], size_cors_srrm4$n[i]))
}

# ---- SRRM3 vs SRRM4 comparison summary ----
cat("\n  SRRM3 vs SRRM4 comparison:\n")
cat(sprintf("    Global r (microexon PSI): SRRM3 = %.3f, SRRM4 = %.3f\n",
            global_cor$estimate, global_cor_srrm4$estimate))
cat(sprintf("    Per-event median |r|:     SRRM3 = %.3f, SRRM4 = %.3f\n",
            median(abs(micro_cor$r_srrm3)),
            median(abs(micro_cor$r_srrm4), na.rm = TRUE)))

cat(sprintf("\nOutputs saved to: %s\n", OUTPUT_DIR))
cat("  - microexon_psi_matrix.csv\n")
cat("  - dosage_correlations.csv\n")
cat("  - dosage_summary_by_size.csv\n")
cat("  - plot1_global_scatter.pdf\n")
cat("  - plot1b_global_scatter_srrm4.pdf\n")
cat("  - plot2_correlation_histogram.pdf\n")
cat("  - plot2b_correlation_histogram_srrm4.pdf\n")
cat("  - plot3_top_correlated_heatmap.pdf\n")
cat("  - plot4_srrm3_vs_srrm4.pdf\n")
cat("  - plot5_dosage_by_size.pdf\n")
cat("  - plot5b_dosage_by_size_srrm4.pdf\n")
cat("\n14e complete.\n")
