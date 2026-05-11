#!/usr/bin/env Rscript
# Novel SRRM3 exon PSI analysis
# Location: SRRM3_novel_exon/scripts/02_analyze_PSI.R

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

# Separate analysis directory
analysis_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069/SRRM3_novel_exon"
results_dir <- file.path(analysis_dir, "results")

cat("=== SRRM3 Novel Exon PSI Analysis ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Read junction counts
junction_file <- file.path(results_dir, "junction_counts_summary.csv")
if (!file.exists(junction_file)) {
    stop("Junction counts file not found. Run 01_extract_junctions.sh first.")
}

junctions <- fread(junction_file)
cat("Loaded junction counts for", nrow(junctions), "samples\n\n")

# Display raw counts
cat("Raw junction counts:\n")
print(junctions)
cat("\n")

# Calculate PSI (Percent Spliced In)
# PSI = Inclusion_reads / (Inclusion_reads + Skipping_reads)
# For cassette exons with two inclusion junctions:
# PSI = (Inc1 + Inc2) / (Inc1 + Inc2 + 2*Skip)
# This weights inclusion and skipping equally

junctions[, `:=`(
    # Average inclusion reads (geometric mean of the two junctions, or arithmetic)
    inc_avg = (inc1_reads + inc2_reads) / 2,
    # Total inclusion evidence
    inc_total = inc1_reads + inc2_reads
)]

# Calculate PSI using the standard formula for cassette exons
# PSI = average(IR1/(IR1+ER), IR2/(IR2+ER)) where IR = inclusion reads, ER = exclusion reads
# Simplified: PSI = (Inc1 + Inc2) / (Inc1 + Inc2 + 2*Skip)
junctions[, PSI := ifelse(
    (inc_total + 2 * skip_reads) > 0,
    inc_total / (inc_total + 2 * skip_reads),
    NA_real_
)]

# Also calculate a simpler PSI using minimum of the two inclusion junctions
# This is more conservative
junctions[, PSI_conservative := ifelse(
    (pmin(inc1_reads, inc2_reads) + skip_reads) > 0,
    pmin(inc1_reads, inc2_reads) / (pmin(inc1_reads, inc2_reads) + skip_reads),
    NA_real_
)]

# Set group factor order
junctions[, group := factor(group, levels = c("Parental", "Neg", "Pos", "KO"))]

cat("PSI values by sample:\n")
print(junctions[, .(sample, group, inc1_reads, inc2_reads, skip_reads, novel_coverage, PSI, PSI_conservative)])
cat("\n")

# Summary statistics by group
cat("Summary statistics by group:\n")
summary_stats <- junctions[, .(
    n = .N,
    mean_PSI = mean(PSI, na.rm = TRUE),
    sd_PSI = sd(PSI, na.rm = TRUE),
    mean_inc1 = mean(inc1_reads),
    mean_inc2 = mean(inc2_reads),
    mean_skip = mean(skip_reads),
    mean_coverage = mean(novel_coverage)
), by = group]
print(summary_stats)
cat("\n")

# Save PSI results
psi_output <- junctions[, .(sample, group, inc1_reads, inc2_reads, skip_reads, novel_coverage, PSI, PSI_conservative)]
fwrite(psi_output, file.path(results_dir, "novel_exon_PSI.csv"))
cat("Saved PSI values to:", file.path(results_dir, "novel_exon_PSI.csv"), "\n\n")

# Statistical tests
cat("=== Statistical Tests ===\n\n")

# Store results
stat_results <- list()

# Test 1: Pos vs Parental
if (sum(junctions$group == "Pos") >= 2 && sum(junctions$group == "Parental") >= 2) {
    test_pos_parent <- wilcox.test(
        junctions[group == "Pos", PSI],
        junctions[group == "Parental", PSI],
        alternative = "greater"
    )
    cat("Pos vs Parental (one-sided, Pos > Parental):\n")
    cat("  Wilcoxon p-value:", format(test_pos_parent$p.value, digits = 4), "\n")
    cat("  Pos mean PSI:", round(mean(junctions[group == "Pos", PSI], na.rm = TRUE), 4), "\n")
    cat("  Parental mean PSI:", round(mean(junctions[group == "Parental", PSI], na.rm = TRUE), 4), "\n\n")
    stat_results[["Pos_vs_Parental"]] <- data.table(
        comparison = "Pos_vs_Parental",
        test = "Wilcoxon",
        alternative = "greater",
        p_value = test_pos_parent$p.value
    )
}

# Test 2: Pos vs Neg
if (sum(junctions$group == "Pos") >= 2 && sum(junctions$group == "Neg") >= 2) {
    test_pos_neg <- wilcox.test(
        junctions[group == "Pos", PSI],
        junctions[group == "Neg", PSI],
        alternative = "greater"
    )
    cat("Pos vs Neg (one-sided, Pos > Neg):\n")
    cat("  Wilcoxon p-value:", format(test_pos_neg$p.value, digits = 4), "\n")
    cat("  Pos mean PSI:", round(mean(junctions[group == "Pos", PSI], na.rm = TRUE), 4), "\n")
    cat("  Neg mean PSI:", round(mean(junctions[group == "Neg", PSI], na.rm = TRUE), 4), "\n\n")
    stat_results[["Pos_vs_Neg"]] <- data.table(
        comparison = "Pos_vs_Neg",
        test = "Wilcoxon",
        alternative = "greater",
        p_value = test_pos_neg$p.value
    )
}

# Test 3: Pos vs KO
if (sum(junctions$group == "Pos") >= 2 && sum(junctions$group == "KO") >= 2) {
    test_pos_ko <- wilcox.test(
        junctions[group == "Pos", PSI],
        junctions[group == "KO", PSI],
        alternative = "greater"
    )
    cat("Pos vs KO (one-sided, Pos > KO):\n")
    cat("  Wilcoxon p-value:", format(test_pos_ko$p.value, digits = 4), "\n")
    cat("  Pos mean PSI:", round(mean(junctions[group == "Pos", PSI], na.rm = TRUE), 4), "\n")
    cat("  KO mean PSI:", round(mean(junctions[group == "KO", PSI], na.rm = TRUE), 4), "\n\n")
    stat_results[["Pos_vs_KO"]] <- data.table(
        comparison = "Pos_vs_KO",
        test = "Wilcoxon",
        alternative = "greater",
        p_value = test_pos_ko$p.value
    )
}

# Test 4: Overall Kruskal-Wallis
kw_test <- kruskal.test(PSI ~ group, data = junctions)
cat("Kruskal-Wallis test (all groups):\n")
cat("  Chi-squared:", round(kw_test$statistic, 4), "\n")
cat("  df:", kw_test$parameter, "\n")
cat("  p-value:", format(kw_test$p.value, digits = 4), "\n\n")
stat_results[["Kruskal_Wallis"]] <- data.table(
    comparison = "All_groups",
    test = "Kruskal-Wallis",
    alternative = "two-sided",
    p_value = kw_test$p.value
)

# Save statistical results
if (length(stat_results) > 0) {
    stat_dt <- rbindlist(stat_results)
    fwrite(stat_dt, file.path(results_dir, "statistical_tests.csv"))
    cat("Saved statistical tests to:", file.path(results_dir, "statistical_tests.csv"), "\n\n")
}

# Create visualization
cat("=== Creating Visualizations ===\n\n")

# Color palette
group_colors <- c(
    "Parental" = "#4DAF4A",  # Green
    "Neg" = "#377EB8",       # Blue
    "Pos" = "#E41A1C",       # Red
    "KO" = "#984EA3"         # Purple
)

# Plot 1: PSI by group (boxplot with points)
p1 <- ggplot(junctions, aes(x = group, y = PSI, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
    scale_fill_manual(values = group_colors) +
    labs(
        title = "SRRM3 Novel Exon Inclusion (PSI)",
        subtitle = "chr5:135,898,574-135,898,652 (79 bp)",
        x = "Sample Group",
        y = "Percent Spliced In (PSI)",
        caption = "PSI = (Inc1 + Inc2) / (Inc1 + Inc2 + 2*Skip)"
    ) +
    theme_bw(base_size = 14) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
    ) +
    ylim(0, 1)

ggsave(file.path(results_dir, "PSI_by_group.pdf"), p1, width = 8, height = 6)
ggsave(file.path(results_dir, "PSI_by_group.png"), p1, width = 8, height = 6, dpi = 150)
cat("Saved PSI boxplot to:", file.path(results_dir, "PSI_by_group.pdf"), "\n")

# Plot 2: Junction read counts
junctions_long <- melt(
    junctions,
    id.vars = c("sample", "group"),
    measure.vars = c("inc1_reads", "inc2_reads", "skip_reads"),
    variable.name = "junction_type",
    value.name = "read_count"
)
junctions_long[, junction_type := factor(junction_type,
    levels = c("inc1_reads", "inc2_reads", "skip_reads"),
    labels = c("Inclusion 1\n(upstream->novel)", "Inclusion 2\n(novel->downstream)", "Skipping\n(upstream->downstream)")
)]

p2 <- ggplot(junctions_long, aes(x = group, y = read_count, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    facet_wrap(~junction_type, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(
        title = "SRRM3 Novel Exon Junction Read Counts",
        x = "Sample Group",
        y = "Read Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

ggsave(file.path(results_dir, "junction_counts_by_group.pdf"), p2, width = 12, height = 5)
ggsave(file.path(results_dir, "junction_counts_by_group.png"), p2, width = 12, height = 5, dpi = 150)
cat("Saved junction counts plot to:", file.path(results_dir, "junction_counts_by_group.pdf"), "\n")

# Plot 3: Novel exon coverage
p3 <- ggplot(junctions, aes(x = group, y = novel_coverage, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
    scale_fill_manual(values = group_colors) +
    labs(
        title = "Reads Overlapping Novel SRRM3 Exon",
        subtitle = "chr5:135,898,574-135,898,652",
        x = "Sample Group",
        y = "Read Count"
    ) +
    theme_bw(base_size = 14) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
    )

ggsave(file.path(results_dir, "novel_exon_coverage.pdf"), p3, width = 8, height = 6)
ggsave(file.path(results_dir, "novel_exon_coverage.png"), p3, width = 8, height = 6, dpi = 150)
cat("Saved coverage plot to:", file.path(results_dir, "novel_exon_coverage.pdf"), "\n")

cat("\n=== Analysis Complete ===\n")
cat("Output files in:", results_dir, "\n")
