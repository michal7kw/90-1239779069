#!/usr/bin/env Rscript
# =============================================================================
# Microexon Functional Annotation
# Project: 90-1239779069
# Analysis: Classify microexons by genomic position (5'UTR, CDS, 3'UTR)
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(RColorBrewer)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
OUTPUT_DIR <- file.path(BASE_DIR, "results/09_microexon_extended/functional_annotation")
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

# Thresholds
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1

# Exon size bins
MICROEXON_MAX <- 30      # 0-30bp: microexons
SMALL_EXON_MAX <- 50     # 30-50bp: small exons

cat("============================================\n")
cat("Microexon Functional Annotation\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("GTF file:", GTF_FILE, "\n")
cat("============================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define comparisons
parental_comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")

# Color palettes
position_colors <- c("5UTR" = "#e74c3c", "CDS" = "#3498db", "3UTR" = "#2ecc71",
                     "UTR" = "#f39c12", "Unknown" = "#95a5a6")
biotype_colors <- c("protein_coding" = "#3498db", "lncRNA" = "#e74c3c",
                    "processed_pseudogene" = "#9b59b6", "other" = "#95a5a6")

# =============================================================================
# Function: Load and parse GTF
# =============================================================================
load_gtf_annotations <- function(gtf_file) {
    cat("Loading GTF annotations...\n")

    if (!file.exists(gtf_file)) {
        stop(paste("GTF file not found:", gtf_file))
    }

    # Import GTF
    gtf <- import(gtf_file)

    # Extract different feature types
    annotations <- list(
        genes = gtf[gtf$type == "gene"],
        transcripts = gtf[gtf$type == "transcript"],
        exons = gtf[gtf$type == "exon"],
        CDS = gtf[gtf$type == "CDS"],
        UTR5 = gtf[gtf$type == "five_prime_UTR" | gtf$type == "5UTR"],
        UTR3 = gtf[gtf$type == "three_prime_UTR" | gtf$type == "3UTR"],
        UTR = gtf[gtf$type == "UTR"]  # Some GTFs use generic UTR
    )

    cat(paste0("  Loaded ", length(annotations$genes), " genes\n"))
    cat(paste0("  Loaded ", length(annotations$transcripts), " transcripts\n"))
    cat(paste0("  Loaded ", length(annotations$exons), " exons\n"))
    cat(paste0("  Loaded ", length(annotations$CDS), " CDS features\n"))
    cat(paste0("  Loaded ", length(annotations$UTR5), " 5'UTR features\n"))
    cat(paste0("  Loaded ", length(annotations$UTR3), " 3'UTR features\n"))

    return(annotations)
}

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

    # Calculate exon length
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
# Function: Annotate exon genomic position
# =============================================================================
annotate_exon_position <- function(se_df, gtf_annotations) {
    cat("Annotating exon genomic positions...\n")

    # Create GRanges from SE events
    # Handle chromosome naming (add 'chr' if needed)
    se_chr <- se_df$chr

    se_gr <- GRanges(
        seqnames = se_chr,
        ranges = IRanges(start = se_df$exonStart_0base + 1,  # Convert 0-based to 1-based
                        end = se_df$exonEnd),
        strand = se_df$strand
    )

    # Initialize position annotation
    se_df$genomic_position <- "Unknown"

    # Check overlap with CDS
    if (length(gtf_annotations$CDS) > 0) {
        cds_overlaps <- overlapsAny(se_gr, gtf_annotations$CDS)
        se_df$genomic_position[cds_overlaps] <- "CDS"
    }

    # Check overlap with 5'UTR
    if (length(gtf_annotations$UTR5) > 0) {
        utr5_overlaps <- overlapsAny(se_gr, gtf_annotations$UTR5)
        se_df$genomic_position[utr5_overlaps & se_df$genomic_position == "Unknown"] <- "5UTR"
    }

    # Check overlap with 3'UTR
    if (length(gtf_annotations$UTR3) > 0) {
        utr3_overlaps <- overlapsAny(se_gr, gtf_annotations$UTR3)
        se_df$genomic_position[utr3_overlaps & se_df$genomic_position == "Unknown"] <- "3UTR"
    }

    # Check generic UTR if specific 5'/3' not available
    if (length(gtf_annotations$UTR) > 0) {
        utr_overlaps <- overlapsAny(se_gr, gtf_annotations$UTR)
        se_df$genomic_position[utr_overlaps & se_df$genomic_position == "Unknown"] <- "UTR"
    }

    # Summary
    cat("  Position annotation summary:\n")
    print(table(se_df$genomic_position))

    return(se_df)
}

# =============================================================================
# Function: Get gene biotype from GTF
# =============================================================================
get_gene_biotype <- function(se_df, gtf_annotations) {
    cat("Annotating gene biotypes...\n")

    genes <- gtf_annotations$genes

    # Extract gene biotypes
    if ("gene_biotype" %in% names(mcols(genes))) {
        biotype_col <- "gene_biotype"
    } else if ("gene_type" %in% names(mcols(genes))) {
        biotype_col <- "gene_type"
    } else {
        warning("No gene biotype column found in GTF")
        se_df$gene_biotype <- "Unknown"
        return(se_df)
    }

    # Create biotype lookup
    gene_biotypes <- setNames(mcols(genes)[[biotype_col]],
                              gsub("\"", "", mcols(genes)$gene_id))

    # Clean gene IDs (remove quotes and version numbers)
    se_df$clean_gene_id <- gsub("\"", "", se_df$GeneID)
    se_df$clean_gene_id <- gsub("\\..*", "", se_df$clean_gene_id)

    # Match biotypes
    se_df$gene_biotype <- gene_biotypes[se_df$clean_gene_id]
    se_df$gene_biotype[is.na(se_df$gene_biotype)] <- "Unknown"

    # Simplify biotypes
    se_df$gene_biotype_simple <- case_when(
        se_df$gene_biotype == "protein_coding" ~ "protein_coding",
        grepl("lncRNA|lincRNA", se_df$gene_biotype) ~ "lncRNA",
        grepl("pseudogene", se_df$gene_biotype) ~ "pseudogene",
        TRUE ~ "other"
    )

    cat("  Gene biotype summary:\n")
    print(table(se_df$gene_biotype_simple))

    return(se_df)
}

# =============================================================================
# Function: Calculate relative position within transcript
# =============================================================================
calculate_relative_position <- function(se_df, gtf_annotations) {
    cat("Calculating relative exon positions...\n")

    transcripts <- gtf_annotations$transcripts

    if (length(transcripts) == 0) {
        se_df$relative_position <- NA
        return(se_df)
    }

    # For each exon, find its relative position in the transcript (0-1)
    # 0 = start of transcript, 1 = end of transcript

    se_df$relative_position <- NA

    for (i in seq_len(nrow(se_df))) {
        # Get gene ID
        gene_id <- gsub("\"", "", se_df$GeneID[i])

        # Find transcripts for this gene
        gene_tx <- transcripts[grepl(gene_id, mcols(transcripts)$gene_id)]

        if (length(gene_tx) > 0) {
            # Use the first/canonical transcript
            tx <- gene_tx[1]
            tx_start <- start(tx)
            tx_end <- end(tx)
            tx_length <- tx_end - tx_start

            exon_mid <- (se_df$exonStart_0base[i] + se_df$exonEnd[i]) / 2

            if (as.character(strand(tx)) == "-") {
                # Reverse for minus strand
                se_df$relative_position[i] <- (tx_end - exon_mid) / tx_length
            } else {
                se_df$relative_position[i] <- (exon_mid - tx_start) / tx_length
            }
        }
    }

    return(se_df)
}

# =============================================================================
# Function: Plot position distribution
# =============================================================================
plot_position_distribution <- function(all_data, output_dir) {
    cat("Generating position distribution plots...\n")

    # 1. Genomic position distribution - all events
    pos_counts_all <- all_data %>%
        group_by(comparison, size_category, genomic_position) %>%
        summarize(count = n(), .groups = "drop")

    p1 <- ggplot(pos_counts_all, aes(x = size_category, y = count, fill = genomic_position)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = position_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Exon Size Category",
            y = "Number of Events",
            title = "Genomic Position Distribution by Exon Size",
            subtitle = "All skipped exon events",
            fill = "Position"
        )

    ggsave(file.path(output_dir, "position_distribution_all.pdf"), p1, width = 14, height = 10)

    # 2. Genomic position - significant only
    pos_counts_sig <- all_data %>%
        filter(significant) %>%
        group_by(comparison, size_category, genomic_position) %>%
        summarize(count = n(), .groups = "drop")

    p2 <- ggplot(pos_counts_sig, aes(x = size_category, y = count, fill = genomic_position)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = position_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Exon Size Category",
            y = "Number of Significant Events",
            title = "Genomic Position Distribution (Significant Events)",
            subtitle = paste0("FDR < ", FDR_THRESHOLD, ", |dPSI| > ", DPSI_THRESHOLD),
            fill = "Position"
        )

    ggsave(file.path(output_dir, "position_distribution_significant.pdf"), p2, width = 14, height = 10)

    # 3. Percentage breakdown for microexons specifically
    microexon_pct <- all_data %>%
        filter(size_category == "Microexon (0-30bp)", significant) %>%
        group_by(comparison, genomic_position) %>%
        summarize(count = n(), .groups = "drop") %>%
        group_by(comparison) %>%
        mutate(percentage = count / sum(count) * 100)

    p3 <- ggplot(microexon_pct, aes(x = comparison, y = percentage, fill = genomic_position)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = position_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Comparison",
            y = "Percentage",
            title = "Genomic Position of Significant Microexons",
            fill = "Position"
        )

    ggsave(file.path(output_dir, "microexon_position_percentage.pdf"), p3, width = 10, height = 8)

    # 4. Relative position histogram
    if ("relative_position" %in% colnames(all_data)) {
        rel_pos_data <- all_data %>%
            filter(!is.na(relative_position), significant)
    
        if (nrow(rel_pos_data) > 0) {
            p4 <- ggplot(rel_pos_data, aes(x = relative_position, fill = size_category)) +
                geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
                facet_wrap(~comparison) +
                scale_fill_brewer(palette = "Set1") +
                theme_bw() +
                labs(
                    x = "Relative Position in Transcript (0 = start, 1 = end)",
                    y = "Count",
                    title = "Distribution of Exon Position Within Transcripts",
                    fill = "Size Category"
                )
    
            ggsave(file.path(output_dir, "relative_position_histogram.pdf"), p4, width = 14, height = 10)
        }
    } else {
        cat("  Skipping relative position plot (column not found)\n")
    }
}

# =============================================================================
# Function: Plot biotype distribution
# =============================================================================
plot_biotype_distribution <- function(all_data, output_dir) {
    cat("Generating biotype distribution plots...\n")

    # 1. Gene biotype - all events
    biotype_counts <- all_data %>%
        group_by(comparison, size_category, gene_biotype_simple) %>%
        summarize(count = n(), .groups = "drop")

    p1 <- ggplot(biotype_counts, aes(x = size_category, y = count, fill = gene_biotype_simple)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = biotype_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Exon Size Category",
            y = "Number of Events",
            title = "Gene Biotype Distribution by Exon Size",
            fill = "Gene Biotype"
        )

    ggsave(file.path(output_dir, "biotype_distribution_all.pdf"), p1, width = 14, height = 10)

    # 2. Significant events only
    biotype_sig <- all_data %>%
        filter(significant) %>%
        group_by(comparison, size_category, gene_biotype_simple) %>%
        summarize(count = n(), .groups = "drop")

    p2 <- ggplot(biotype_sig, aes(x = size_category, y = count, fill = gene_biotype_simple)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~comparison) +
        scale_fill_manual(values = biotype_colors) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
            x = "Exon Size Category",
            y = "Number of Significant Events",
            title = "Gene Biotype Distribution (Significant Events)",
            fill = "Gene Biotype"
        )

    ggsave(file.path(output_dir, "biotype_distribution_significant.pdf"), p2, width = 14, height = 10)
}

# =============================================================================
# Function: Summarize functional annotation
# =============================================================================
summarize_functional_annotation <- function(all_data, output_dir) {
    cat("Generating summary statistics...\n")

    # 1. Overall summary by position
    position_summary <- all_data %>%
        filter(significant) %>%
        group_by(comparison, size_category, genomic_position) %>%
        summarize(
            n_events = n(),
            mean_dPSI = mean(IncLevelDifference),
            median_dPSI = median(IncLevelDifference),
            pct_inclusion = sum(IncLevelDifference > 0) / n() * 100,
            .groups = "drop"
        )

    write.csv(position_summary, file.path(output_dir, "position_summary_stats.csv"), row.names = FALSE)

    # 2. Biotype summary
    biotype_summary <- all_data %>%
        filter(significant) %>%
        group_by(comparison, size_category, gene_biotype_simple) %>%
        summarize(
            n_events = n(),
            n_genes = n_distinct(GeneID),
            mean_dPSI = mean(IncLevelDifference),
            .groups = "drop"
        )

    write.csv(biotype_summary, file.path(output_dir, "biotype_summary_stats.csv"), row.names = FALSE)

    # 3. Cross-tabulation: position x biotype for microexons
    microexon_cross <- all_data %>%
        filter(significant, size_category == "Microexon (0-30bp)") %>%
        group_by(comparison, genomic_position, gene_biotype_simple) %>%
        summarize(n_events = n(), .groups = "drop") %>%
        pivot_wider(names_from = gene_biotype_simple, values_from = n_events, values_fill = 0)

    write.csv(microexon_cross, file.path(output_dir, "microexon_position_biotype_cross.csv"), row.names = FALSE)

    # 4. Detailed annotation for all significant events
    annotated_events <- all_data %>%
        filter(significant) %>%
        select(comparison, GeneID, geneSymbol, chr, exonStart_0base, exonEnd, exon_length,
               size_category, genomic_position, gene_biotype_simple,
               IncLevelDifference, FDR) %>%
        arrange(comparison, FDR)

    write.csv(annotated_events, file.path(output_dir, "microexon_genomic_annotation.csv"), row.names = FALSE)

    return(list(position = position_summary, biotype = biotype_summary))
}

# =============================================================================
# Function: Enrichment test for position
# =============================================================================
test_position_enrichment <- function(all_data, output_dir) {
    cat("Testing position enrichment...\n")

    results <- list()

    for (comp in unique(all_data$comparison)) {
        comp_data <- all_data %>% filter(comparison == comp)

        # Compare microexon position distribution to regular exons
        micro_pos <- table(comp_data$genomic_position[comp_data$size_category == "Microexon (0-30bp)"])
        regular_pos <- table(comp_data$genomic_position[comp_data$size_category == "Regular (>50bp)"])

        # Ensure same categories
        all_pos <- union(names(micro_pos), names(regular_pos))
        micro_pos <- micro_pos[all_pos]
        regular_pos <- regular_pos[all_pos]
        micro_pos[is.na(micro_pos)] <- 0
        regular_pos[is.na(regular_pos)] <- 0

        # Chi-square test
        if (sum(micro_pos) > 0 && sum(regular_pos) > 0) {
            test_mat <- rbind(micro_pos, regular_pos)
            tryCatch({
                chi_test <- chisq.test(test_mat)
                results[[comp]] <- data.frame(
                    comparison = comp,
                    chi_sq = chi_test$statistic,
                    p_value = chi_test$p.value,
                    df = chi_test$parameter
                )
            }, error = function(e) {
                results[[comp]] <- data.frame(
                    comparison = comp,
                    chi_sq = NA,
                    p_value = NA,
                    df = NA
                )
            })
        }
    }

    if (length(results) > 0) {
        enrichment_df <- do.call(rbind, results)
        write.csv(enrichment_df, file.path(output_dir, "position_enrichment_test.csv"), row.names = FALSE)
    }
}

# =============================================================================
# Main Execution
# =============================================================================

# 1. Load GTF annotations
cat("\n=== Step 1: Loading GTF Annotations ===\n")
gtf_annotations <- tryCatch({
    load_gtf_annotations(GTF_FILE)
}, error = function(e) {
    warning(paste("Could not load GTF:", e$message))
    warning("Proceeding without genomic position annotation")
    NULL
})

# 2. Load SE data for all comparisons
cat("\n=== Step 2: Loading SE Events ===\n")
all_data <- list()
for (comp in parental_comparisons) {
    cat(paste0("Loading ", comp, "...\n"))
    df <- load_se_with_size(comp, SPLICING_DIR)
    if (!is.null(df)) {
        all_data[[comp]] <- df
    }
}
combined_data <- do.call(rbind, all_data)
cat(paste0("Total SE events loaded: ", nrow(combined_data), "\n"))

# 3. Annotate genomic positions
cat("\n=== Step 3: Annotating Genomic Positions ===\n")
if (!is.null(gtf_annotations)) {
    combined_data <- annotate_exon_position(combined_data, gtf_annotations)
    combined_data <- get_gene_biotype(combined_data, gtf_annotations)
    # Note: relative position calculation can be slow for large datasets
    # combined_data <- calculate_relative_position(combined_data, gtf_annotations)
} else {
    combined_data$genomic_position <- "Unknown"
    combined_data$gene_biotype <- "Unknown"
    combined_data$gene_biotype_simple <- "Unknown"
}

# 4. Generate plots
cat("\n=== Step 4: Generating Plots ===\n")
plot_position_distribution(combined_data, OUTPUT_DIR)
plot_biotype_distribution(combined_data, OUTPUT_DIR)

# 5. Generate summaries
cat("\n=== Step 5: Generating Summaries ===\n")
summaries <- summarize_functional_annotation(combined_data, OUTPUT_DIR)

# 6. Enrichment tests
cat("\n=== Step 6: Enrichment Tests ===\n")
test_position_enrichment(combined_data, OUTPUT_DIR)

# =============================================================================
# Final Summary
# =============================================================================

cat("\n============================================\n")
cat("Functional Annotation Summary\n")
cat("============================================\n")

for (comp in parental_comparisons) {
    comp_data <- combined_data %>% filter(comparison == comp, significant)

    cat(paste0("\n", comp, ":\n"))
    cat(paste0("  Total significant SE: ", nrow(comp_data), "\n"))

    # Microexon position breakdown
    micro_pos <- comp_data %>%
        filter(size_category == "Microexon (0-30bp)") %>%
        group_by(genomic_position) %>%
        summarize(n = n(), .groups = "drop")

    if (nrow(micro_pos) > 0) {
        cat("  Microexon positions:\n")
        for (i in seq_len(nrow(micro_pos))) {
            cat(paste0("    ", micro_pos$genomic_position[i], ": ", micro_pos$n[i], "\n"))
        }
    }
}

cat("\n============================================\n")
cat("Functional Annotation Complete!\n")
cat("============================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - microexon_genomic_annotation.csv: Full annotation table\n")
cat("  - position_summary_stats.csv: Position statistics\n")
cat("  - biotype_summary_stats.csv: Gene biotype statistics\n")
cat("  - position_distribution_*.pdf: Position plots\n")
cat("  - biotype_distribution_*.pdf: Biotype plots\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
