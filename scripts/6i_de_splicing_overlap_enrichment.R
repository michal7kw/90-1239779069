#!/usr/bin/env Rscript
# =============================================================================
# GO/GSEA Enrichment for Genes with BOTH DE and Differential Splicing
# Project: 90-1239779069
# Analysis: Enrichment on genes significant in BOTH differential expression
#           AND alternative splicing (overlap genes)
# =============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(DOSE)
    library(VennDiagram)
    library(grid)
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
DESEQ_DIR <- file.path(BASE_DIR, "results/04_deseq2")
EVENTS_DIR <- file.path(BASE_DIR, "results/10_significant_events")
OUTPUT_DIR <- file.path(BASE_DIR, "results/13_de_splicing_overlap_enrichment")

# Thresholds (DE significance already filtered in input files)
# Splicing events already filtered: FDR < 0.05, |dPSI| > 0.1
MIN_GENES_GO <- 5
MIN_GENES_GSEA <- 15

cat("============================================\n")
cat("DE + Splicing Overlap Enrichment Analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "gene_lists"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "overlap_GO"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "overlap_GSEA"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "comparison_plots"), showWarnings = FALSE)

# Define comparisons (only vs Parental)
comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")

# =============================================================================
# Function: Load DE significant genes
# =============================================================================
load_de_significant <- function(comparison, deseq_dir) {
    file_path <- file.path(deseq_dir, paste0(comparison, "_significant.csv"))
    if (!file.exists(file_path)) {
        warning(paste("DE file not found:", file_path))
        return(NULL)
    }

    df <- read.csv(file_path, stringsAsFactors = FALSE)
    # Clean gene IDs (remove version numbers if present)
    df$gene_id_clean <- gsub("\\..*", "", df$gene_id)

    cat(paste0("  Loaded ", nrow(df), " significant DE genes\n"))
    return(df)
}

# =============================================================================
# Function: Load DE all genes (for GSEA ranking)
# =============================================================================
load_de_all <- function(comparison, deseq_dir) {
    file_path <- file.path(deseq_dir, paste0(comparison, "_all.csv"))
    if (!file.exists(file_path)) {
        warning(paste("DE all file not found:", file_path))
        return(NULL)
    }

    df <- read.csv(file_path, stringsAsFactors = FALSE)
    df$gene_id_clean <- gsub("\\..*", "", df$gene_id)

    return(df)
}

# =============================================================================
# Function: Load splicing significant events
# =============================================================================
load_splicing_significant <- function(comparison, events_dir) {
    file_path <- file.path(events_dir, paste0(comparison, "_significant_events.csv"))
    if (!file.exists(file_path)) {
        warning(paste("Splicing file not found:", file_path))
        return(NULL)
    }

    df <- read.csv(file_path, stringsAsFactors = FALSE)
    # Clean gene IDs
    df$GeneID_clean <- gsub("\\..*", "", df$GeneID)

    # Get unique genes with summary statistics
    gene_summary <- df %>%
        group_by(GeneID_clean, GeneSymbol) %>%
        summarize(
            n_events = n(),
            event_types = paste(unique(EventType), collapse = ","),
            max_abs_dPSI = max(abs(DeltaPSI)),
            best_dPSI = DeltaPSI[which.max(abs(DeltaPSI))],
            min_FDR = min(FDR),
            .groups = "drop"
        )

    cat(paste0("  Loaded ", nrow(df), " splicing events from ",
               nrow(gene_summary), " unique genes\n"))
    return(gene_summary)
}

# =============================================================================
# Function: Identify overlap genes
# =============================================================================
identify_overlap_genes <- function(de_df, splicing_df, comparison) {
    cat(paste0("  Identifying overlap genes...\n"))

    # Inner join on cleaned gene IDs
    overlap <- inner_join(
        de_df %>% select(gene_id_clean, gene_symbol, log2FoldChange, padj, baseMean),
        splicing_df %>% select(GeneID_clean, GeneSymbol, n_events, event_types,
                               max_abs_dPSI, best_dPSI, min_FDR),
        by = c("gene_id_clean" = "GeneID_clean")
    )

    # Rename columns for clarity
    overlap <- overlap %>%
        rename(
            gene_id = gene_id_clean,
            DE_log2FC = log2FoldChange,
            DE_padj = padj,
            DE_baseMean = baseMean,
            splicing_n_events = n_events,
            splicing_event_types = event_types,
            splicing_max_abs_dPSI = max_abs_dPSI,
            splicing_best_dPSI = best_dPSI,
            splicing_min_FDR = min_FDR
        )

    cat(paste0("  Found ", nrow(overlap), " overlap genes\n"))
    return(overlap)
}

# =============================================================================
# Function: Convert gene IDs to Entrez
# =============================================================================
convert_to_entrez <- function(ensembl_ids) {
    # Remove version numbers
    clean_ids <- gsub("\\..*", "", ensembl_ids)

    # Map to Entrez
    entrez_ids <- mapIds(org.Mm.eg.db,
                         keys = clean_ids,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

    return(entrez_ids)
}

# =============================================================================
# Function: Run GO enrichment for overlap genes
# =============================================================================
run_overlap_go_enrichment <- function(gene_ids, comparison, output_dir) {
    cat(paste0("  Running GO enrichment...\n"))

    # Convert to Entrez IDs
    entrez_ids <- convert_to_entrez(gene_ids)
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]

    n_mapped <- length(entrez_ids)
    n_input <- length(gene_ids)
    cat(paste0("    Mapped ", n_mapped, " of ", n_input, " genes to Entrez IDs\n"))

    if (n_mapped < MIN_GENES_GO) {
        warning(paste("Not enough genes with Entrez IDs for GO enrichment:", n_mapped))
        return(NULL)
    }

    results <- list()

    # GO Biological Process
    go_bp <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
        results$BP <- go_bp

        # Save results
        write.csv(as.data.frame(go_bp),
                 file.path(output_dir, "overlap_GO",
                          paste0(comparison, "_GO_BP.csv")),
                 row.names = FALSE)

        # Dot plot
        if (nrow(as.data.frame(go_bp)) >= 5) {
            p <- dotplot(go_bp, showCategory = 20,
                        title = paste0("GO BP (DE+Splicing Overlap): ", comparison))
            ggsave(file.path(output_dir, "overlap_GO",
                            paste0(comparison, "_GO_BP_dotplot.pdf")),
                   p, width = 10, height = 12)

            # Enrichment map
            tryCatch({
                go_bp_sim <- pairwise_termsim(go_bp)
                p2 <- emapplot(go_bp_sim, showCategory = 30)
                ggsave(file.path(output_dir, "overlap_GO",
                                paste0(comparison, "_GO_BP_emap.pdf")),
                       p2, width = 12, height = 12)
            }, error = function(e) {
                warning(paste("Could not create emapplot for", comparison))
            })
        }

        cat(paste0("    GO BP: ", nrow(as.data.frame(go_bp)), " significant terms\n"))
    } else {
        cat(paste0("    GO BP: No significant terms found\n"))
    }

    # GO Molecular Function
    go_mf <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Mm.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
        results$MF <- go_mf

        write.csv(as.data.frame(go_mf),
                 file.path(output_dir, "overlap_GO",
                          paste0(comparison, "_GO_MF.csv")),
                 row.names = FALSE)

        if (nrow(as.data.frame(go_mf)) >= 5) {
            p <- dotplot(go_mf, showCategory = 20,
                        title = paste0("GO MF (DE+Splicing Overlap): ", comparison))
            ggsave(file.path(output_dir, "overlap_GO",
                            paste0(comparison, "_GO_MF_dotplot.pdf")),
                   p, width = 10, height = 12)
        }

        cat(paste0("    GO MF: ", nrow(as.data.frame(go_mf)), " significant terms\n"))
    }

    # GO Cellular Component
    go_cc <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Mm.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
        results$CC <- go_cc

        write.csv(as.data.frame(go_cc),
                 file.path(output_dir, "overlap_GO",
                          paste0(comparison, "_GO_CC.csv")),
                 row.names = FALSE)

        if (nrow(as.data.frame(go_cc)) >= 5) {
            p <- dotplot(go_cc, showCategory = 20,
                        title = paste0("GO CC (DE+Splicing Overlap): ", comparison))
            ggsave(file.path(output_dir, "overlap_GO",
                            paste0(comparison, "_GO_CC_dotplot.pdf")),
                   p, width = 10, height = 12)
        }

        cat(paste0("    GO CC: ", nrow(as.data.frame(go_cc)), " significant terms\n"))
    }

    return(results)
}

# =============================================================================
# Function: Run GSEA ranked by log2FoldChange
# =============================================================================
run_overlap_gsea <- function(overlap_genes, de_all_df, comparison, output_dir) {
    cat(paste0("  Running GSEA (ranked by log2FC)...\n"))

    # Get overlap gene IDs
    overlap_gene_ids <- overlap_genes$gene_id

    if (length(overlap_gene_ids) < MIN_GENES_GSEA) {
        warning(paste("Not enough overlap genes for GSEA:", length(overlap_gene_ids)))
        return(NULL)
    }

    # Filter DE all results to overlap genes only
    de_overlap <- de_all_df %>%
        filter(gene_id_clean %in% overlap_gene_ids) %>%
        filter(!is.na(log2FoldChange))

    # Convert to Entrez
    de_overlap$entrez <- convert_to_entrez(de_overlap$gene_id_clean)
    de_overlap <- de_overlap %>%
        filter(!is.na(entrez))

    if (nrow(de_overlap) < MIN_GENES_GSEA) {
        warning(paste("Not enough genes with Entrez IDs for GSEA:", nrow(de_overlap)))
        return(NULL)
    }

    # Create ranked gene list by log2FoldChange
    gene_list <- de_overlap$log2FoldChange
    names(gene_list) <- de_overlap$entrez
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Remove duplicates
    gene_list <- gene_list[!duplicated(names(gene_list))]

    cat(paste0("    ", length(gene_list), " genes for GSEA\n"))

    # Run GSEA with GO BP
    tryCatch({
        gsea_result <- gseGO(geneList = gene_list,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP",
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 0.05,
                            verbose = FALSE)

        if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
            # Save results
            write.csv(as.data.frame(gsea_result),
                     file.path(output_dir, "overlap_GSEA",
                              paste0(comparison, "_GSEA_GO_BP.csv")),
                     row.names = FALSE)

            cat(paste0("    GSEA: ", nrow(as.data.frame(gsea_result)), " significant terms\n"))

            # GSEA dotplot
            if (nrow(as.data.frame(gsea_result)) >= 5) {
                p <- dotplot(gsea_result, showCategory = 20,
                            title = paste0("GSEA (DE+Splicing Overlap): ", comparison),
                            split = ".sign") +
                    facet_grid(.~.sign)
                ggsave(file.path(output_dir, "overlap_GSEA",
                                paste0(comparison, "_GSEA_dotplot.pdf")),
                       p, width = 14, height = 12)

                # Ridge plot
                p2 <- ridgeplot(gsea_result, showCategory = 15) +
                    labs(title = paste0("GSEA Ridge Plot (Overlap): ", comparison))
                ggsave(file.path(output_dir, "overlap_GSEA",
                                paste0(comparison, "_GSEA_ridgeplot.pdf")),
                       p2, width = 12, height = 10)
            }

            return(gsea_result)
        } else {
            cat(paste0("    GSEA: No significant terms found\n"))
        }
    }, error = function(e) {
        warning(paste("GSEA failed for", comparison, ":", e$message))
    })

    return(NULL)
}

# =============================================================================
# Function: Generate Venn diagram
# =============================================================================
generate_venn_diagram <- function(de_genes, splicing_genes, comparison, output_dir) {
    cat(paste0("  Generating Venn diagram...\n"))

    de_set <- unique(de_genes$gene_id_clean)
    splicing_set <- unique(splicing_genes$GeneID_clean)

    # Create Venn diagram
    pdf(file.path(output_dir, "comparison_plots",
                  paste0(comparison, "_venn_overlap.pdf")),
        width = 8, height = 8)

    venn <- draw.pairwise.venn(
        area1 = length(de_set),
        area2 = length(splicing_set),
        cross.area = length(intersect(de_set, splicing_set)),
        category = c("DE genes", "Splicing genes"),
        fill = c("#E41A1C", "#377EB8"),
        alpha = 0.5,
        cat.pos = c(-20, 20),
        cat.dist = 0.05,
        cat.cex = 1.2,
        cex = 1.5,
        fontfamily = "sans",
        cat.fontfamily = "sans",
        main = comparison
    )

    grid.text(comparison, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

    dev.off()

    return(list(
        de_only = length(setdiff(de_set, splicing_set)),
        splicing_only = length(setdiff(splicing_set, de_set)),
        overlap = length(intersect(de_set, splicing_set))
    ))
}

# =============================================================================
# Main Execution
# =============================================================================

all_results <- list()
venn_stats <- list()

for (comp in comparisons) {
    cat(paste0("\n=== Processing: ", comp, " ===\n"))

    # 1. Load data
    de_sig <- load_de_significant(comp, DESEQ_DIR)
    de_all <- load_de_all(comp, DESEQ_DIR)
    splicing_sig <- load_splicing_significant(comp, EVENTS_DIR)

    if (is.null(de_sig) || is.null(splicing_sig)) {
        warning(paste("Missing data for", comp, "- skipping"))
        next
    }

    # 2. Generate Venn diagram
    venn_stats[[comp]] <- generate_venn_diagram(de_sig, splicing_sig, comp, OUTPUT_DIR)

    # 3. Identify overlap genes
    overlap <- identify_overlap_genes(de_sig, splicing_sig, comp)

    if (nrow(overlap) == 0) {
        warning(paste("No overlap genes found for", comp))
        next
    }

    # 4. Save gene lists
    # Overlap genes
    write.csv(overlap,
             file.path(OUTPUT_DIR, "gene_lists",
                      paste0(comp, "_overlap_genes.csv")),
             row.names = FALSE)

    # DE-only genes
    de_only <- de_sig %>%
        filter(!gene_id_clean %in% overlap$gene_id) %>%
        select(gene_id_clean, gene_symbol, log2FoldChange, padj)
    write.csv(de_only,
             file.path(OUTPUT_DIR, "gene_lists",
                      paste0(comp, "_de_only_genes.csv")),
             row.names = FALSE)

    # Splicing-only genes
    splicing_only <- splicing_sig %>%
        filter(!GeneID_clean %in% overlap$gene_id) %>%
        select(GeneID_clean, GeneSymbol, n_events, event_types, max_abs_dPSI)
    write.csv(splicing_only,
             file.path(OUTPUT_DIR, "gene_lists",
                      paste0(comp, "_splicing_only_genes.csv")),
             row.names = FALSE)

    # 5. Run GO enrichment
    go_results <- run_overlap_go_enrichment(overlap$gene_id, comp, OUTPUT_DIR)

    # 6. Run GSEA
    gsea_results <- run_overlap_gsea(overlap, de_all, comp, OUTPUT_DIR)

    # Store results
    all_results[[comp]] <- list(
        overlap = overlap,
        go = go_results,
        gsea = gsea_results,
        venn = venn_stats[[comp]]
    )
}

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n=== Generating Summary ===\n")

# Create summary data frame
summary_df <- data.frame(
    Comparison = character(),
    DE_genes = integer(),
    Splicing_genes = integer(),
    Overlap_genes = integer(),
    Overlap_percent_DE = numeric(),
    Overlap_percent_splicing = numeric(),
    GO_BP_terms = integer(),
    GO_MF_terms = integer(),
    GO_CC_terms = integer(),
    GSEA_terms = integer(),
    stringsAsFactors = FALSE
)

for (comp in names(all_results)) {
    res <- all_results[[comp]]

    de_count <- res$venn$de_only + res$venn$overlap
    splicing_count <- res$venn$splicing_only + res$venn$overlap

    go_bp_terms <- if (!is.null(res$go$BP)) nrow(as.data.frame(res$go$BP)) else 0
    go_mf_terms <- if (!is.null(res$go$MF)) nrow(as.data.frame(res$go$MF)) else 0
    go_cc_terms <- if (!is.null(res$go$CC)) nrow(as.data.frame(res$go$CC)) else 0
    gsea_terms <- if (!is.null(res$gsea)) nrow(as.data.frame(res$gsea)) else 0

    summary_df <- rbind(summary_df, data.frame(
        Comparison = comp,
        DE_genes = de_count,
        Splicing_genes = splicing_count,
        Overlap_genes = res$venn$overlap,
        Overlap_percent_DE = round(100 * res$venn$overlap / de_count, 1),
        Overlap_percent_splicing = round(100 * res$venn$overlap / splicing_count, 1),
        GO_BP_terms = go_bp_terms,
        GO_MF_terms = go_mf_terms,
        GO_CC_terms = go_cc_terms,
        GSEA_terms = gsea_terms,
        stringsAsFactors = FALSE
    ))
}

write.csv(summary_df,
         file.path(OUTPUT_DIR, "gene_lists", "all_comparisons_overlap_summary.csv"),
         row.names = FALSE)

# Summary bar plot
if (nrow(summary_df) > 0) {
    plot_df <- summary_df %>%
        select(Comparison, DE_only = DE_genes, Splicing_only = Splicing_genes,
               Overlap = Overlap_genes) %>%
        pivot_longer(cols = -Comparison, names_to = "Category", values_to = "Count")

    # Adjust for actual counts
    for (i in seq_len(nrow(summary_df))) {
        comp <- summary_df$Comparison[i]
        plot_df$Count[plot_df$Comparison == comp & plot_df$Category == "DE_only"] <-
            all_results[[comp]]$venn$de_only
        plot_df$Count[plot_df$Comparison == comp & plot_df$Category == "Splicing_only"] <-
            all_results[[comp]]$venn$splicing_only
    }

    p <- ggplot(plot_df, aes(x = Comparison, y = Count, fill = Category)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("DE_only" = "#E41A1C",
                                     "Splicing_only" = "#377EB8",
                                     "Overlap" = "#4DAF4A")) +
        labs(title = "Gene Counts by Category",
             subtitle = "DE-only, Splicing-only, and Overlap genes",
             x = "Comparison",
             y = "Number of Genes") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(OUTPUT_DIR, "comparison_plots", "summary_overlap_counts.pdf"),
           p, width = 10, height = 6)
}

# =============================================================================
# Final Summary
# =============================================================================

cat("\n============================================\n")
cat("DE + Splicing Overlap Enrichment Complete!\n")
cat("============================================\n\n")

cat("Summary by Comparison:\n")
print(summary_df)

cat("\n\nOutput directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - gene_lists/: Overlap, DE-only, and splicing-only gene lists\n")
cat("  - overlap_GO/: GO enrichment (BP, MF, CC) for overlap genes\n")
cat("  - overlap_GSEA/: GSEA ranked by log2FoldChange\n")
cat("  - comparison_plots/: Venn diagrams and summary plots\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
