#!/usr/bin/env Rscript
# =============================================================================
# GO/GSEA Enrichment Analysis for Differentially Spliced Genes
# Project: 90-1239779069
# Analysis: Gene Ontology and Gene Set Enrichment for DS genes
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
})

# Configuration
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Elisa_top/90-1239779069"
SPLICING_DIR <- file.path(BASE_DIR, "results/05_splicing")
DESEQ_DIR <- file.path(BASE_DIR, "results/04_deseq2")
OUTPUT_DIR <- file.path(BASE_DIR, "results/08_enrichment_analysis")

# Thresholds
FDR_THRESHOLD <- 0.05
DPSI_THRESHOLD <- 0.1
PADJ_THRESHOLD <- 0.05  # For gene expression
LFC_THRESHOLD <- 1      # log2FC threshold for DE genes

cat("============================================\n")
cat("GO/GSEA Enrichment Analysis\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Splicing thresholds: FDR <", FDR_THRESHOLD, ", |dPSI| >", DPSI_THRESHOLD, "\n")
cat("Expression thresholds: padj <", PADJ_THRESHOLD, ", |log2FC| >", LFC_THRESHOLD, "\n")
cat("============================================\n\n")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "splicing_GO"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "splicing_GSEA"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "expression_GO"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "expression_GSEA"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "comparison"), showWarnings = FALSE)

# Define comparisons
parental_comparisons <- c("Neg_vs_Parental", "Pos_vs_Parental", "KO_vs_Parental")
event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

# =============================================================================
# Function: Load rMATS results
# =============================================================================
load_rmats_results <- function(comparison, event_type, splicing_dir) {
    file_path <- file.path(splicing_dir, comparison, paste0(event_type, ".MATS.JC.txt"))
    if (!file.exists(file_path)) {
        return(NULL)
    }
    df <- read.delim(file_path, stringsAsFactors = FALSE)
    df$event_type <- event_type
    df$comparison <- comparison
    
    # Calculate average coverage (IJC + SJC)
    # Define parsing function inside or check if it exists globally
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
# Function: Get differentially spliced genes
# =============================================================================
get_ds_genes <- function(splicing_dir, comparison, event_types) {
    cat(paste0("Extracting DS genes for ", comparison, "...\n"))

    all_sig <- list()

    for (event in event_types) {
        df <- load_rmats_results(comparison, event, splicing_dir)
        if (!is.null(df) && nrow(df) > 0) {
            sig <- df %>%
                filter(FDR < FDR_THRESHOLD, abs(IncLevelDifference) > DPSI_THRESHOLD) %>%
                select(GeneID, geneSymbol, IncLevelDifference, FDR, event_type)
            if (nrow(sig) > 0) {
                all_sig[[event]] <- sig
            }
        }
    }

    if (length(all_sig) == 0) {
        return(NULL)
    }

    combined <- do.call(rbind, all_sig)

    # Get unique genes with best dPSI
    gene_summary <- combined %>%
        group_by(GeneID, geneSymbol) %>%
        summarize(
            max_abs_dPSI = max(abs(IncLevelDifference)),
            best_dPSI = IncLevelDifference[which.max(abs(IncLevelDifference))],
            min_FDR = min(FDR),
            n_events = n(),
            event_types = paste(unique(event_type), collapse = ","),
            .groups = "drop"
        ) %>%
        arrange(min_FDR)

    cat(paste0("  Found ", nrow(gene_summary), " unique DS genes\n"))

    return(gene_summary)
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
# Function: Run GO enrichment
# =============================================================================
run_go_enrichment <- function(gene_ids, comparison, output_dir, prefix = "splicing") {
    cat(paste0("Running GO enrichment for ", comparison, "...\n"))

    # Convert to Entrez IDs
    entrez_ids <- convert_to_entrez(gene_ids)
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]

    if (length(entrez_ids) < 5) {
        warning(paste("Not enough genes with Entrez IDs for", comparison))
        return(NULL)
    }

    cat(paste0("  ", length(entrez_ids), " genes with Entrez IDs\n"))

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
                 file.path(output_dir, paste0(prefix, "_GO"),
                          paste0(comparison, "_GO_BP.csv")),
                 row.names = FALSE)

        # Dot plot
        if (nrow(as.data.frame(go_bp)) >= 5) {
            p <- dotplot(go_bp, showCategory = 20,
                        title = paste0("GO Biological Process: ", comparison))
            ggsave(file.path(output_dir, paste0(prefix, "_GO"),
                            paste0(comparison, "_GO_BP_dotplot.pdf")),
                   p, width = 10, height = 12)

            # Enrichment map
            tryCatch({
                go_bp_sim <- pairwise_termsim(go_bp)
                p2 <- emapplot(go_bp_sim, showCategory = 30)
                ggsave(file.path(output_dir, paste0(prefix, "_GO"),
                                paste0(comparison, "_GO_BP_emap.pdf")),
                       p2, width = 12, height = 12)
            }, error = function(e) {
                warning(paste("Could not create emapplot for", comparison))
            })
        }
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
                 file.path(output_dir, paste0(prefix, "_GO"),
                          paste0(comparison, "_GO_MF.csv")),
                 row.names = FALSE)

        if (nrow(as.data.frame(go_mf)) >= 5) {
            p <- dotplot(go_mf, showCategory = 20,
                        title = paste0("GO Molecular Function: ", comparison))
            ggsave(file.path(output_dir, paste0(prefix, "_GO"),
                            paste0(comparison, "_GO_MF_dotplot.pdf")),
                   p, width = 10, height = 12)
        }
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
                 file.path(output_dir, paste0(prefix, "_GO"),
                          paste0(comparison, "_GO_CC.csv")),
                 row.names = FALSE)

        if (nrow(as.data.frame(go_cc)) >= 5) {
            p <- dotplot(go_cc, showCategory = 20,
                        title = paste0("GO Cellular Component: ", comparison))
            ggsave(file.path(output_dir, paste0(prefix, "_GO"),
                            paste0(comparison, "_GO_CC_dotplot.pdf")),
                   p, width = 10, height = 12)
        }
    }

    return(results)
}

# =============================================================================
# Function: Run GSEA ranked by dPSI
# =============================================================================
run_gsea_dpsi <- function(gene_data, comparison, output_dir) {
    cat(paste0("Running GSEA for ", comparison, "...\n"))

    # Convert to Entrez and create ranked list
    gene_data$entrez <- convert_to_entrez(gene_data$GeneID)
    gene_data <- gene_data[!is.na(gene_data$entrez), ]

    if (nrow(gene_data) < 50) {
        warning(paste("Not enough genes for GSEA in", comparison))
        return(NULL)
    }

    # Create ranked gene list (by dPSI)
    gene_list <- gene_data$best_dPSI
    names(gene_list) <- gene_data$entrez
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Remove duplicates (keep first/highest)
    gene_list <- gene_list[!duplicated(names(gene_list))]

    cat(paste0("  ", length(gene_list), " genes for GSEA\n"))

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
                     file.path(output_dir, "splicing_GSEA",
                              paste0(comparison, "_GSEA_GO_BP.csv")),
                     row.names = FALSE)

            # GSEA plot for top pathways
            if (nrow(as.data.frame(gsea_result)) >= 5) {
                p <- dotplot(gsea_result, showCategory = 20,
                            title = paste0("GSEA GO BP: ", comparison),
                            split = ".sign") +
                    facet_grid(.~.sign)
                ggsave(file.path(output_dir, "splicing_GSEA",
                                paste0(comparison, "_GSEA_dotplot.pdf")),
                       p, width = 14, height = 12)

                # Ridge plot
                p2 <- ridgeplot(gsea_result, showCategory = 15) +
                    labs(title = paste0("GSEA Ridge Plot: ", comparison))
                ggsave(file.path(output_dir, "splicing_GSEA",
                                paste0(comparison, "_GSEA_ridgeplot.pdf")),
                       p2, width = 12, height = 10)
            }

            return(gsea_result)
        }
    }, error = function(e) {
        warning(paste("GSEA failed for", comparison, ":", e$message))
    })

    return(NULL)
}

# =============================================================================
# Function: Run GO/GSEA for gene expression
# =============================================================================
run_expression_enrichment <- function(deseq_dir, comparison, output_dir) {
    cat(paste0("Running enrichment for gene expression: ", comparison, "...\n"))

    # Load DESeq2 results
    file_path <- file.path(deseq_dir, paste0(comparison, "_all.csv"))
    if (!file.exists(file_path)) {
        warning(paste("DESeq2 results not found:", file_path))
        return(NULL)
    }

    de_results <- read.csv(file_path)

    # Get significant DE genes
    sig_genes <- de_results %>%
        filter(!is.na(padj), padj < PADJ_THRESHOLD, abs(log2FoldChange) > LFC_THRESHOLD)

    cat(paste0("  ", nrow(sig_genes), " significant DE genes\n"))

    if (nrow(sig_genes) < 5) {
        warning(paste("Not enough DE genes for", comparison))
        return(NULL)
    }

    # Convert to Entrez
    entrez_ids <- convert_to_entrez(sig_genes$gene_id)
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]

    if (length(entrez_ids) < 5) {
        return(NULL)
    }

    # GO enrichment
    go_bp <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)

    if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
        write.csv(as.data.frame(go_bp),
                 file.path(output_dir, "expression_GO",
                          paste0(comparison, "_DE_GO_BP.csv")),
                 row.names = FALSE)

        if (nrow(as.data.frame(go_bp)) >= 5) {
            p <- dotplot(go_bp, showCategory = 20,
                        title = paste0("GO BP (DE genes): ", comparison))
            ggsave(file.path(output_dir, "expression_GO",
                            paste0(comparison, "_DE_GO_BP_dotplot.pdf")),
                   p, width = 10, height = 12)
        }
    }

    # GSEA with all genes
    de_results$entrez <- convert_to_entrez(de_results$gene_id)
    de_ranked <- de_results %>%
        filter(!is.na(entrez), !is.na(log2FoldChange)) %>%
        arrange(desc(log2FoldChange))

    gene_list <- de_ranked$log2FoldChange
    names(gene_list) <- de_ranked$entrez
    gene_list <- gene_list[!duplicated(names(gene_list))]

    tryCatch({
        gsea_result <- gseGO(geneList = gene_list,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP",
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 0.05,
                            verbose = FALSE)

        if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
            write.csv(as.data.frame(gsea_result),
                     file.path(output_dir, "expression_GSEA",
                              paste0(comparison, "_DE_GSEA.csv")),
                     row.names = FALSE)

            if (nrow(as.data.frame(gsea_result)) >= 5) {
                p <- dotplot(gsea_result, showCategory = 20,
                            title = paste0("GSEA (DE genes): ", comparison),
                            split = ".sign") +
                    facet_grid(.~.sign)
                ggsave(file.path(output_dir, "expression_GSEA",
                                paste0(comparison, "_DE_GSEA_dotplot.pdf")),
                       p, width = 14, height = 12)
            }
        }
    }, error = function(e) {
        warning(paste("Expression GSEA failed for", comparison))
    })

    return(list(GO = go_bp))
}

# =============================================================================
# Function: Compare splicing vs expression enrichment
# =============================================================================
compare_enrichment <- function(splicing_go, expression_go, comparison, output_dir) {
    cat(paste0("Comparing splicing vs expression enrichment: ", comparison, "...\n"))

    if (is.null(splicing_go$BP) || is.null(expression_go$GO)) {
        return(NULL)
    }

    # Get top terms from each
    splicing_terms <- head(as.data.frame(splicing_go$BP), 30)$Description
    expression_terms <- head(as.data.frame(expression_go$GO), 30)$Description

    # Find overlap
    common_terms <- intersect(splicing_terms, expression_terms)

    if (length(common_terms) > 0) {
        cat(paste0("  Common GO terms: ", length(common_terms), "\n"))

        # Create comparison data frame
        comparison_df <- data.frame(
            Term = c(splicing_terms, expression_terms),
            Source = c(rep("Splicing", length(splicing_terms)),
                      rep("Expression", length(expression_terms)))
        )

        write.csv(comparison_df,
                 file.path(output_dir, "comparison",
                          paste0(comparison, "_enrichment_comparison.csv")),
                 row.names = FALSE)
    }

    return(common_terms)
}

# =============================================================================
# Main Execution
# =============================================================================

all_splicing_go <- list()
all_splicing_gsea <- list()
all_expression_go <- list()

for (comp in parental_comparisons) {
    cat(paste0("\n=== Processing: ", comp, " ===\n"))

    # 1. Get DS genes
    ds_genes <- get_ds_genes(SPLICING_DIR, comp, event_types)

    if (!is.null(ds_genes) && nrow(ds_genes) > 0) {
        # Save DS genes list
        write.csv(ds_genes,
                 file.path(OUTPUT_DIR, paste0(comp, "_DS_genes.csv")),
                 row.names = FALSE)

        # 2. GO enrichment for DS genes
        splicing_go <- run_go_enrichment(ds_genes$GeneID, comp, OUTPUT_DIR, "splicing")
        all_splicing_go[[comp]] <- splicing_go

        # 3. GSEA ranked by dPSI
        splicing_gsea <- run_gsea_dpsi(ds_genes, comp, OUTPUT_DIR)
        all_splicing_gsea[[comp]] <- splicing_gsea
    }

    # 4. Expression enrichment
    expression_go <- run_expression_enrichment(DESEQ_DIR, comp, OUTPUT_DIR)
    all_expression_go[[comp]] <- expression_go

    # 5. Compare splicing vs expression
    if (!is.null(all_splicing_go[[comp]]) && !is.null(expression_go)) {
        compare_enrichment(all_splicing_go[[comp]], expression_go, comp, OUTPUT_DIR)
    }
}

# =============================================================================
# Combined Analysis
# =============================================================================

cat("\n=== Combined Analysis ===\n")

# Get all DS genes across comparisons
all_ds_genes <- list()
for (comp in parental_comparisons) {
    ds_genes <- get_ds_genes(SPLICING_DIR, comp, event_types)
    if (!is.null(ds_genes)) {
        all_ds_genes[[comp]] <- ds_genes
    }
}

# Find common DS genes
if (length(all_ds_genes) >= 2) {
    common_genes <- Reduce(intersect, lapply(all_ds_genes, function(x) x$GeneID))
    cat(paste0("Common DS genes across all comparisons: ", length(common_genes), "\n"))

    if (length(common_genes) >= 5) {
        # Run GO on common DS genes
        common_go <- run_go_enrichment(common_genes, "common_DS_genes", OUTPUT_DIR, "splicing")

        # Save common DS genes
        common_df <- all_ds_genes[[1]] %>%
            filter(GeneID %in% common_genes)
        write.csv(common_df,
                 file.path(OUTPUT_DIR, "common_DS_genes.csv"),
                 row.names = FALSE)
    }
}

# =============================================================================
# Summary
# =============================================================================

cat("\n============================================\n")
cat("Enrichment Analysis Summary\n")
cat("============================================\n")

for (comp in parental_comparisons) {
    cat(paste0("\n", comp, ":\n"))

    # DS genes count
    if (comp %in% names(all_ds_genes)) {
        cat(paste0("  DS genes: ", nrow(all_ds_genes[[comp]]), "\n"))
    }

    # Splicing GO
    if (comp %in% names(all_splicing_go) && !is.null(all_splicing_go[[comp]]$BP)) {
        cat(paste0("  Splicing GO BP terms: ",
                  nrow(as.data.frame(all_splicing_go[[comp]]$BP)), "\n"))
    }

    # Expression GO
    if (comp %in% names(all_expression_go) && !is.null(all_expression_go[[comp]]$GO)) {
        cat(paste0("  Expression GO BP terms: ",
                  nrow(as.data.frame(all_expression_go[[comp]]$GO)), "\n"))
    }
}

cat("\n============================================\n")
cat("Enrichment Analysis Complete!\n")
cat("============================================\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("\nGenerated outputs:\n")
cat("  - splicing_GO/: GO enrichment for differentially spliced genes\n")
cat("  - splicing_GSEA/: GSEA ranked by dPSI\n")
cat("  - expression_GO/: GO enrichment for DE genes\n")
cat("  - expression_GSEA/: GSEA for gene expression\n")
cat("  - comparison/: Splicing vs expression comparison\n")
cat("\nEnd time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================\n")
