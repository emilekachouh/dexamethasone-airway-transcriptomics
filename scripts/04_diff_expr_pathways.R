# 04_diff_expr_pathways.R - Differential expression and pathway analysis
# DESeq2 workflow with pathway enrichment

# Load required libraries
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)

# Create output directories
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

# Load filtered data
cat("Loading filtered airway dataset...\n")
airway_filtered <- readRDS("data/processed/airway_filtered.rds")

cat("=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
cat("Design: Dexamethasone treatment vs control\n")
cat("Method: DESeq2 with negative binomial modeling\n")

# =============================================================================
# DESEQ2 ANALYSIS
# =============================================================================

# Create DESeq2 object with design formula
# We'll test for treatment effect while accounting for cell line differences
dds <- DESeqDataSet(airway_filtered, design = ~ cell + dex)

cat("DESeq2 object created with design: ~ cell + dex\n")
cat("Genes:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Get results for dexamethasone effect
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Summary of results
cat("\nDESeq2 results summary:\n")
summary(res)

# =============================================================================
# QUALITY CONTROL OF DE RESULTS
# =============================================================================

cat("\n=== DIFFERENTIAL EXPRESSION QC ===\n")

# MA plot
png("figures/13_ma_plot.png", width = 8, height = 6, units = "in", res = 300)
plotMA(res, main = "MA Plot: Dexamethasone vs Control")
dev.off()
cat("Saved: figures/13_ma_plot.png\n")

# Dispersion plot
png("figures/14_dispersion_plot.png", width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()
cat("Saved: figures/14_dispersion_plot.png\n")

# =============================================================================
# VOLCANO PLOT
# =============================================================================

cat("\n=== VOLCANO PLOT ===\n")

# Prepare data for volcano plot
volcano_data <- data.frame(
  gene = rownames(res),
  log2FC = res$log2FoldChange,
  log10_pval = -log10(res$pvalue),
  padj = res$padj,
  significant = res$padj < 0.05 & !is.na(res$padj),
  highly_significant = res$padj < 0.01 & !is.na(res$padj) & abs(res$log2FoldChange) > 1
)

# Remove rows with missing values
volcano_data <- volcano_data[complete.cases(volcano_data[, c("log2FC", "log10_pval")]), ]

# Count significant genes
n_sig <- sum(volcano_data$significant, na.rm = TRUE)
n_up <- sum(volcano_data$significant & volcano_data$log2FC > 0, na.rm = TRUE)
n_down <- sum(volcano_data$significant & volcano_data$log2FC < 0, na.rm = TRUE)

cat("Significant genes (FDR < 0.05):", n_sig, "\n")
cat("Upregulated:", n_up, "\n")
cat("Downregulated:", n_down, "\n")

# Create volcano plot
p1 <- ggplot(volcano_data, aes(x = log2FC, y = log10_pval)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Volcano Plot: Dexamethasone vs Control",
       subtitle = paste("Significant genes:", n_sig, "(", n_up, "up,", n_down, "down)"),
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "FDR < 0.05") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/15_volcano_plot.png", p1, width = 10, height = 8)
cat("Saved: figures/15_volcano_plot.png\n")

# =============================================================================
# TOP DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

cat("\n=== TOP DIFFERENTIALLY EXPRESSED GENES ===\n")

# Get top genes by adjusted p-value
res_ordered <- res[order(res$padj), ]
top_genes <- head(res_ordered, 20)

# Add gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = rownames(top_genes),
                       column = "SYMBOL", 
                       keytype = "ENSEMBL",
                       multiVals = "first")

top_genes_table <- data.frame(
  ensembl_id = rownames(top_genes),
  gene_symbol = gene_symbols,
  baseMean = round(top_genes$baseMean, 1),
  log2FoldChange = round(top_genes$log2FoldChange, 2),
  pvalue = format(top_genes$pvalue, scientific = TRUE, digits = 3),
  padj = format(top_genes$padj, scientific = TRUE, digits = 3)
)

cat("Top 20 differentially expressed genes:\n")
print(top_genes_table)

# Save results table
write.csv(top_genes_table, "results/top_de_genes.csv", row.names = FALSE)
cat("Saved: results/top_de_genes.csv\n")

# =============================================================================
# HEATMAP OF TOP GENES
# =============================================================================

cat("\n=== HEATMAP OF TOP GENES ===\n")

# Get normalized counts for top 50 genes
top50_genes <- rownames(head(res_ordered, 50))
normalized_counts <- counts(dds, normalized = TRUE)
top50_counts <- normalized_counts[top50_genes, ]

# Log transform for heatmap
log_counts <- log2(top50_counts + 1)

# Sample annotation for heatmap
col_annotation <- data.frame(
  Treatment = colData(dds)$dex,
  Cell_Line = colData(dds)$cell
)
rownames(col_annotation) <- colnames(log_counts)

# Create heatmap
png("figures/16_top_genes_heatmap.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(log_counts,
         annotation_col = col_annotation,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes")
dev.off()
cat("Saved: figures/16_top_genes_heatmap.png\n")

# =============================================================================
# PATHWAY ENRICHMENT ANALYSIS
# =============================================================================

cat("\n=== PATHWAY ENRICHMENT ANALYSIS ===\n")

# Prepare gene list for GSEA
# Remove genes with NA values
res_clean <- res[!is.na(res$stat), ]

# Create ranked gene list (rank by test statistic)
gene_list <- res_clean$stat
names(gene_list) <- rownames(res_clean)
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Gene list for GSEA:", length(gene_list), "genes\n")

# Load pathway databases
# We'll use Hallmark gene sets from MSigDB
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}
library(msigdbr)

# Get Hallmark pathways
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$ensembl_gene, hallmark_sets$gs_name)

cat("Hallmark pathways loaded:", length(hallmark_list), "pathways\n")

# Run GSEA
set.seed(123)
fgsea_results <- fgsea(pathways = hallmark_list,
                       stats = gene_list,
                       minSize = 15,
                       maxSize = 500,
                       nperm = 10000)

# Order by significance
fgsea_results <- fgsea_results[order(fgsea_results$padj), ]

cat("Significant pathways (FDR < 0.05):", sum(fgsea_results$padj < 0.05, na.rm = TRUE), "\n")

# Display top pathways
cat("\nTop 10 enriched pathways:\n")
top_pathways <- fgsea_results[1:10, c("pathway", "pval", "padj", "ES", "NES", "size")]
print(top_pathways)

# Save pathway results (convert from data.table first)
fgsea_df <- as.data.frame(fgsea_results)

# Create a clean version for CSV (remove list columns)
fgsea_for_csv <- fgsea_df[, !names(fgsea_df) %in% "leadingEdge"]
fgsea_for_csv$pval <- format(fgsea_for_csv$pval, scientific = TRUE, digits = 3)
fgsea_for_csv$padj <- format(fgsea_for_csv$padj, scientific = TRUE, digits = 3)

write.csv(fgsea_for_csv, "results/pathway_enrichment.csv", row.names = FALSE)
cat("Saved: results/pathway_enrichment.csv\n")

# Save complete results (including leadingEdge) as RDS
saveRDS(fgsea_results, "results/fgsea_complete_results.rds")
# =============================================================================
# PATHWAY VISUALIZATION
# =============================================================================

cat("\n=== PATHWAY VISUALIZATION ===\n")

# Create enrichment score plot for top pathways
top10_pathways <- head(fgsea_results$pathway, 10)

# Enrichment score plot
enrichment_data <- fgsea_results[1:20, ] %>%
  mutate(pathway_short = gsub("HALLMARK_", "", pathway),
         pathway_short = gsub("_", " ", pathway_short),
         direction = ifelse(NES > 0, "Upregulated", "Downregulated"))

p2 <- ggplot(enrichment_data, aes(x = reorder(pathway_short, NES), y = NES, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "Top 20 Enriched Pathways",
       subtitle = "Hallmark Gene Sets (FDR < 0.05)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)",
       fill = "Direction") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

ggsave("figures/17_pathway_enrichment.png", p2, width = 12, height = 10)
cat("Saved: figures/17_pathway_enrichment.png\n")

# =============================================================================
# SAVE FINAL RESULTS
# =============================================================================

cat("\n=== SAVING FINAL RESULTS ===\n")

# Compile all results
final_results <- list(
  deseq_results = res,
  top_genes = top_genes_table,
  pathway_results = fgsea_results,
  summary = list(
    total_genes_tested = nrow(res),
    significant_genes = n_sig,
    upregulated = n_up,
    downregulated = n_down,
    significant_pathways = sum(fgsea_results$padj < 0.05, na.rm = TRUE)
  ),
  dds_object = dds
)

saveRDS(final_results, "results/differential_expression_results.rds")

cat("Differential expression and pathway analysis complete!\n")
cat("Created 5 new figures (13-17)\n")
cat("Results saved to: results/\n")
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Significant genes:", n_sig, "\n")
cat("Upregulated genes:", n_up, "\n")
cat("Downregulated genes:", n_down, "\n")
cat("Significant pathways:", sum(fgsea_results$padj < 0.05, na.rm = TRUE), "\n")

