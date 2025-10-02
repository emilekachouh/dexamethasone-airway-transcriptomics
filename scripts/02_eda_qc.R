# 02_eda_qc.R - Quality Control for Airway Dataset
# Dexamethasone treatment in airway epithelial cells

# Load required libraries
library(airway)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)

# Create output directories
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("data/processed")) dir.create("data/processed")

# Load the airway dataset
cat("Loading airway dataset...\n")
data(airway)

# Basic data structure
cat("=== BASIC DATA STRUCTURE ===\n")
cat("Genes:", nrow(airway), "\n")
cat("Samples:", ncol(airway), "\n")
print(airway)

# Extract components
expr_matrix <- assay(airway)  # Raw count matrix
sample_info <- colData(airway)

cat("Expression data type:", class(expr_matrix), "\n")
cat("Expression range:", min(expr_matrix), "to", max(expr_matrix), "\n")

# Show sample information
cat("\n=== SAMPLE INFORMATION ===\n")
print(sample_info)
cat("Treatment groups:\n")
print(table(sample_info$dex))
cat("Cell lines:\n")
print(table(sample_info$cell))

# =============================================================================
# QUALITY CONTROL CHECKS
# =============================================================================

cat("\n=== SAMPLE QUALITY CONTROL ===\n")

# 1. Library size distribution (total counts per sample)
library_sizes <- colSums(expr_matrix)
cat("Library size range:", min(library_sizes), "to", max(library_sizes), "\n")

# Create library size plot
library_df <- data.frame(
  sample = colnames(expr_matrix),
  library_size = library_sizes,
  treatment = sample_info$dex,
  cell_line = sample_info$cell
)

p1 <- ggplot(library_df, aes(x = treatment, y = library_size, fill = treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Library Sizes by Treatment",
       x = "Dexamethasone Treatment",
       y = "Total Counts per Sample") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("figures/01_library_sizes.png", p1, width = 8, height = 6)
cat("Saved: figures/01_library_sizes.png\n")

# 2. Number of detected genes per sample
detected_genes <- colSums(expr_matrix > 0)
cat("Detected genes per sample - range:", min(detected_genes), "to", max(detected_genes), "\n")

detected_df <- data.frame(
  sample = colnames(expr_matrix),
  detected_genes = detected_genes,
  treatment = sample_info$dex
)

p2 <- ggplot(detected_df, aes(x = treatment, y = detected_genes, fill = treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2)) +
  labs(title = "Number of Detected Genes per Sample",
       x = "Dexamethasone Treatment",
       y = "Genes with Counts > 0") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("figures/02_detected_genes.png", p2, width = 8, height = 6)
cat("Saved: figures/02_detected_genes.png\n")

# 3. Sample correlation analysis
cat("\n=== SAMPLE CORRELATION ANALYSIS ===\n")

# Use variance stabilizing transformation for correlation
vst_data <- vst(expr_matrix)
sample_cors <- cor(vst_data)

# Create correlation heatmap with sample annotations
annotation_col <- data.frame(
  Treatment = sample_info$dex,
  Cell_Line = sample_info$cell
)
rownames(annotation_col) <- colnames(expr_matrix)

png("figures/03_sample_correlation_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(sample_cors,
         annotation_col = annotation_col,
         main = "Sample-Sample Correlations",
         show_rownames = TRUE,
         show_colnames = TRUE)
dev.off()
cat("Saved: figures/03_sample_correlation_heatmap.png\n")

# =============================================================================
# GENE-LEVEL QUALITY CONTROL
# =============================================================================

cat("\n=== GENE-LEVEL QUALITY CONTROL ===\n")

# 1. Gene expression distribution
gene_means <- rowMeans(expr_matrix)
gene_vars <- apply(expr_matrix, 1, var)

cat("Gene expression means - range:", round(min(gene_means), 2), "to", round(max(gene_means), 2), "\n")

# Create mean expression plot
mean_df <- data.frame(gene_mean = gene_means)

p3 <- ggplot(mean_df, aes(x = gene_mean)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  scale_x_log10() +
  labs(title = "Distribution of Gene Expression Means",
       x = "Mean Count (log10 scale)",
       y = "Number of Genes") +
  theme_minimal()

ggsave("figures/04_gene_means.png", p3, width = 8, height = 6)
cat("Saved: figures/04_gene_means.png\n")

# 2. Mean-variance relationship (important for count data)
mv_df <- data.frame(
  mean = gene_means,
  variance = gene_vars
) %>%
  filter(mean > 0, variance > 0)  # Remove zeros for log scale

p4 <- ggplot(mv_df, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "loess", color = "red") +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  labs(title = "Mean-Variance Relationship",
       subtitle = "Blue line: Poisson expectation (variance = mean)",
       x = "Mean Count (log10)",
       y = "Variance (log10)") +
  theme_minimal()

ggsave("figures/05_mean_variance.png", p4, width = 8, height = 6)
cat("Saved: figures/05_mean_variance.png\n")

# 3. Filter low-expression genes
min_samples <- 3  # At least 3 samples (we only have 8 total)
min_count <- 10   # At least 10 counts

keep_genes <- rowSums(expr_matrix >= min_count) >= min_samples
cat("Genes passing filter (count >=", min_count, "in >=", min_samples, "samples):", sum(keep_genes), "\n")
cat("Genes removed:", sum(!keep_genes), "\n")
cat("Proportion kept:", round(mean(keep_genes), 3), "\n")

# =============================================================================
# STATISTICAL POWER ANALYSIS
# =============================================================================

cat("\n=== POWER ANALYSIS ===\n")

# Simple simulation: what fold-changes can we detect with n=4 per group?
set.seed(123)
n_per_group <- 4
base_mean <- 100
fold_changes <- c(1.2, 1.5, 2.0, 3.0, 4.0)
n_simulations <- 1000

power_results <- sapply(fold_changes, function(fc) {
  p_values <- replicate(n_simulations, {
    # Simulate counts under alternative hypothesis
    group1 <- rnbinom(n_per_group, mu = base_mean, size = 10)
    group2 <- rnbinom(n_per_group, mu = base_mean * fc, size = 10)
    
    # t-test on log-transformed data (quick approximation)
    if(all(group1 > 0) & all(group2 > 0)) {
      t.test(log(group1), log(group2))$p.value
    } else {
      1  # Conservative: assign p=1 if zeros present
    }
  })
  
  # Power = proportion of p-values < 0.05
  mean(p_values < 0.05)
})

power_df <- data.frame(
  fold_change = fold_changes,
  power = power_results
)

p5 <- ggplot(power_df, aes(x = fold_change, y = power)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(title = "Statistical Power Analysis",
       subtitle = "Power to detect fold-changes with n=4 per group",
       x = "Fold Change",
       y = "Power (1 - β)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()

ggsave("figures/06_power_analysis.png", p5, width = 8, height = 6)
cat("Saved: figures/06_power_analysis.png\n")

# =============================================================================
# SAVE PROCESSED DATA
# =============================================================================

cat("\n=== SAVING PROCESSED DATA ===\n")

# Create filtered dataset
airway_filtered <- airway[keep_genes, ]
cat("Filtered dataset:", nrow(airway_filtered), "genes x", ncol(airway_filtered), "samples\n")

# Save quality control results
qc_results <- list(
  library_sizes = library_sizes,
  detected_genes = detected_genes,
  gene_means = gene_means,
  gene_vars = gene_vars,
  sample_correlations = sample_cors,
  power_analysis = power_df,
  filter_summary = list(
    original_genes = nrow(airway),
    filtered_genes = nrow(airway_filtered),
    removed_genes = sum(!keep_genes),
    filter_criteria = paste("count >=", min_count, "in >=", min_samples, "samples")
  )
)

saveRDS(qc_results, "data/processed/qc_results.rds")
saveRDS(airway_filtered, "data/processed/airway_filtered.rds")

cat("\nQuality control complete!\n")
cat("Dataset: Airway epithelial cells treated with dexamethasone\n")
cat("Design: 4 cell lines x 2 treatments (control vs dexamethasone)\n")
cat("Results saved to: data/processed/\n")
cat("Figures saved to: figures/\n")
