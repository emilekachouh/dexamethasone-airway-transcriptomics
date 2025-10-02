# 03_high_dim_analysis.R - High-dimensional analysis of airway data
# PCA, clustering, and dimensionality reduction
install.packages("mclust")
# Load required libraries
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(umap)
library(cluster)
library(dplyr)
library(mclust)
# Create output directories
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

# Load filtered data
cat("Loading filtered airway dataset...\n")
airway_filtered <- readRDS("data/processed/airway_filtered.rds")

cat("Working with", nrow(airway_filtered), "genes and", ncol(airway_filtered), "samples\n")

# Extract sample information
sample_info <- colData(airway_filtered)
expr_matrix <- assay(airway_filtered)

# =============================================================================
# DATA TRANSFORMATION FOR HIGH-DIMENSIONAL ANALYSIS
# =============================================================================

cat("\n=== DATA TRANSFORMATION ===\n")

# Variance Stabilizing Transformation (VST)
# This is crucial for PCA with count data
cat("Applying variance stabilizing transformation...\n")
vst_data <- vst(expr_matrix)

cat("Original data range:", round(min(expr_matrix), 2), "to", round(max(expr_matrix), 2), "\n")
cat("VST data range:", round(min(vst_data), 2), "to", round(max(vst_data), 2), "\n")

# =============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA)
# =============================================================================

cat("\n=== PRINCIPAL COMPONENT ANALYSIS ===\n")

# Perform PCA on most variable genes
ntop <- 1000  # Use top 1000 most variable genes
rv <- rowVars(vst_data)
select_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca_data <- vst_data[select_genes, ]

cat("PCA on top", length(select_genes), "most variable genes\n")

# Perform PCA
pca_result <- prcomp(t(pca_data), center = TRUE, scale. = FALSE)

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cat("Variance explained by PC1:", round(var_explained[1] * 100, 1), "%\n")
cat("Variance explained by PC2:", round(var_explained[2] * 100, 1), "%\n")
cat("Cumulative variance (PC1-PC2):", round(sum(var_explained[1:2]) * 100, 1), "%\n")

# Create PCA data frame for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = colnames(airway_filtered),
  Treatment = sample_info$dex,
  Cell_Line = sample_info$cell
)

# PCA plot colored by treatment
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Cell_Line)) +
  geom_point(size = 4) +
  labs(title = "PCA: Samples Colored by Treatment",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "% variance)")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/07_pca_treatment.png", p1, width = 10, height = 8)
cat("Saved: figures/07_pca_treatment.png\n")

# PCA plot colored by cell line
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell_Line, shape = Treatment)) +
  geom_point(size = 4) +
  labs(title = "PCA: Samples Colored by Cell Line",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "% variance)")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/08_pca_cellline.png", p2, width = 10, height = 8)
cat("Saved: figures/08_pca_cellline.png\n")

# Variance explained plot (scree plot)
var_df <- data.frame(
  PC = 1:8,
  Variance_Explained = var_explained * 100,
  Cumulative = cumsum(var_explained) * 100
)

p3 <- ggplot(var_df, aes(x = PC, y = Variance_Explained)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = Cumulative), color = "red", size = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 2) +
  labs(title = "PCA Variance Explained",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  scale_x_continuous(breaks = 1:8) +
  theme_minimal()

ggsave("figures/09_pca_variance.png", p3, width = 8, height = 6)
cat("Saved: figures/09_pca_variance.png\n")

# =============================================================================
# GENE LOADINGS ANALYSIS
# =============================================================================

cat("\n=== GENE LOADINGS ANALYSIS ===\n")

# Extract loadings for PC1 and PC2
loadings_df <- data.frame(
  Gene = rownames(pca_result$rotation),
  PC1_loading = pca_result$rotation[, 1],
  PC2_loading = pca_result$rotation[, 2]
)

# Find genes with strongest loadings
pc1_top <- loadings_df %>% 
  arrange(desc(abs(PC1_loading))) %>% 
  head(10)

pc2_top <- loadings_df %>% 
  arrange(desc(abs(PC2_loading))) %>% 
  head(10)

cat("Top 10 genes driving PC1:\n")
print(pc1_top$Gene)

cat("Top 10 genes driving PC2:\n")
print(pc2_top$Gene)

# =============================================================================
# UMAP NON-LINEAR DIMENSIONALITY REDUCTION
# =============================================================================

cat("\n=== UMAP ANALYSIS ===\n")

# UMAP for comparison with PCA
set.seed(123)
umap_result <- umap(t(pca_data), n_neighbors = 4, min_dist = 0.3)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Sample = colnames(airway_filtered),
  Treatment = sample_info$dex,
  Cell_Line = sample_info$cell
)

p4 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Treatment, shape = Cell_Line)) +
  geom_point(size = 4) +
  labs(title = "UMAP: Non-linear Dimensionality Reduction",
       x = "UMAP1",
       y = "UMAP2") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/10_umap.png", p4, width = 10, height = 8)
cat("Saved: figures/10_umap.png\n")

# =============================================================================
# HIERARCHICAL CLUSTERING
# =============================================================================

cat("\n=== HIERARCHICAL CLUSTERING ===\n")

# Distance matrix and clustering
sample_distances <- dist(t(pca_data))
hc <- hclust(sample_distances, method = "complete")

# Create dendrogram with annotations
annotation_df <- data.frame(
  Treatment = sample_info$dex,
  Cell_Line = sample_info$cell
)
rownames(annotation_df) <- colnames(airway_filtered)

# Plot dendrogram
png("figures/11_dendrogram.png", width = 10, height = 6, units = "in", res = 300)
plot(hc, main = "Hierarchical Clustering of Samples", xlab = "Samples", sub = "")
dev.off()
cat("Saved: figures/11_dendrogram.png\n")

# =============================================================================
# CLUSTERING VALIDATION
# =============================================================================

cat("\n=== CLUSTERING VALIDATION ===\n")

# K-means clustering
set.seed(123)
k2 <- kmeans(t(pca_data), centers = 2, nstart = 25)
k4 <- kmeans(t(pca_data), centers = 4, nstart = 25)

# Silhouette analysis
sil_k2 <- silhouette(k2$cluster, sample_distances)
sil_k4 <- silhouette(k4$cluster, sample_distances)

cat("K=2 clustering - Average silhouette width:", round(mean(sil_k2[, 3]), 3), "\n")
cat("K=4 clustering - Average silhouette width:", round(mean(sil_k4[, 3]), 3), "\n")

# Compare clustering with known groups
treatment_clusters <- ifelse(sample_info$dex == "trt", 1, 2)
cluster_agreement_k2 <- adjustedRandIndex(k2$cluster, treatment_clusters)
cat("K=2 vs Treatment agreement (ARI):", round(cluster_agreement_k2, 3), "\n")

# Add clustering results to PCA plot
pca_df$K2_cluster <- as.factor(k2$cluster)
pca_df$K4_cluster <- as.factor(k4$cluster)

p5 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = K2_cluster, shape = Treatment)) +
  geom_point(size = 4) +
  labs(title = "PCA with K-means Clustering (K=2)",
       x = paste0("PC1 (", round(var_explained[1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(var_explained[2] * 100, 1), "% variance)"),
       color = "Cluster") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/12_pca_clusters.png", p5, width = 10, height = 8)
cat("Saved: figures/12_pca_clusters.png\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save analysis results
high_dim_results <- list(
  pca_result = pca_result,
  variance_explained = var_explained,
  pca_coordinates = pca_df,
  umap_result = umap_result,
  umap_coordinates = umap_df,
  gene_loadings = loadings_df,
  clustering = list(
    k2 = k2,
    k4 = k4,
    silhouette_k2 = mean(sil_k2[, 3]),
    silhouette_k4 = mean(sil_k4[, 3]),
    treatment_agreement = cluster_agreement_k2
  ),
  vst_data = vst_data
)

saveRDS(high_dim_results, "results/high_dim_analysis.rds")

cat("High-dimensional analysis complete!\n")
cat("Created 6 new figures (07-12)\n")
cat("Results saved to: results/high_dim_analysis.rds\n")
