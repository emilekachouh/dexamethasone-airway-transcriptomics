# 01_data_prep.R - Load and prepare airway dataset
# Dexamethasone treatment in human airway epithelial cells

# Load required libraries
library(airway)
library(SummarizedExperiment)
library(dplyr)

# Create data directory if it doesn't exist
if (!dir.exists("data/raw")) {
  dir.create("data/raw", recursive = TRUE)
}

# Load the airway dataset
cat("Loading airway dataset...\n")
data(airway)

# Save the raw data for consistency
saveRDS(airway, file = "data/raw/airway_dataset.rds")

# Basic exploration
cat("\n=== DATASET OVERVIEW ===\n")
cat("Dataset: Human airway epithelial cells treated with dexamethasone\n")
cat("Source: Himes et al. (2014) PLoS One\n")
print(airway)

cat("\n=== EXPERIMENTAL DESIGN ===\n")
sample_info <- colData(airway)
print(sample_info)

cat("\n=== DESIGN SUMMARY ===\n")
cat("Total samples:", ncol(airway), "\n")
cat("Total genes:", nrow(airway), "\n")

cat("Treatment groups:\n")
treatment_table <- table(sample_info$dex)
print(treatment_table)

cat("Cell lines:\n")
cell_table <- table(sample_info$cell)
print(cell_table)

cat("Design: Cell line × Treatment interaction\n")
design_table <- table(sample_info$cell, sample_info$dex)
print(design_table)

cat("\n=== EXPRESSION DATA DETAILS ===\n")
expr_matrix <- assay(airway)
cat("Data type:", class(expr_matrix), "\n")
cat("Data range:", min(expr_matrix), "to", max(expr_matrix), "\n")
cat("Total counts:", sum(expr_matrix), "\n")

# Check for missing data
missing_counts <- sum(is.na(expr_matrix))
cat("Missing values:", missing_counts, "\n")

# Basic count statistics
cat("Mean counts per gene:", round(mean(rowSums(expr_matrix)), 1), "\n")
cat("Mean counts per sample:", round(mean(colSums(expr_matrix)), 1), "\n")

cat("\n=== RESEARCH QUESTIONS ===\n")
cat("1. Which genes respond to dexamethasone treatment?\n")
cat("2. Are there cell-line-specific responses?\n")
cat("3. What biological pathways are affected?\n")
cat("4. Can we identify glucocorticoid receptor targets?\n")

cat("\nData preparation complete!\n")
cat("Raw data saved to: data/raw/airway_dataset.rds\n")
cat("Ready for quality control and analysis.\n")
