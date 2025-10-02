# =============================================================================
# SETUP SCRIPT - Install Required Packages
# Run this once before starting the analysis
# =============================================================================

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "SummarizedExperiment", 
  "airway",
  "fgsea",
  "org.Hs.eg.db"
)

BiocManager::install(bioc_packages, update = FALSE)

# CRAN packages
cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "patchwork",
  "cowplot",
  "gridExtra",
  "scales",
  "stringr",
  "umap",
  "ggdendro"
)

install.packages(cran_packages)

# Verify installation
cat("\n=== Checking package installation ===\n")
all_packages <- c(bioc_packages, cran_packages)

for(pkg in all_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "installed successfully\n")
  } else {
    cat("✗", pkg, "FAILED to install\n")
  }
}

cat("\nSetup complete! You're ready to run the analysis.\n")
