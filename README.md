# Dexamethasone Transcriptomic Analysis

**Systematic Transcriptomic Analysis of Dexamethasone-Treated Human Airway Epithelial Cells**

## Overview
This project contains a complete RNA-seq analysis pipeline investigating dexamethasone treatment effects in human airway epithelial cells. The analysis identifies genome-wide transcriptional responses to glucocorticoid treatment and characterizes affected biological pathways.

## Key Findings
- **4,099 differentially expressed genes** (24.7% of transcriptome; FDR < 0.05)
- **18 significantly enriched biological pathways** including anti-inflammatory and metabolic processes
- **Novel identification of SPARCL1** as major glucocorticoid-responsive gene (4.57-fold upregulation)
- Treatment status explains **42.8% of transcriptional variance** (PC1)

## Project Structure
```
breast_cancer_multiomics/
├── data/
│   ├── raw/              # Original datasets (not uploaded - see DATA_README.md)
│   └── processed/        # Processed data (not uploaded - see DATA_README.md)
├── scripts/              # R analysis pipeline (numbered workflow)
├── results/              # Analysis outputs
├── figures/              # Publication-quality figures
├── docs/                 # Manuscript and reports
└── setup.R               # Package installation script
```

## Analysis Pipeline

The analysis follows a systematic 5-step workflow:

1. **Data Preparation** (`01_data_prep.R`) - Load and preprocess RNA-seq data
2. **Quality Control** (`02_eda_qc.R`) - QC metrics, filtering, validation
3. **Dimensionality Reduction** (`03_high_dim_analysis.R`) - PCA, UMAP, clustering
4. **Differential Expression** (`04_diff_expr_pathways.R`) - DESeq2 analysis and GSEA
5. **Publication Figures** (`05_publication_figures.R`) - Generate manuscript figures

## Requirements

- R >= 4.0
- Bioconductor packages: DESeq2, SummarizedExperiment, fgsea
- CRAN packages: ggplot2, dplyr, pheatmap, patchwork

See `setup.R` for complete package list and installation instructions.

## Usage

```r
# Install all required packages
source("setup.R")

# Run complete analysis pipeline
source("scripts/01_data_prep.R")
source("scripts/02_eda_qc.R")
source("scripts/03_high_dim_analysis.R")
source("scripts/04_diff_expr_pathways.R")
source("scripts/05_publication_figures.R")
```

## Results Summary

- **Upregulated pathways**: Adipogenesis, TNF-α signaling, oxidative phosphorylation
- **Downregulated pathways**: P53 signaling, E2F targets, interferon response
- **Top upregulated genes**: SPARCL1, KLF15, DUSP1, CACNB2
- **Top downregulated genes**: VCAM1, KCTD12, ADAM12

## Data Availability

Raw data files are not included due to GitHub size limits. See `DATA_README.md` for instructions on obtaining the original datasets.

## Citation

If you use this analysis pipeline, please cite:

```
Kachouh, E. (2025). Systematic Transcriptomic Analysis Reveals Coordinated 
Anti-Inflammatory and Metabolic Reprogramming in Dexamethasone-Treated 
Human Airway Epithelial Cells. [Manuscript in preparation]
```

## Contact

For questions or data requests, please open an issue or contact [emilekachouh1@gmail.com].


