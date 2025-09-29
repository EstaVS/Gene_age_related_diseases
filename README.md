# Gene Age and Age-Related Diseases

This repository contains code and data for the analysis presented in the study: **"Investigating the Relationship Between Gene Evolutionary Age and Susceptibility to Age-Related Diseases"**.

**Author:** Esta V. Shrewsbury, 
**Supervisor:** Dr Dan V. Nicolau

## Overview

This project aims to explore the evolutionary trajectory of genes involved in age-related diseases. The central hypothesis is that genes associated with age-related diseases (specifically cancer) may have distinct evolutionary ages compared to genes not linked to such conditions. I analysed whether these disease-associated genes originated at specific points in evolutionary history (e.g. metazoan, vertebrate-specific, or are ancient) and what functional implications this may have.

## Key Findings (Summary)

- Gene evolutionary age (phylostrata) contributes modestly but significantly to cancer transcriptomic variance, as shown by ANOVA across 29 cancer types from TCGA and OncoDB.

- Mid-age phylostrata (PS9–PS12) — associated with vertebrate development and immune function — are consistently upregulated across many cancers, suggesting a conserved transcriptional program.

- GSEA revealed lineage-specific enrichment patterns:

    - KIRC: Enriched in PS10–PS12 genes involved in immune signaling and angiogenesis (e.g. *VEGFA*, *CXCL10*).
    
    - LIHC: Showed dual enrichment — ancient PS2 genes (e.g. *UBE2C*, *MAGEA3*) and younger *PS14* genes (e.g. *GPC3*).
    
    - BRCA: Dataset-specific patterns — OncoDB highlighted ancient mitotic genes (PS3–PS4), while TCGA showed repression of developmental genes (PS7).

- PCA of phylostrata enrichment scores revealed clustering of cancer types by evolutionary expression profiles, highlighting shared oncogenic mechanisms.

- Phylostratigraphy provides a novel lens for understanding cancer heterogeneity, revealing how tumours co-opt evolutionarily conserved gene modules to support proliferation, immune evasion, and dedifferentiation.

## Prerequisites
Of course! Here is a clear and concise **Prerequisites** section for your GitHub README, detailing the software, data, and computational resources needed to run the analysis pipeline.

---

## 📋 Prerequisites

Before running the analysis pipeline, ensure you have the following installed and configured:

### 1. Software & Programming Language

*   **R** (Version ≥ 4.1.0)
*   **RStudio** (Recommended, for interactive analysis)

### 2. R Packages

The core analysis relies on the following R packages, which can be installed from CRAN and Bioconductor:

*   **Data Download & Wrangling:** `TCGAbiolinks`, `dplyr`, `readr`, `stringr`, `tibble`
*   **Differential Expression:** `DESeq2`, `limma`
*   **Statistical Modeling:** `aov` (base R)
*   **Enrichment Analysis:** `fgsea`
*   **Visualization:** `ggplot2`, `factoextra`, `patchwork`, `ComplexHeatmap`

### 3. Containerization (Optional but Recommended)

*   **Docker**: For a fully reproducible environment, a Docker snapshot based on **Bioconductor version 3.21** is provided. This ensures all package versions and dependencies are consistent.

### 4. High-Performance Computing (HPC) Environment

*   **SLURM Workload Manager**: The pipeline is designed to be executed on an HPC cluster using SLURM for job scheduling and parallel processing.
*   **Bash Shell**: Required for executing master scripts and job submission.

### 5. Data Sources

The analysis requires transcriptomic data from the following sources. Scripts are provided to facilitate download and organization.

*   **The Cancer Genome Atlas (TCGA):** Raw RNA-seq count data for tumour and normal samples. Downloaded via the `TCGAbiolinks` R package.
*   **OncoDB:** Pre-processed differential expression data. (Note: This dataset includes only the top ~1,000 differentially expressed genes per cancer type).

### 6. Phylostratigraphic Annotation File

A pre-compiled gene-to-phylostrata mapping file is required for evolutionary age annotation. This file is included in the repository in the `data/` directory.

## Data

### Input Data Sources

The analysis relies on the following key datasets. Due to licensing, the raw data files are not stored in this repository but can be downloaded from the original sources.

1.  **Phylostrata Dataframe:** Privately produced, available in 'Phylostrata/gene_phylostrata.txt'
2.  **OncoDB:** Differentially Expressed Genes from Various Cancers available at *https://oncodb.org/data_downloads.html*
3.  **TCGA:** Database query and direct download to HPC environment via *TCGAbiolinks* R package. 


### Processed Data

The processed TCGA datasets used for the main analysis are located in `data/TCGA_DE_results`.

## Citation

If you use this code or data in your research, please cite our publication: (In Preparation)

> E. V. Shrewsbury *et al.* (2025-2026). Investigating the Relationship Between Gene Evolutionary Age and Susceptibility to Age-Related Diseases. *[Unknown]*. [DOI/Link]


## Acknowledgments

- We thank the developers of the public databases used in this study (OncoDB, The Cancer Genome Atlas).
- This work was supported by King's College London for the completion of the MSc in Applied Bioinformatics.
