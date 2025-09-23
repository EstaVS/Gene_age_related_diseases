# Gene Age and Age-Related Diseases

This repository contains code and data for the analysis presented in the study: **"Investigating the Relationship Between Gene Evolutionary Age and Susceptibility to Age-Related Diseases"**.

**Authors:** Esta V. Shrewsbury
**Supervisor:** Dr Dan V. Nicolau
**In affiliation with...* King's College London, for the completion of the MSc in Applied Bioinformatics

## Overview

This project aims to explore the evolutionary trajectory of genes involved in age-related diseases. The central hypothesis is that genes associated with age-related diseases (specifically cancer) may have distinct evolutionary ages compared to genes not linked to such conditions. I analysed whether these disease-associated genes originated at specific points in evolutionary history (e.g. metazoan, vertebrate-specific, or are ancient) and what functional implications this may have.

## Repository Structure

```
Gene_age_related_diseases/
├── data/
│   ├── raw/                   # Original, immutable data
│   │   ├── gene_ages/         # Gene age datasets (e.g., from GenTree, EggNOG)
│   │   ├── disease_genes/     # Gene-disease associations (e.g., from DisGeNET, OMIM)
│   │   └── functional_annotations/ # GO, KEGG, Pathway annotations
│   └── processed/             # Cleaned, analysis-ready datasets
│       ├── merged_dataset.csv
│       └── analysis_subset.tsv
├── notebooks/
│   ├── 01_data_preprocessing.ipynb
│   ├── 02_statistical_analysis.ipynb
│   └── 03_visualization.ipynb
├── scripts/
│   ├── parse_gene_ages.py
│   ├── enrichment_analysis.R
│   └── functions.py           # Helper functions
├── results/
│   ├── figures/               # Generated plots (PDF/PNG)
│   │   ├── figure1_age_distribution.pdf
│   │   └── figure2_enrichment_plot.pdf
│   └── tables/                # Statistical results
│       └── table1_summary_stats.csv
├── docs/                      # Supplementary documentation
├── environment.yml            # Conda environment for reproducibility
└── README.md                  # This file
```

*(Note: The structure above is a suggestion. Please adjust the folder and file names to match your actual repository.)*

## Key Findings (Summary)

*A brief summary of your main results. For example:*
- Genes associated with age-related diseases are significantly enriched for [e.g., ancient evolutionary origin / specific evolutionary age categories].
- The functional profile of disease-associated genes differs based on their evolutionary age.
- [Other key finding 1]
- [Other key finding 2]

## Prerequisites

Before running the code, ensure you have the following installed:

- **Python 3.8+**
- **R 4.0+** (if R scripts are included)
- Package managers: `pip` and/or `conda`

## Installation & Setup

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/EstaVS/Gene_age_related_diseases.git
    cd Gene_age_related_diseases
    ```

2.  **Set up the computational environment.**

    **Using Conda (Recommended):**
    The `environment.yml` file contains all necessary dependencies.
    ```bash
    conda env create -f environment.yml
    conda activate gene_age_diseases
    ```

    *Example `environment.yml` content:*
    ```yaml
    name: gene_age_diseases
    channels:
      - conda-forge
      - bioconda
      - defaults
    dependencies:
      - python=3.9
      - pandas
      - numpy
      - scipy
      - matplotlib
      - seaborn
      - jupyter
      - r-base=4.1
      - r-ggplot2
      - r-dplyr
    ```

    **Using pip:**
    ```bash
    pip install -r requirements.txt
    ```

## Data

### Input Data Sources

The analysis relies on the following key datasets. Due to licensing, the raw data files are not stored in this repository but can be downloaded from the original sources.

1.  **Gene Evolutionary Ages:** Obtained from [GenTree](https://gentree.ioz.ac.cn/), [EggNOG](http://eggnog5.embl.de/), or similar databases.
2.  **Age-Related Disease Genes:** Curated from [DisGeNET](https://www.disgenet.org/), [OMIM](https://www.omim.org/), or literature review.
3.  **Functional Annotations:** Gene Ontology (GO) and KEGG pathways from [Ensembl BioMart](https://www.ensembl.org/biomart/martview) or similar.

*Instructions for downloading the raw data should be placed in `data/raw/README.md`.*

### Processed Data

The final, merged dataset used for the main analysis is located in `data/processed/merged_dataset.csv`.

## Usage

The analysis is divided into logical steps, typically run in order.

### 1. Data Preprocessing

Run the Jupyter notebook to clean and merge the raw data files.
```bash
jupyter notebook notebooks/01_data_preprocessing.ipynb
```
*This script will:*
- Load gene age and disease association data.
- Map gene identifiers to a common standard (e.g., Ensembl ID).
- Merge datasets into a single analysis table.

### 2. Statistical Analysis

Run the main analysis notebook to perform statistical tests.
```bash
jupyter notebook notebooks/02_statistical_analysis.ipynb
```
*This script will:*
- Test for enrichment of age-related disease genes in different evolutionary age categories (e.g., using Chi-squared or Fisher's exact tests).
- Perform functional enrichment analysis (e.g., using GOseq or clusterProfiler).

### 3. Generate Figures and Results

Run the visualization notebook to recreate all publication-ready figures.
```bash
jupyter notebook notebooks/03_visualization.ipynb
```

## Citation

If you use this code or data in your research, please cite our publication:

> [Author Names]. (Year). Investigating the Relationship Between Gene Evolutionary Age and Susceptibility to Age-Related Diseases. *[Journal Name]*. [DOI/Link]

*(Replace this with your actual citation once available)*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. (*Choose an appropriate license: MIT, GPL, etc.*)

## Acknowledgments

- We thank the developers of the public databases used in this study (DisGeNET, GenTree, etc.).
- This work was supported by [Your Funding Source].
