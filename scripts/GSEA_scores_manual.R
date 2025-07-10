setwd("~/Desktop/gene_age_project/")

library(readr)
library(dplyr)
library(stringr)

# List your OncoDB GSEA files
gsea_files <- list.files("stats/TCGA/GSEA", pattern = "_gsea_results.csv$", full.names = TRUE)

# Loop over them
for (file in gsea_files) {
  cancer <- str_extract(basename(file), "^[^_]+")
  gsea <- read_csv(file, show_col_types = FALSE)
  
  strongest <- gsea %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    select(pathway, NES, padj)
  
  cat("\n======", cancer, "======\n")
  if (nrow(strongest) == 0) {
    cat("No significant pathways (padj < 0.05)\n")
  } else {
    print(strongest)
  }
}

