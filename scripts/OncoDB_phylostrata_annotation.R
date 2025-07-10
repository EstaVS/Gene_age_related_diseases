### Adding phylostrata column to OncoDB DE raw data files ###

setwd("~/Desktop/gene_age_project/")

library(readr)
library(dplyr)
library(stringr)
library(purrr)

# Load and clean phylostrata
phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID) %>%
  mutate(
    gene_name = toupper(trimws(gene_name))
  ) %>%
  distinct(gene_name, .keep_all = TRUE)

# Create output folder
dir.create("OncoDB_data_annotated", showWarnings = FALSE)

# List DE files
files <- list.files("OncoDB_data", pattern = "^expression_diff_.*\\.txt$", full.names = TRUE)

for (file in files) {
  # Extract project ID
  project_id <- str_extract(basename(file), "(?<=expression_diff_)[^\\.]+")
  cancer_type_name <- project_info %>%
    filter(project_id == !!project_id) %>%
    pull(cancer_type)
  
  cat("Processing:", project_id, "-", cancer_type_name, "\n")
  
  # Read and clean expression data
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    setNames(c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")) %>%
    mutate(
      gene_name = toupper(trimws(gene_name)),
      padj = as.numeric(gsub("q<=|Q<=|Q=|q=", "", padj))
    )
  
  # Join with phylostrata
  annotated <- data %>%
    inner_join(phylostrata, by = "gene_name") %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange),
      !is.na(Phylostrata)
    )
  
  cat("  → Matched genes:", nrow(annotated), "\n")
  
  # Save annotated data
  write.csv(annotated,
              file = paste0("OncoDB_data_annotated/expression_diff_", project_id, ".csv"),
              row.names = FALSE, quote = FALSE)
}
