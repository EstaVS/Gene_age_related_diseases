### Master Script: TCGA-Wide Phylostrata Differential Expression Analysis ###

setwd("/workdir")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(umap)

# Load phylostrata table
phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID)

# Prepare output folders
dir.create("DE_results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
dir.create("stats", showWarnings = FALSE)

# Get TCGA project IDs
projects <- TCGAbiolinks::getGDCprojects()$project_id
tcga_projects <- projects[grepl("^TCGA-", projects)]

# Store all DE results
de_results_list <- list()
skipped_projects <- c()

### Loop over TCGA projects ###
for (project_id in tcga_projects) {
  tryCatch({
    cat("\nProcessing:", project_id, "\n")
    
    query <- GDCquery(
      project = project_id,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    GDCdownload(query)
    data <- GDCprepare(query)
    
    sample_types <- table(data$sample_type)
    tumor_count <- sum(data$sample_type == "Primary Tumor")
    normal_count <- sum(data$sample_type == "Solid Tissue Normal")
    
    if (tumor_count < 3 | normal_count < 3) {
      cat("Skipping", project_id, "due to insufficient tumor or normal samples\n")
      skipped_projects <- c(skipped_projects, project_id)
      next
    }
    
    combined_data <- data[, data$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
    
    gene_mapping <- as.data.frame(rowData(combined_data)) %>%
      mutate(gene_id = rownames(.)) %>%
      select(gene_id, gene_name)
    
    combined_data <- combined_data[rownames(combined_data) %in% gene_mapping$gene_id, ]
    gene_mapping <- gene_mapping %>% filter(gene_id %in% rownames(combined_data))
    
    dds <- DESeqDataSetFromMatrix(countData = assay(combined_data),
                                  colData = colData(combined_data),
                                  design = ~ sample_type)
    
    # Filter out low-count genes (less than 10 counts in fewer than 3 samples)
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    
    dds <- DESeq(dds)
    
    # Normalised counts
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Extract tumor and normal sample names
    tumor_samples <- colnames(dds)[dds$sample_type == "Primary Tumor"]
    normal_samples <- colnames(dds)[dds$sample_type == "Solid Tissue Normal"]
    
    # Compute gene-wise medians
    tumor_medians <- apply(norm_counts[, tumor_samples], 1, median, na.rm = TRUE)
    normal_medians <- apply(norm_counts[, normal_samples], 1, median, na.rm = TRUE)
    
    # Add medians to DE results
    res_df$cancer_median <- tumor_medians[rownames(res_df)]
    res_df$normal_median <- normal_medians[rownames(res_df)]
    
    # Continue DE
    res_df <- as.data.frame(results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))) %>%
      mutate(gene_id = rownames(.)) %>%
      left_join(gene_mapping, by = "gene_id") %>%
      mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
      select(gene_name, everything()) %>%
      left_join(phylostrata, by = "gene_name") %>%
      filter(!is.na(Phylostrata))
    
    # Save DE
    write_csv(res_df, paste0("DE_results/", project_id, "_DE.csv"))
    
    # Store for combined analysis
    res_df$project_id <- project_id
    de_results_list[[project_id]] <- res_df
    
  }, error = function(e) {
    cat("Error in", project_id, ":", e$message, "\n")
    skipped_projects <- c(skipped_projects, project_id)
  })
}

# Save skipped projects list
write_lines(skipped_projects, "DE_results/skipped_TCGAprojects_due_to_sample_size.txt")

### Combine All Results ###
all_de <- bind_rows(de_results_list)
write_csv(all_de, "DE_results/all_TCGAprojects_combined_DE.csv")