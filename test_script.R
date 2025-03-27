setwd("~/Desktop/gene_age_project/")

### Download data ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

breast_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
) # Construct a query to fetch data based on specific filters

GDCdownload(breast_query) # Download the data

breast_data <- GDCprepare(breast_query) # Load the data in R

### Process data ###

# Install SummarizedExperiment
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

# Load the package
library(SummarizedExperiment)

# View sample metadata
metadata <- colData(breast_data)
print(head(metadata))

# Check the unique values in the 'sample_type' column
table(metadata$sample_type)

# Subset tumor samples
tumor_samples <- breast_data[, breast_data$sample_type == "Primary Tumor"]

# Subset normal samples
normal_samples <- breast_data[, breast_data$sample_type == "Solid Tissue Normal"]

# Combine tumor and normal samples into a single object
combined_data <- cbind(tumor_samples, normal_samples)

### Data Exploration & Cleaning ###

library(dplyr)

# Extract gene metadata
gene_info <- as.data.frame(rowData(breast_data))

# View the first few rows
head(gene_info)

# Check if one gene_id maps to multiple gene_names
dup_ids <- gene_info %>%
  group_by(gene_id) %>%
  summarize(unique_gene_names = n_distinct(gene_name)) %>%
  filter(unique_gene_names > 1)

print(dup_ids)  # If empty, IDs and names are perfectly correlated

# Check if one gene_name maps to multiple gene_ids
dup_names <- gene_info %>%
  group_by(gene_name) %>%
  summarize(unique_gene_ids = n_distinct(gene_id)) %>%
  filter(unique_gene_ids > 1)

print(dup_names)  # If empty, no duplicate mappings exist

colnames(rowData(combined_data))

# Convert rowData to a data frame
gene_mapping <- as.data.frame(rowData(combined_data))

# Add row names (Ensembl gene IDs) as a separate column
gene_mapping$gene_id <- rownames(gene_mapping)

# Select only relevant columns
gene_mapping <- gene_mapping %>%
  select(gene_id, gene_name)

### Differential Expression ###

# Load DESeq2
library(DESeq2)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = assay(combined_data),  # Expression data
  colData = colData(combined_data),  # Sample metadata
  design = ~ sample_type             # Design formula (tumor vs normal)
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
results <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))

# View the results
head(results)

# Filter for significant genes
significant_genes <- results[!is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) > 1, ]

# View significant genes
print(significant_genes)

# Adding gene names to significant genes dataframe
res_df <- as.data.frame(significant_genes) %>%
  mutate(gene_id = rownames(.))  # Ensure gene IDs are a column

# Merge to add gene names
res_annotated <- merge(res_df, gene_mapping, by = "gene_id", all.x = TRUE)

# Move gene_name to the first column
res_annotated <- res_annotated %>%
  select(gene_name, everything()) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))  # Fill missing names with gene IDs

head(res_annotated)

# All genes
all_genes <- as.data.frame(results) %>% 
  mutate(gene_id = rownames(.))
all_genes <- merge(all_genes, gene_mapping, by = "gene_id", all.x = TRUE)
all_genes <- all_genes %>%
  select(gene_name, everything()) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

### Saving the all_genes file for other scripts ###
write_csv(all_genes, "differential_expression_data.csv")

### Visualising the data ### Pre-Phylostrata
library(ggplot2)

# Sort by absolute log2FoldChange and select top 20
top_20_genes <- res_annotated %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)

# Create the bar plot
ggplot(top_20_genes, aes(x = gene_name, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("blue", "red")) + # Colors for up/down regulation
  theme_minimal() +
  labs(title = "Log2FoldChange for Selected Genes",
       x = "Gene",
       y = "Log2FoldChange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

### PCA plot ###

#  Extract normalized expression matrix
expression_matrix <- counts(dds, normalized = TRUE)

# Check the dimensions
dim(expression_matrix) # Should be genes (rows) x samples (columns)

# Filter lowly expressed genes (e.g., keep genes with counts > 10 in at least 10% of samples)
keep_genes <- rowSums(expression_matrix > 10) >= 0.1 * ncol(expression_matrix)
expression_matrix <- expression_matrix[keep_genes, ]

# Log-transform the data (add a pseudocount to avoid log(0))
expression_matrix <- log2(expression_matrix + 1)

# Transpose the matrix for PCA
expression_matrix_t <- t(expression_matrix)

# Perform PCA
pca_result <- prcomp(expression_matrix_t, scale. = TRUE)

# View PCA results
summary(pca_result)
head(pca_result)

# Create a data frame for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)

# Calculate variance explained
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Plot PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Breast Cancer Gene Expression Data",
       x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)"))

# Extract clinical data
clinical_data <- colData(combined_data)

# Add clinical data to the PCA data frame
pca_df$Group <- clinical_data$sample_type # Example: 'sample_type' column

# Plot with colors
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Breast Cancer Gene Expression Data",
       x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)"))

### Proportion of variance explained by each Principal Component ###

# Create a data frame for plotting
variance_df <- data.frame(
  PC = 1:length(variance_explained),
  Variance = variance_explained
)

# Scree plot
ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Scree Plot: Proportion of Variance Explained by Each PC",
       x = "Principal Component (PC)",
       y = "Proportion of Variance Explained")

### Loadings Plot ###

library(ggrepel)

# Variance Stabilizing Transformation (for better PCA performance)
vsd <- vst(dds, blind = TRUE)

# Extract transformed count matrix
pca_data <- assay(vsd)

# Perform PCA
pca_result <- prcomp(t(pca_data), scale. = TRUE)

# Extract loadings
loadings <- as.data.frame(pca_result$rotation)  # Rotation matrix contains loadings

# Add gene names for easier interpretation
loadings$gene_id <- rownames(loadings)

# Merge with gene names
loadings_annotated <- merge(loadings, gene_mapping, by = "gene_id", all.x = TRUE)

# Sort by absolute contribution to PC1 or PC2
top_genes <- loadings_annotated %>%
  arrange(desc(abs(PC1))) %>%
  head(30)  # Select top 30 contributing genes

ggplot(top_genes, aes(x = PC1, y = PC2, label = gene_name)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_text_repel(size = 4, max.overlaps = 20) +  # Avoid label overlap
  theme_minimal() +
  ggtitle("PCA Loading Plot: Top Genes") +
  xlab("PC1 Loadings") +
  ylab("PC2 Loadings")

### Adding Phylostrata column to the Raw Data ###

phylostrata <- read.table("phylostratum_database/gene_phylostrata.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(phylostrata)
# Renaming GeneID column name to gene_name
phylostrata <- phylostrata %>%
  mutate(gene_name = GeneID) %>%
  select(gene_name, everything()) %>%
  select(-GeneID)

# Merging raw data & phylostrata by gene name
all_genes <- all_genes %>%
  left_join(phylostrata, by = "gene_name")
head(all_genes)

### Visualising the data ### Post Phylostrata
library(tidyverse)

numeric_data <- all_genes %>%
  select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  drop_na()  # remove rows with NAs if needed

pca <- prcomp(numeric_data, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  mutate(
    gene_name = all_genes$gene_name[!is.na(all_genes$padj)],
    phylostrata = all_genes$Phylostrata[!is.na(all_genes$padj)]
  )

ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(phylostrata))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA of DE Stats Colored by Phylostrata",
    color = "Phylostrata"
  ) +
  theme_minimal()

# With NAs removed before PCA
all_genes_filtered <- all_genes %>%
  filter(!is.na(Phylostrata))

numeric_data <- all_genes_filtered %>%
  select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  drop_na()  # remove rows with NAs if needed

pca <- prcomp(numeric_data, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  mutate(
    gene_name = all_genes_filtered$gene_name[!is.na(all_genes_filtered$padj)],
    phylostrata = all_genes_filtered$Phylostrata[!is.na(all_genes_filtered$padj)]
  )

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(phylostrata))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA of Differential Expression results for Breast Cancer samples",
    color = "Phylostrata"
  ) +
  labs(color = "Phylostrata") +
  theme_classic()

ggsave("breast_cancer_pca.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

### Subsetted for Phyostrata 1-4 ###
pca_df_subset <- pca_df %>%
  filter(phylostrata %in% 1:4)

pca_plot_2 <- ggplot(pca_df_subset, aes(x = PC1, y = PC2, color = as.factor(phylostrata))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA of Differential Expression results for Breast Cancer samples - Only phylostrata 1-4",
    color = "Phylostrata"
  ) +
  labs(color = "Phylostrata") +
  theme_classic()

ggsave("breast_cancer_pca_subset.png", plot = pca_plot_2, width = 8, height = 6, dpi = 300)

### Opposite subset ###
pca_df_subset_2 <- pca_df %>%
  filter(phylostrata %in% 10:11)

pca_plot_3 <- ggplot(pca_df_subset_2, aes(x = PC1, y = PC2, color = as.factor(phylostrata))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA of Differential Expression results for Breast Cancer samples - Only phylostrata 10-11",
    color = "Phylostrata"
  ) +
  labs(color = "Phylostrata") +
  theme_classic()

ggsave("breast_cancer_pca_subset_2.png", plot = pca_plot_3, width = 8, height = 6, dpi = 300)

# Violin plots
all_genes_filtered_clean <- all_genes_filtered %>%
  filter(!is.na(log2FoldChange), !is.na(Phylostrata))

violin_plot <- ggplot(all_genes_filtered_clean, aes(x = as.factor(Phylostrata), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_violin(trim = FALSE, scale = "width", color = "gray30", alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Breast Cancer Log2FoldChange Results by Phylostrata",
    x = "Phylostrata",
    y = "log2(Fold Change)",
    fill = "Phylostrata"
  ) +
  theme(legend.position = "none")

ggsave("breast_cancer_violin.png", plot = violin_plot, width = 8, height = 6, dpi = 300)

### ANOVA across Phylostrata ###
kruskal.test(log2FoldChange ~ as.factor(Phylostrata), data = all_genes_filtered)

