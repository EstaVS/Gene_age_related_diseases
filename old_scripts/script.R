### Setup ###
setwd("~/Desktop/gene_age_project")

library(SummarizedExperiment)
library(tidyverse)
library(stringr)

### Automating Naming ###

# Define cancer types and associated project IDs
project_info <- tribble(
  ~project_id, ~cancer_type,
  "ACC",  "Adrenocortical Carcinoma",
  "BLCA", "Bladder Urothelial Carcinoma",
  "BRCA", "Breast Invasive Carcinoma",
  "CESC", "Cervical Squamous Cell Carcinoma & Endocervical Adenocarcinoma",
  "CHOL", "Cholangiocarcinoma",
  "COAD", "Colon Adenocarcinoma",
  "GBM",  "Glioblastoma Multiforme",
  "HNSC", "Head & Neck Squamous Cell Carcinoma",
  "KICH", "Kidney Chromophobe",
  "KIRC", "Kidney Renal Clear Cell Carcinoma",
  "KIRP", "Kidney Renal Papillary Cell Carcinoma",
  "LGG",  "Low Grade Glioma",
  "LIHC", "Liver Hepatocellular Carcinoma",
  "LUAD", "Lung Adenocarcinoma",
  "LUSC", "Lung Squamous Cell Carcinoma",
  "OV",   "Ovarian Serous Cystadenocarcinoma",
  "PAAD", "Pancreatic Adenocarcinoma",
  "PCPG", "Pheochromocytoma & Paraganglioma",
  "PRAD", "Prostate Adenocarcinoma",
  "READ", "Rectum Adenocarcinoma",
  "SARC", "Sarcoma",
  "SKCM", "Skin Cutaneous Melanoma",
  "TGCT", "Testicular Germ Cell Tumors",
  "THCA", "Thyroid Carcinoma",
  "THYM", "Thymoma",
  "UCEC", "Uterine Corpus Endometrial Carcinoma",
  "UCS",  "Uterine Carcinosarcoma"
)

# Define the file you want to process
project_path <- "OncoDB_data/expression_diff_BRCA.txt"

# Extract project ID from filename
project_id <- str_extract(project_path, "(?<=expression_diff_)[^\\.]+")

# Match cancer type
cancer_type_name <- project_info %>%
  filter(project_id == !!project_id) %>%
  pull(cancer_type)

# Preview
cat("Project ID:", project_id, "\nCancer Type:", cancer_type_name, "\n")

### Data Loading ##

data <- read.table(project_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID)

# Keep only genes that are found in Phylostrata
genes_in_phylostrata <- data %>%
  filter(Gene %in% phylostrata$gene_name)

# Filter data for only genes found in Phylostrata
data <- data[data$Gene %in% genes_in_phylostrata$Gene, ]

### Data Cleaning ###

# Changing column names
colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")

# Remove 'q<=' strings and convert to numeric
data$padj <- gsub("q<=|Q<=|Q=|q=", "", data$padj)
data$padj <- as.numeric(data$padj)

# Remove NAs
data_clean <- data %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# Check structure
str(data_clean)
summary(data_clean$padj)
summary(data_clean$log2FoldChange)

### Visualization: Top Genes (FDR < 0.001) ###
top_20_genes <- data %>%
  filter(!is.na(log2FoldChange), padj < 0.001) %>%  # Add FDR filtering
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)

top_20_genes <- top_20_genes %>%
  mutate(Expression = ifelse(log2FoldChange > 0, "Overexpressed", "Underexpressed"))

top_DE_plot <- ggplot(top_20_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = Expression)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Underexpressed" = "blue", "Overexpressed" = "red")) +
  coord_flip() +
  theme_minimal() +
  labs(title = paste0("Top 20 Differentially Expressed Genes (FDR < 0.001) in ", cancer_type_name," (OncoDB)"),
       x = "Gene",
       y = "Log2 Fold Change",
       fill = NULL) +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12)
    )
show(top_DE_plot)

### PCA ###

# Add Phylostrata info to Data
data <- data %>%
  left_join(phylostrata, by = "gene_name")

# Keep only numeric data
numeric_data <- data %>%
  select(where(is.numeric))  # keeps only numeric columns

# Now PCA
pca_phylostrata <- prcomp(numeric_data, center = TRUE, scale. = TRUE)
pca_phyl_df <- as.data.frame(pca_phylostrata$x) %>%
  mutate(gene_name = data$gene_name,
         phylostrata = as.factor(data$Phylostrata))

variance_explained <- pca_phylostrata$sdev^2 / sum(pca_phylostrata$sdev^2)

# PCA plot colored by Phylostrata
pca_phyl_plot <- ggplot(pca_phyl_df, aes(x = PC1, y = PC2, color = phylostrata)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(title = paste0("PCA of Differential Expression Results for ", cancer_type_name, " (OncoDB)"),
       x = paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)"),
       color = "Phylostrata")
show(pca_phyl_plot)

### Save plots ###
ggsave(paste0("plots/OncoDB/Barplot/OncoDB", project_id, "_top_20_genes.png"), plot = top_DE_plot, width = 8, height = 6)
ggsave(paste0("plots/OncoDB/PCA/OncoDB", project_id, "_pca_phylostrata.png"), plot = pca_phyl_plot, width = 8, height = 6)

### Violin Plot of log2FoldChange by Phylostrata (FDR < 0.001) ###

# Clean the data for plotting
data <- data %>%
  filter(!is.na(log2FoldChange), !is.na(Phylostrata), padj < 0.001)  # Add FDR filter

# Create the violin plot
violin_plot <- ggplot(data, aes(x = as.factor(Phylostrata), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_violin(trim = FALSE, scale = "width", color = "gray30", alpha = 0.8) +
  theme_classic() +
  labs(
    title = paste("Log2FoldChange Distribution by Phylostrata (FDR < 0.001) in", cancer_type_name, "(OncoDB)"),
    x = "Phylostrata",
    y = "Log2(Fold Change)",
    fill = "Phylostrata"
  ) +
  theme(legend.position = "none")
show(violin_plot)

# Save the plot
ggsave(paste0("plots/OncoDB/Violin/OncoDB", project_id, "_violin.png"), plot = violin_plot, width = 8, height = 6, dpi = 300)

### Kruskal-Wallis Test ###
kruskal_test_result <- kruskal.test(log2FoldChange ~ as.factor(Phylostrata), data = res_df_clean)
message("Kruskal-Wallis p-value: ", kruskal_test_result$p.value)