


setwd("~/Downloads/gene_age_project/raw_data/expression/")

library(dplyr) 

BLCA <- read.table("BLCA_Differential_Gene_Expression_Table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BLCA

library(ggplot2)

BLCA$Gene.symbol <- as.factor(BLCA$Gene.symbol)  # Convert to factor

# Keep only top 20 genes with highest absolute fold change
BLCA_top <- BLCA %>% 
  arrange(desc(abs(log2.fold.change))) %>% 
  head(20)

ggplot(BLCA_top, aes(x = reorder(Gene.symbol, log2.fold.change), y = log2.fold.change, fill = Gene.symbol)) +
  geom_col() +  # Works better for negative values
  coord_flip() +  # Flip for better readability
  #scale_fill_brewer(palette = "Set3") +  # Uses a good categorical palette - doesn't work
  labs(title = "Top 20 Gene Expression Changes",
       x = "Gene",
       y = "log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

phylostrata <- read.table("~/Downloads/gene_age_project/phylostratum_database/gene_phylostrata.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
phylostrata

### Making new dataframe

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Set the directory containing the files
data_dir <- "/Users/estashrewsbury/Downloads/gene_age_project/raw_data/expression/"

# List all files in the directory
file_list <- list.files(path = data_dir, pattern = "*.txt", full.names = TRUE)

# Initialize an empty list to store dataframes
df_list <- list()

# Loop through each file and process it
for (file in file_list) {
  # Extract cancer type from filename (modify as needed based on your file naming)
  cancer_type <- str_extract(basename(file), "^[A-Za-z0-9]+")
  
  # Read the file
  df <- read_delim(file, delim = "\t")
  
  # Keep only relevant columns
  df <- df %>% select(`Gene symbol`, `log2 fold change`)
  
  # Rename the log2FoldChange column to the cancer type
  colnames(df)[2] <- cancer_type
  
  # Store in the list
  df_list[[cancer_type]] <- df
}

# Merge all dataframes by Gene symbol
df_combined <- Reduce(function(x, y) full_join(x, y, by = "Gene symbol"), df_list)

# Convert tibble to a data frame
df_combined <- as.data.frame(df_combined)

# Verify that "Gene symbol" exists
if("Gene symbol" %in% colnames(df_combined)) {
  # Set row names using 'Gene symbol'
  rownames(df_combined) <- df_combined$`Gene symbol`
  
  # Remove 'Gene symbol' column since it's now the row names
  df_combined <- df_combined[, !(colnames(df_combined) %in% "Gene symbol")]
} else {
  stop("Error: 'Gene symbol' column not found. Check column names with colnames(df_combined).")
}

# Check if row names are set correctly
head(df_combined)

### Check if gene symbols in df match phylostrata



# Extract gene symbols from Phylostrata file
phylostrata_genes <- unique(phylostrata$`Gene symbol`)  # Update column name if needed

# Extract row names from df_combined (our expression data)
expression_genes <- rownames(df_combined)

# Find genes that are in both datasets
common_genes <- intersect(expression_genes, phylostrata_genes)

# Find genes that are missing in the Phylostrata file
missing_in_phylostrata <- setdiff(expression_genes, phylostrata_genes)

# Find genes that are in Phylostrata but not in the expression dataset
missing_in_expression <- setdiff(phylostrata_genes, expression_genes)

# Print summary
cat("Number of matching genes:", length(common_genes), "\n")
cat("Number of genes in expression data but not in Phylostrata:", length(missing_in_phylostrata), "\n")
cat("Number of genes in Phylostrata but not in expression data:", length(missing_in_expression), "\n")

# View the missing genes
missing_in_phylostrata

### Ignore for now

### Violin plots for Log2FoldChange

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Convert row names (genes) into a column
df_combined$Gene <- rownames(df_combined)

# Ensure Gene column is treated as character
df_combined$Gene <- as.character(df_combined$Gene)

# Compute the top 20 genes with the highest absolute log2FoldChange across any cancer type
top_genes <- df_combined %>%
  mutate(max_log2FC = apply(select(., -Gene), 1, function(x) max(abs(x), na.rm = TRUE))) %>%
  arrange(desc(max_log2FC)) %>%
  slice(1:20) %>%
  pull(Gene)  # Extract top 20 gene names

# Subset the dataset to include only the top 10 genes
df_top <- df_combined %>% filter(Gene %in% top_genes)

# Convert wide format to long format using pivot_longer (better than melt)
df_long <- df_top %>%
  pivot_longer(cols = -Gene, names_to = "Cancer_Type", values_to = "log2FoldChange") %>%
  drop_na()  # Remove missing values

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)

# Define a color palette with 50+ colors
palette_50 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(df_long$Gene)))

# Ensure Gene is a factor to apply colors correctly
df_long$Gene <- factor(df_long$Gene, levels = unique(df_long$Gene))

ggplot(df_long, aes(x = Gene, y = log2FoldChange, fill = Gene)) +
  geom_violin(trim = FALSE, width = 1, scale = "width", color = "black") +  # Adjust width and scaling
  geom_jitter(width = 0.15, alpha = 0.6, size = 1, color = "black") +  # Reduce point spread
  scale_fill_manual(values = palette_50) +  # Apply custom color palette
  labs(x = "Gene Symbol", y = "log2 Fold Change", title = "Violin Plots of log2 Fold Change for Many Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")  # Hide legend
