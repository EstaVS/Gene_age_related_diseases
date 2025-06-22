setwd("~/Desktop/gene_age_project/")

library(tidyr)

### Download data ###

library(TCGAbiolinks)

brain_query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
) # Construct a query to fetch data based on specific filters

GDCdownload(brain_query) # Download the data

brain_data <- GDCprepare(brain_query) # Load the data in R

### Process data ###

# Load the package
library(SummarizedExperiment)

# View sample metadata
metadata <- colData(brain_data)
head(metadata)

# Check the unique values in the 'sample_type' column
table(metadata$sample_type)

# Subset both tumor and normal tissue samples from TCGA
tcga_subset <- brain_data[, brain_data$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]

# Check result
table(colData(tcga_subset)$sample_type)

# Get raw counts for Differential Expression
tcga_counts <- assay(tcga_subset, "unstranded")

### Load Libraries ###
library(dplyr)
library(readr)
library(data.table)

### Load and Process GTEx Cortex Data ###
gtex_cortex <- fread("GTEx_data/gene_reads_v10_brain_cortex.gct.gz", skip = 2)
gtex_cortex$Name <- sub("\\..*$", "", gtex_cortex$Name)  # Remove Ensembl ID version numbers
gtex_clean <- gtex_cortex[, lapply(.SD, sum), by = Name, .SDcols = -(1:2)]  # Sum duplicates
gtex_counts <- as.data.frame(gtex_clean)
rownames(gtex_counts) <- gtex_counts$Name
gtex_counts <- gtex_counts[, -which(names(gtex_counts) == "Name")]

### Load and Process TCGA Counts ###
rownames(tcga_counts) <- sub("\\..*$", "", rownames(tcga_counts))  # Remove Ensembl ID version numbers

### Combine GTEx and TCGA ###
common_genes <- intersect(rownames(gtex_counts), rownames(tcga_counts)) # By Ensembl ID
gtex_mat <- gtex_counts[common_genes, ]
tcga_mat <- tcga_counts[common_genes, ]
combined_data <- cbind(tcga_mat, gtex_mat)

### Create Sample Metadata ###

### --- SAMPLE METADATA (TCGA + GTEx) --- ###

# TCGA metadata
tcga_metadata <- colData(tcga_subset) %>%
  as.data.frame() %>%
  mutate(
    sample_type = sample_type,
    source = "TCGA",
    type = ifelse(sample_type == "Primary Tumor", "Tumor", "Normal"),
    batch = "TCGA"
  ) %>%
  select(sample_type, source, type, batch)

# GTEx metadata
gtex_metadata <- data.frame(
  sample_type = rep("Normal", ncol(gtex_mat)),
  source = "GTEx",
  type = "Normal",
  batch = "GTEx",
  row.names = colnames(gtex_mat)
)

# Combine metadata and subset to used columns
metadata <- rbind(
  tcga_metadata[colnames(tcga_mat), ],
  gtex_metadata
)

### --- GENE-LEVEL ANNOTATION WITH PHYLOSTRATA --- ###

# Read phylostrata table
phylostrata <- read.table("phylostrata/phylostrata_ensembl.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Match to genes in combined expression matrix
phylostrata_matched <- phylostrata[match(rownames(combined_data), phylostrata$EnsemblID), ]

# Merge with expression matrix
combined_data_phy <- cbind(Phylostrata = phylostrata_matched$Phylostrata, combined_data)

# Filter out rows with missing phylostrata
filtered_data <- combined_data_phy[!is.na(combined_data_phy$Phylostrata), ]

# Extract count matrix (remove Phylostrata column)
count_data <- filtered_data[, -which(colnames(filtered_data) == "Phylostrata")]

# Align metadata to count matrix
metadata <- metadata[colnames(count_data), ]

### --- GENE NAME MAPPING --- ###

gene_mapping <- phylostrata %>%
  filter(EnsemblID %in% rownames(filtered_data)) %>%
  select(EnsemblID, GeneID) %>%
  rename(gene_id = EnsemblID, gene_name = GeneID)

### --- DESEQ2 ANALYSIS --- ###

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = metadata,
  design = ~ batch + type
)

# Removing low count genes
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("type", "Tumor", "Normal")) %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(.))

# Filter significant DE genes
significant_genes <- res %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 2)

# Annotate with gene names
res_annotated <- significant_genes %>%
  left_join(gene_mapping, by = "gene_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
  relocate(gene_name)

# Save full DE results
all_genes <- res %>%
  left_join(gene_mapping, by = "gene_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
  relocate(gene_name)

write_csv(all_genes, "Brain_DE.csv")

### Visualize Top DE Genes ###
library(ggplot2)
 
# Function to plot top DE genes with varying n
plot_top_de_genes <- function(n, res_annotated, phylostrata) {
  top_genes <- res_annotated %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = n) %>%
    left_join(phylostrata, by = c("gene_id" = "EnsemblID"))
  
  ggplot(top_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange,
                        fill = as.factor(Phylostrata))) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d(name = "Phylostrata", option = "C") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_flip() +
    labs(title = paste("Top", n, "Differentially Expressed Genes"),
         x = "Gene",
         y = "Log2 Fold Change")
}

plot_top_de_genes(20, res_annotated, phylostrata) # Top 20
ggsave("plots/Brain_poster/top_20_DE_genes.png", width = 8, height = 6, dpi = 300)

### PCA: Expression-Based ###

## Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Extract transformed expression matrix
vst_expr <- assay(vsd)

# Remove batch effects with limma
library(limma)
vst_corrected <- removeBatchEffect(vst_expr, batch = metadata$batch)

# PCA
pca_expr <- prcomp(t(vst_corrected), scale. = TRUE)

# Build PCA dataframe
pca_df <- as.data.frame(pca_expr$x) %>%
  mutate(
    Sample = colnames(vst_corrected),
    Group = metadata$type,
    Batch = metadata$batch
  )

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "PCA of Normalised Expression by Tissue Status",
    x = paste0("PC1 (", round(100 * summary(pca_expr)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_expr)$importance[2, 2], 1), "%)"),
    color = "Group"
  )
ggsave("plots/Brain_poster/PCA_norm_exp.png", width = 8, height = 6, dpi = 300)

### UMAP: Expression Based ###
library(umap)

# Run UMAP on vst-corrected expression
umap_result <- umap(t(vst_corrected))  # transpose: genes x samples -> samples x genes

# Build UMAP dataframe
umap_df <- as.data.frame(umap_result$layout) %>%
  mutate(
    Sample = colnames(vst_corrected),
    Group = metadata$type
  )

# Plot UMAP
ggplot(umap_df, aes(x = V1, y = V2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "UMAP of Normalised Expression",
       x = "UMAP1", y = "UMAP2") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("plots/Brain_poster/UMAP_norm_exp.png", width = 8, height = 6, dpi = 300)

### PCA: DE Stats with Phylostrata ###
# VST not needed here as it is based on my DE analysis which already factors in batch effects

# Merge phylostrata into all_genes
all_genes <- left_join(all_genes, phylostrata, by = c("gene_id" = "EnsemblID"))

all_genes_filtered <- all_genes %>%
  filter(!is.na(Phylostrata)) %>%
  drop_na(log2FoldChange, padj)

numeric_stats <- all_genes_filtered %>%
  select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

pca_de <- prcomp(numeric_stats, center = TRUE, scale. = TRUE)
pca_de_df <- as.data.frame(pca_de$x) %>%
  mutate(gene_name = all_genes_filtered$gene_name,
         phylostrata = all_genes_filtered$Phylostrata)

pca_var <- pca_de$sdev^2  # variance per principal component
pca_var_explained <- round(100 * pca_var / sum(pca_var), 1)  # percentage

ggplot(pca_de_df, aes(x = PC1, y = PC2, color = as.factor(phylostrata))) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "PCA of DE Statistics by Phylostrata",
    x = paste0("PC1 (", pca_var_explained[1], "%)"),
    y = paste0("PC2 (", pca_var_explained[2], "%)"),
    color = "Phylostrata"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("plots/Brain_poster/PCA_DE.png", width = 8, height = 6, dpi = 300)

### Subset PCA by Phylostrata ###
library(patchwork)
library(scales)

# Set color palette for phylostrata 1–16
phylostrata_colors <- hue_pal()(16)
names(phylostrata_colors) <- as.character(1:16)

# Get unique phylostrata
unique_phylostrata <- sort(unique(pca_de_df$phylostrata))

# Global axis limits for PC1 and PC2 - keep axes the same for all plots
x_limits <- range(pca_de_df$PC1, na.rm = TRUE)
y_limits <- range(pca_de_df$PC2, na.rm = TRUE)

# Store plots
plot_list <- list()

for (ps in unique_phylostrata) {
  subset_df <- filter(pca_de_df, phylostrata == ps)
  
  p <- ggplot(subset_df, aes(x = PC1, y = PC2, color = as.factor(phylostrata))) +
    geom_point(alpha = 0.7, size = 2) +
    labs(title = paste("Phylostrata", ps),
         x = paste0("PC1 (", pca_var_explained[1], "%)"),
         y = paste0("PC2 (", pca_var_explained[2], "%)"),
         ) +
    scale_color_manual(values = phylostrata_colors) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  plot_list[[as.character(ps)]] <- p
}

# Combine in a grid with a title
combined_plot <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = "PCA of DE Statistics by Phylostrata",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
ggsave("plots/Brain_poster/PCA_DE_all_phylostrata_colored.png",
       plot = combined_plot, width = 16, height = 12, dpi = 300)

### Violin Plot of Log2FC by Phylostrata ###
ggplot(all_genes_filtered, aes(x = as.factor(Phylostrata), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_violin(trim = FALSE, scale = "width", color = "gray30", alpha = 0.8) +
  theme_classic() +
  labs(title = "Log2 Fold Change by Phylostrata",
       x = "Phylostrata",
       y = "Log2(Fold Change)") +
  theme(legend.position = "none")

ggsave("brain_cancer_violin.png", width = 8, height = 6, dpi = 300)

### Kruskal-Wallis Test ###
# To identify if there is AT LEAST one phylostrata that is significantly differentially expressed

# Kruskal-Willis, non-parametric alternative to ANOVA
kruskal.test(log2FoldChange ~ as.factor(Phylostrata), data = all_genes_filtered)

# Exploring phylostrata relationships Pairwise-Wilcox test
pw_result <- pairwise.wilcox.test(
  all_genes_filtered$log2FoldChange,
  as.factor(all_genes_filtered$Phylostrata),
  p.adjust.method = "BH"
)

pmat <- pw_result$p.value
# Fill lower triangle to make it symmetric
pw_df <- pmat # Redundant but kept in for now due to time constraints for editting
pw_df[lower.tri(pw_df)] <- t(pw_df)[lower.tri(pw_df)]

library(reshape2)
melted <- reshape2::melt(pmat, na.rm = TRUE)
colnames(melted) <- c("Phylostrata1", "Phylostrata2", "p_adj")

# Changing Phylostrata to factors for plotting
melted$Phylostrata1 <- factor(melted$Phylostrata1, levels = sort(unique(as.numeric(melted$Phylostrata1))))
melted$Phylostrata2 <- factor(melted$Phylostrata2, levels = sort(unique(as.numeric(melted$Phylostrata2))))

ggplot(melted, aes(x = Phylostrata1, y = Phylostrata2, fill = p_adj)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(trans = "log10", name = "Adjusted p-value") +
  theme_minimal() +
  labs(title = "Pairwise Wilcoxon Test (BH-adjusted)",
       x = "Phylostrata", y = "Phylostrata") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    axis.text.y = element_text(angle = 0),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("plots/Brain_poster/heatmap.png", width = 8, height = 6, dpi = 300)

# Testing if there's a bias towards specific phylostrata for the top DE genes

# Contingency table construction
# Define top N genes
top_n <- 20

# Get top N genes with phylostrata
top_genes <- res_annotated %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n) %>%
  left_join(phylostrata, by = c("gene_id" = "EnsemblID"))

# All genes with phylostrata
all_genes_with_ps <- res_annotated %>%
  left_join(phylostrata, by = c("gene_id" = "EnsemblID")) %>%
  filter(!is.na(Phylostrata))

# Count phylostrata in top and background
top_counts <- table(top_genes$Phylostrata)
background_counts <- table(all_genes_with_ps$Phylostrata)

# Ensure both tables have the same factor levels
all_levels <- union(names(top_counts), names(background_counts))
top_counts <- top_counts[all_levels]
background_counts <- background_counts[all_levels]
top_counts[is.na(top_counts)] <- 0
background_counts[is.na(background_counts)] <- 0

# Fisher test
contingency <- rbind(
  Top = top_counts,
  Background = background_counts - top_counts
)

# Run test
fisher.test(contingency, simulate.p.value = TRUE, B = 10000)
