setwd("~/Desktop/gene_age_project/")

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(pheatmap)
library(umap)
library(scales)

# Load phylostrata mapping
phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID)

# Load project info
project_info <- read_csv("metadata/project_info.csv")

# Ensure output directories exist
dir.create("plots/OncoDB/Barplot", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/PCA", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/UMAP", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/Violin", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/Wilcox", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/Enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("stats", showWarnings = FALSE, recursive = TRUE)

# Choose one OncoDB differential expression file to run manually
project_path <- "OncoDB_data/expression_diff_BRCA.txt"  # Change this to your desired project
project_id <- str_extract(basename(project_path), "(?<=expression_diff_)[^\\.]+")
cancer_type_name <- project_info %>%
  filter(project_id == !!project_id) %>%
  pull(cancer_type)

cat("\nProcessing:", project_id, "-", cancer_type_name, "\n")

# Read and clean data
data <- read.table(project_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")
data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))

data <- data %>% filter(gene_name %in% phylostrata$gene_name)
data <- data %>% left_join(phylostrata, by = "gene_name")
data_clean <- data %>% filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))

# Fisher enrichment test on top 120
top_n <- 120
top_genes <- data_clean %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n)

top_counts <- table(factor(top_genes$Phylostrata, levels = 1:16))
background_counts <- table(factor(data_clean$Phylostrata, levels = 1:16))

contingency <- rbind(
  Top = top_counts,
  Background = background_counts - top_counts
)

fisher_result <- fisher.test(contingency, simulate.p.value = TRUE, B = 10000)
print(fisher_result)

# Top 100 Bar Plot - colored by Phylostrata
top_100_genes <- data_clean %>%
  filter(padj < 0.001) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(100)

top_DE_plot <- ggplot(top_100_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  scale_fill_viridis_d(name = "Phylostrata") +
  labs(title = paste0("Top 100 DE Genes (FDR < 0.001) in ", cancer_type_name, " (OncoDB)"),
       x = "Gene", y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 6), plot.title = element_text(size = 12))
show(top_DE_plot)

ggsave(paste0("plots/OncoDB/Barplot/OncoDB_", project_id, "_top_100_phylostrata.png"), plot = top_DE_plot, width = 8, height = 6)

# Enrichment plot
phylostrata_levels <- sort(unique(data_clean$Phylostrata))
enrichment_pvals <- sapply(phylostrata_levels, function(ps) {
  in_top <- sum(top_100_genes$Phylostrata == ps)
  in_bg <- sum(data_clean$Phylostrata == ps) - in_top
  not_in_ps_top <- top_n - in_top
  not_in_ps_bg <- nrow(data_clean) - top_n - in_bg
  matrix <- matrix(c(in_top, in_bg, not_in_ps_top, not_in_ps_bg), nrow = 2)
  fisher.test(matrix)$p.value
})

adjusted_pvals <- p.adjust(enrichment_pvals, method = "BH")
enrichment_df <- data.frame(
  Phylostrata = phylostrata_levels,
  p_adj = adjusted_pvals
)

enrichment_plot <- ggplot(enrichment_df, aes(x = as.factor(Phylostrata), y = -log10(p_adj), fill = as.factor(Phylostrata))) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  theme_bw() +
  labs(title = "Phylostrata Enrichment (Fisher Test, BH-adjusted)",
       x = "Phylostrata", y = "-log10(adjusted p-value)") +
  theme(legend.position = "none")

show(enrichment_plot)

# PCA
pca_input <- data_clean %>% mutate(log_padj = -log10(padj)) %>% select(log2FoldChange, log_padj)
pca_result <- prcomp(scale(pca_input), center = TRUE)
pca_df <- as.data.frame(pca_result$x) %>%
  mutate(phylostrata = as.factor(data_clean$Phylostrata))

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = phylostrata)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = paste0("PCA of DE Genes in ", cancer_type_name, " (OncoDB)"))
show(pca_plot)

# UMAP
umap_result <- umap(scale(pca_input))
umap_df <- as.data.frame(umap_result$layout) %>%
  mutate(phylostrata = as.factor(data_clean$Phylostrata))

umap_plot <- ggplot(umap_df, aes(x = V1, y = V2, color = phylostrata)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = paste0("UMAP of DE Genes in ", cancer_type_name, " (OncoDB)"))
show(umap_plot)

# Violin Plot
violin_plot <- data_clean %>%
  filter(padj < 0.001) %>%
  ggplot(aes(x = as.factor(Phylostrata), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_violin(trim = FALSE, scale = "width", color = "gray30", alpha = 0.8) +
  theme_bw() +
  labs(title = paste0("Log2FC by Phylostrata (FDR < 0.001) in ", cancer_type_name),
       x = "Phylostrata", y = "Log2(Fold Change)") +
  theme(legend.position = "none")
show(violin_plot)

# Pairwise Wilcoxon Test with filled matrix
all_genes_filtered <- data_clean %>% filter(padj < 0.05)
pw_result <- pairwise.wilcox.test(
  all_genes_filtered$log2FoldChange,
  as.factor(all_genes_filtered$Phylostrata),
  p.adjust.method = "BH"
)

pmat <- pw_result$p.value
pw_df <- pmat
pw_df[lower.tri(pw_df)] <- t(pmat)[lower.tri(pmat)]

melted <- reshape2::melt(pmat, na.rm = TRUE)
colnames(melted) <- c("Phylostrata1", "Phylostrata2", "p_adj")
melted$Phylostrata1 <- factor(melted$Phylostrata1, levels = sort(unique(as.numeric(melted$Phylostrata1))))
melted$Phylostrata2 <- factor(melted$Phylostrata2, levels = sort(unique(as.numeric(melted$Phylostrata2))))

heatmap_plot <- ggplot(melted, aes(x = Phylostrata1, y = Phylostrata2, fill = p_adj)) +
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
show(heatmap_plot)

# Save plots
ggsave(paste0("plots/OncoDB/Barplot/OncoDB_", project_id, "_top_100_genes.png"), plot = top_DE_plot, width = 8, height = 6)
ggsave(paste0("plots/OncoDB/PCA/OncoDB_", project_id, "_pca.png"), plot = pca_plot, width = 8, height = 6)
ggsave(paste0("plots/OncoDB/UMAP/OncoDB_", project_id, "_umap.png"), plot = umap_plot, width = 8, height = 6)
ggsave(paste0("plots/OncoDB/Violin/OncoDB_", project_id, "_violin.png"), plot = violin_plot, width = 8, height = 6)
ggsave(paste0("plots/OncoDB/Wilcox/OncoDB_", project_id, "_wilcox_heatmap.png"), plot = heatmap_plot, width = 8, height = 6)
ggsave(paste0("plots/OncoDB/Enrichment/OncoDB_", project_id, "_phylostrata_enrichment.png"), plot = enrichment_plot, width = 8, height = 6)

# Save stats
if (!is.null(fisher_result)) {
  capture.output(fisher_result, file = paste0("stats/OncoDB_", project_id, "_fisher.txt"))
}
capture.output(pw_result, file = paste0("stats/OncoDB_", project_id, "_wilcox.txt"))
