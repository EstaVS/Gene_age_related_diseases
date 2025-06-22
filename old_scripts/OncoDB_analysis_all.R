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
project_info <- read_csv("metadata/project_info.csv")  # assuming this exists and maps project_id to cancer_type

# Ensure output directories exist
dir.create("plots/OncoDB/Barplot", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/PCA", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/UMAP", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/Violin", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OncoDB/Wilcox", recursive = TRUE, showWarnings = FALSE)
dir.create("stats", showWarnings = FALSE, recursive = TRUE)

# List all OncoDB differential expression files
files <- list.files("OncoDB_data", pattern = "expression_diff_.*\\.txt", full.names = TRUE)

for (project_path in files) {
  project_id <- str_extract(basename(project_path), "(?<=expression_diff_)[^\\.]+")
  cancer_type_name <- project_info %>%
    filter(project_id == !!project_id) %>%
    pull(cancer_type)
  
  cat("\nProcessing:", project_id, "-", cancer_type_name, "\n")
  
  data <- read.table(project_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")
  data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))
  
  data <- data %>% filter(gene_name %in% phylostrata$gene_name)
  data <- data %>% left_join(phylostrata, by = "gene_name")
  data_clean <- data %>% filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))
  
  # Top 20 Bar Plot
  top_20_genes <- data_clean %>%
    filter(padj < 0.001) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(20) %>%
    mutate(Expression = ifelse(log2FoldChange > 0, "Overexpressed", "Underexpressed"))
  
  top_DE_plot <- ggplot(top_20_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = Expression)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Underexpressed" = "blue", "Overexpressed" = "red")) +
    coord_flip() +
    theme_bw() +
    labs(title = paste0("Top 20 DE Genes (FDR < 0.001) in ", cancer_type_name, " (OncoDB)"),
         x = "Gene", y = "Log2 Fold Change") +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 12))
  
  # PCA with refined inputs
  pca_input <- data_clean %>%
    mutate(log_padj = -log10(padj)) %>%
    select(log2FoldChange, log_padj)
  pca_result <- prcomp(scale(pca_input), center = TRUE)
  pca_df <- as.data.frame(pca_result$x) %>%
    mutate(phylostrata = as.factor(data_clean$Phylostrata))
  
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = phylostrata)) +
    geom_point(alpha = 0.7) +
    theme_bw() +
    labs(title = paste0("PCA of DE Genes in ", cancer_type_name, " (OncoDB)"))
  
  # UMAP with refined inputs
  umap_result <- umap(scale(pca_input))
  umap_df <- as.data.frame(umap_result$layout) %>%
    mutate(phylostrata = as.factor(data_clean$Phylostrata))
  
  umap_plot <- ggplot(umap_df, aes(x = V1, y = V2, color = phylostrata)) +
    geom_point(alpha = 0.7) +
    theme_bw() +
    labs(title = paste0("UMAP of DE Genes in ", cancer_type_name, " (OncoDB)"))
  
  # Violin Plot
  violin_plot <- data_clean %>%
    filter(padj < 0.001) %>%
    ggplot(aes(x = as.factor(Phylostrata), y = log2FoldChange, fill = as.factor(Phylostrata))) +
    geom_violin(trim = FALSE, scale = "width", color = "gray30", alpha = 0.8) +
    theme_bw() +
    labs(title = paste0("Log2FC by Phylostrata (FDR < 0.001) in ", cancer_type_name),
         x = "Phylostrata", y = "Log2(Fold Change)") +
    theme(legend.position = "none")
  
  # Statistical tests
  data_clean <- data_clean %>%
    mutate(direction = case_when(
      log2FoldChange >= 2 & padj < 0.05 ~ "Up",
      log2FoldChange <= -2 & padj < 0.05 ~ "Down",
      TRUE ~ "NS"
    ))
  
  fisher_table <- table(data_clean$direction, data_clean$Phylostrata)
  fisher_result <- tryCatch(fisher.test(fisher_table), error = function(e) NULL)
  
  wilcox_data <- data_clean %>% filter(padj < 0.05)
  wilcox_test <- pairwise.wilcox.test(
    wilcox_data$log2FoldChange,
    as.factor(wilcox_data$Phylostrata),
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  wilcox_matrix <- as.matrix(as.data.frame(wilcox_test$p.value))
  for (i in 1:nrow(wilcox_matrix)) {
    for (j in 1:ncol(wilcox_matrix)) {
      if (is.na(wilcox_matrix[i,j])) {
        wilcox_matrix[i,j] <- wilcox_matrix[j,i]
      }
    }
  }
  wilcox_matrix[is.na(wilcox_matrix)] <- 1
  rownames(wilcox_matrix) <- colnames(wilcox_matrix) <- colnames(wilcox_test$p.value)
  
  # Create ggplot-style heatmap from Wilcoxon matrix
  wilcox_df <- as.data.frame(as.table(wilcox_matrix))
  colnames(wilcox_df) <- c("Phylostrata1", "Phylostrata2", "p_adj")
  wilcox_df$Phylostrata1 <- as.numeric(as.character(wilcox_df$Phylostrata1))
  wilcox_df$Phylostrata2 <- as.numeric(as.character(wilcox_df$Phylostrata2))
  wilcox_df$p_adj <- as.numeric(wilcox_df$p_adj)
  
  gg_wilcox <- ggplot(wilcox_df, aes(x = Phylostrata2, y = Phylostrata1, fill = p_adj)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("purple", "darkblue", "green", "yellow"),
                         trans = "log10",
                         name = "Adjusted p-value",
                         breaks = c(1e-9, 1e-6, 1e-3),
                         labels = scales::trans_format("log10", math_format(10^.x))) +
    coord_fixed() +
    theme_bw() +
    labs(title = "Pairwise Wilcoxon Test (BH-adjusted)",
         x = "Phylostrata", y = "Phylostrata") +
    theme(axis.text = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"))
  
  # Save plots
  ggsave(paste0("plots/OncoDB/Barplot/OncoDB_", project_id, "_top_20_genes.png"), plot = top_DE_plot, width = 8, height = 6)
  ggsave(paste0("plots/OncoDB/PCA/OncoDB_", project_id, "_pca.png"), plot = pca_plot, width = 8, height = 6)
  ggsave(paste0("plots/OncoDB/UMAP/OncoDB_", project_id, "_umap.png"), plot = umap_plot, width = 8, height = 6)
  ggsave(paste0("plots/OncoDB/Violin/OncoDB_", project_id, "_violin.png"), plot = violin_plot, width = 8, height = 6)
  ggsave(paste0("plots/OncoDB/Wilcox/OncoDB_", project_id, "_wilcox_heatmap.png"), plot = gg_wilcox, width = 8, height = 6)
  
  # Save stats
  if (!is.null(fisher_result)) {
    capture.output(fisher_result, file = paste0("stats/OncoDB_", project_id, "_fisher.txt"))
  }
  capture.output(wilcox_test, file = paste0("stats/OncoDB_", project_id, "_wilcox.txt"))
}
