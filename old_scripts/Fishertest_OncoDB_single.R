setwd("~/Desktop/gene_age_project/")

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(pheatmap)
library(umap)
library(scales)
library(tidyr)

# Load phylostrata mapping
phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID)

# Load project info
project_info <- read_csv("metadata/project_info.csv")

# Runs all OncoDB projects file-by-file
project_path <- "OncoDB_data/expression_diff_ACC.txt"
project_id <- str_extract(basename(project_path), "(?<=expression_diff_)[^\\.]+")
cancer_type_name <- project_info %>%
  filter(project_id == !!project_id) %>%
  pull(cancer_type)

cat("\nProcessing:", project_id, "-", cancer_type_name, "\n")

# Create a directory of the cancer for testing
dir.create(paste0("plots/OncoDB/", project_id), recursive = TRUE, showWarnings = FALSE)

# Read and clean data
data <- read.table(project_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")
data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))

data <- data %>% filter(gene_name %in% phylostrata$gene_name)
data <- data %>% left_join(phylostrata, by = "gene_name")
data_clean <- data %>% filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))

# Automatically find the optimal top_n that gives the lowest Fisher p-value
search_ns <- seq(20, 200, by = 10)
results <- data.frame(top_n = numeric(), p_value = numeric())

for (n in search_ns) {
  top_genes_test <- data_clean %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = n)
  
  top_counts_test <- table(factor(top_genes_test$Phylostrata, levels = 1:16))
  background_counts_test <- table(factor(data_clean$Phylostrata, levels = 1:16))
  contingency_test <- rbind(Top = top_counts_test, Background = background_counts_test - top_counts_test)
  
  p_val <- fisher.test(contingency_test, simulate.p.value = TRUE, B = 10000)$p.value
  results <- rbind(results, data.frame(top_n = n, p_value = p_val))
}

best_n <- results %>% arrange(p_value) %>% slice(1) %>% pull(top_n)
cat("\nBest top_n based on Fisher p-value:", best_n, "\n")

# Fisher enrichment test using best_n
top_n <- best_n
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

# Save to file
output_path <- paste0("stats/OncoDB_", project_id, "_fisher_test.txt")
capture.output(fisher_result, file = output_path)

# Top N Bar Plot - colored by Phylostrata
top_n_genes <- data_clean %>%
  filter(padj < 0.001) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(top_n)

top_DE_plot <- ggplot(top_n_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  scale_fill_viridis_d(name = "Phylostrata") +
  labs(title = paste0("Top ", top_n, " DE Genes (FDR < 0.001) in ", cancer_type_name, " (OncoDB)"),
       x = "Gene", y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 6), plot.title = element_text(size = 12))
show(top_DE_plot)

# Enrichment plot
phylostrata_levels <- sort(unique(data_clean$Phylostrata))
enrichment_pvals <- sapply(phylostrata_levels, function(ps) {
  in_top <- sum(top_n_genes$Phylostrata == ps)
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
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = Inf, y = -log10(0.05), label = "p = 0.05", hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
  scale_fill_viridis_d() +
  theme_bw() +
  labs(title = paste0("Phylostrata Enrichment (Fisher Test, BH-adjusted) in ", cancer_type_name, " (OncoDB)"),
       x = "Phylostrata", y = "-log10(adjusted p-value)") +
  theme(legend.position = "none")
show(enrichment_plot)

# Bootstrapped Downsampling Analysis
set.seed(42)
boot_df <- data.frame()
for (i in 1:100) {
  sampled <- data_clean %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_sample(n = top_n)
  
  for (ps in phylostrata_levels) {
    in_top <- sum(sampled$Phylostrata == ps)
    in_bg <- sum(data_clean$Phylostrata == ps) - in_top
    not_in_ps_top <- top_n - in_top
    not_in_ps_bg <- nrow(data_clean) - top_n - in_bg
    mat <- matrix(c(in_top, in_bg, not_in_ps_top, not_in_ps_bg), nrow = 2)
    test <- suppressWarnings(fisher.test(mat))
    boot_df <- rbind(boot_df, data.frame(Boot = i, Phylostrata = ps, OR = test$estimate, pval = test$p.value))
  }
}

boot_summary <- boot_df %>%
  group_by(Phylostrata) %>%
  summarise(median_OR = median(OR, na.rm = TRUE), prop_sig = mean(pval < 0.05))

boot_plot <- ggplot(boot_summary, aes(x = as.factor(Phylostrata), y = median_OR, fill = prop_sig)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f\n(%.2f)", median_OR, prop_sig)), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_viridis_c(name = "Proportion p < 0.05") +
  theme_bw() +
  labs(title = paste0("Bootstrapped Phylostrata Enrichment in ", cancer_type_name, " (OncoDB)"),
       x = "Phylostrata", y = "Median Odds Ratio")
show(boot_plot)

ggsave(paste0("plots/OncoDB/", project_id, "/OncoDB_", project_id, "_top_", top_n, "_phylostrata.png"), plot = top_DE_plot, width = 8, height = 6, dpi = 300)
ggsave(paste0("plots/OncoDB/", project_id, "/OncoDB_", project_id, "_enrichment.png"), plot = enrichment_plot, width = 8, height = 6, dpi = 300)
ggsave(paste0("plots/OncoDB/", project_id, "/OncoDB_", project_id, "_bootstrap.png"), plot = boot_plot, width = 8, height = 6, dpi = 300)
