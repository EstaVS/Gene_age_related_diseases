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
library(metap)
library(scales)
library(viridis)

# Load phylostrata mapping
phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID)

# Load project info
project_info <- read_csv("metadata/project_info.csv")

# Loop through all expression files
files <- list.files("OncoDB_data", pattern = "expression_diff_.*\\.txt", full.names = TRUE)

for (project_path in files) {
  project_id <- str_extract(basename(project_path), "(?<=expression_diff_)[^\\.]+")
  cancer_type_name <- project_info %>%
    filter(project_id == !!project_id) %>%
    pull(cancer_type)
  
  cat("\nProcessing:", project_id, "-", cancer_type_name, "\n")
  
  dir.create(paste0("plots/OncoDB/", project_id), recursive = TRUE, showWarnings = FALSE)
  
  # Read and clean data
  data <- read.table(project_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")
  data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))
  
  data <- data %>% filter(gene_name %in% phylostrata$gene_name)
  data <- data %>% left_join(phylostrata, by = "gene_name")
  data_clean <- data %>% filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))
  
  # Automatically find best top_n (starting from 1000 down to 20)
  search_ns <- seq(1000, 20, by = -10)
  best_n <- NA
  
  for (n in search_ns) {
    top_genes_test <- data_clean %>%
      arrange(desc(abs(log2FoldChange))) %>%
      slice_head(n = n)
    
    top_counts_test <- table(factor(top_genes_test$Phylostrata, levels = 1:16))
    background_counts_test <- table(factor(data_clean$Phylostrata, levels = 1:16))
    
    contingency_test <- rbind(
      Top = top_counts_test,
      Background = background_counts_test - top_counts_test
    )
    
    if (all(rowSums(contingency_test) > 0)) {
      p_val <- fisher.test(contingency_test, simulate.p.value = TRUE, B = 10000)$p.value
    } else {
      p_val <- NA
    }
    
    # Stop when we first find a significant value
    if (!is.na(p_val) && p_val < 0.05) {
      best_n <- n
      break
    }
  }
  
  # Fallback n = 20 in case none were significant
  if (is.na(best_n)) {
    best_n <- 20
  }
  
  top_n <- best_n
  cat("\nSelected top_n based on Fisher p < 0.05:", top_n, "\n")
  
  # Fisher test and save
  top_genes <- data_clean %>% arrange(desc(abs(log2FoldChange))) %>% slice_head(n = top_n)
  top_counts <- table(factor(top_genes$Phylostrata, levels = 1:16))
  background_counts <- table(factor(data_clean$Phylostrata, levels = 1:16))
  contingency <- rbind(Top = top_counts, Background = background_counts - top_counts)
  fisher_result <- fisher.test(contingency, simulate.p.value = TRUE, B = 10000)

  # Save Fisher test result
  capture.output(fisher_result, file = paste0("stats/OncoDB/fisher/OncoDB_", project_id, "_fisher_test.txt"))
  
  # Append the top_n used
  cat(paste0("Top N genes used: ", top_n, "\n"),
      file = paste0("stats/OncoDB/fisher/OncoDB_", project_id, "_fisher_test.txt"),
      append = TRUE)
  # Barplot of DE genes
  top_n_genes <- data_clean %>%
    filter(padj < 0.001) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(top_n)
  
  # Adjust label size based on number of top genes
  y_label_size <- if (top_n > 140) {
    4
  } else if (top_n > 100) {
    6
  } else {
    8
  }
  
  top_DE_plot <- ggplot(top_n_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = as.factor(Phylostrata))) +
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    scale_fill_viridis_d(name = "Phylostrata") +
    labs(title = paste0("Top ", top_n, " DE Genes (FDR < 0.001) in ", cancer_type_name," (OncoDB)"),
         x = "Gene", y = "Log2 Fold Change") +
    theme(axis.text.y = element_text(size = y_label_size), plot.title = element_text(size = 12))
  ggsave(paste0("plots/OncoDB/", project_id, "/OncoDB_", project_id, "_top_", top_n, "_phylostrata.png"), plot = top_DE_plot, width = 10, height = 8, dpi = 500)
  
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
  enrichment_df <- data.frame(Phylostrata = phylostrata_levels, p_adj = adjusted_pvals)
  enrichment_plot <- ggplot(enrichment_df, aes(x = as.factor(Phylostrata), y = -log10(p_adj), fill = as.factor(Phylostrata))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = Inf, y = -log10(0.05), label = "p = 0.05", hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
    scale_fill_viridis_d() + theme_bw() +
    labs(title = paste0("Phylostrata Enrichment (Fisher Test, BH-adjusted) in ", cancer_type_name," (OncoDB)"),
         x = "Phylostrata", y = "-log10(adjusted p-value)") + theme(legend.position = "none")
  ggsave(paste0("plots/OncoDB/", project_id, "/OncoDB_", project_id, "_enrichment.png"), plot = enrichment_plot, width = 8, height = 6, dpi = 300)
  
  # Bootstrapping
  set.seed(42)
  boot_df <- data.frame()
  for (i in 1:100) {
    sampled <- data_clean %>% arrange(desc(abs(log2FoldChange))) %>% slice_sample(n = top_n)
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
  ggsave(paste0("plots/OncoDB/", project_id, "/OncoDB_", project_id, "_bootstrap.png"), plot = boot_plot, width = 8, height = 6, dpi = 300)
} # End of loop

### === Combined Fisher Test Across All Cancers === ###
# Combined Fisher Test Analysis Across All OncoDB Cancers
data_all <- list()

# Collect top DE genes across all projects
for (project_path in files) {
  data <- read.table(project_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange")
  data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))
  data <- data %>% filter(gene_name %in% phylostrata$gene_name) %>% left_join(phylostrata, by = "gene_name")
  data_clean <- data %>% filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))
  data_all[[basename(project_path)]] <- data_clean
}

# Combine all DE results
data_combined <- bind_rows(data_all)

# Define top N
top_n_combined <- 160

# Get top N DE genes across all cancers
combined_top_genes <- data_combined %>%
  filter(padj < 0.001) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n_combined)

# Fisher test
top_counts_all <- table(factor(combined_top_genes$Phylostrata, levels = 1:16))
background_counts_all <- table(factor(data_combined$Phylostrata, levels = 1:16))
contingency_all <- rbind(Top = top_counts_all, Background = background_counts_all - top_counts_all)
fisher_all <- fisher.test(contingency_all, simulate.p.value = TRUE, B = 10000)
print(fisher_all)
capture.output(fisher_all, file = "stats/OncoDB/OncoDB_combined_fisher.txt")

# Top N Bar Plot
combined_barplot <- ggplot(combined_top_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = as.factor(Phylostrata))) +
  geom_bar(stat = "identity") + coord_flip() +
  theme_bw() + scale_fill_viridis_d(name = "Phylostrata") +
  labs(title = paste0("Top ", top_n_combined, " DE Genes Across All Cancers (FDR < 0.001) (OncoDB)"),
       x = "Gene", y = "Log2 Fold Change") +
  theme(axis.text.y = element_text(size = 4), plot.title = element_text(size = 12))

ggsave("plots/OncoDB/OncoDB_combined_top_n_barplot.png", plot = combined_barplot, width = 10, height = 8, dpi = 500)

# Enrichment Plot
phylostrata_levels <- sort(unique(data_combined$Phylostrata))
enrichment_pvals <- sapply(phylostrata_levels, function(ps) {
  in_top <- sum(combined_top_genes$Phylostrata == ps)
  in_bg <- sum(data_combined$Phylostrata == ps) - in_top
  not_in_ps_top <- top_n_combined - in_top
  not_in_ps_bg <- nrow(data_combined) - top_n_combined - in_bg
  matrix <- matrix(c(in_top, in_bg, not_in_ps_top, not_in_ps_bg), nrow = 2)
  fisher.test(matrix)$p.value
})
adjusted_pvals <- p.adjust(enrichment_pvals, method = "BH")
enrichment_df <- data.frame(Phylostrata = phylostrata_levels, p_adj = adjusted_pvals)

enrichment_plot <- ggplot(enrichment_df, aes(x = as.factor(Phylostrata), y = -log10(p_adj), fill = as.factor(Phylostrata))) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = Inf, y = -log10(0.05), label = "p = 0.05", hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
  scale_fill_viridis_d() + theme_bw() +
  labs(title = "Phylostrata Enrichment Across All Cancers (BH-adjusted) (OncoDB)",
       x = "Phylostrata", y = "-log10(adjusted p-value)") + theme(legend.position = "none")

ggsave("plots/OncoDB/OncoDB_combined_enrichment.png", plot = enrichment_plot, width = 8, height = 6, dpi = 300)

# Bootstrapped Analysis
set.seed(42)
boot_df <- data.frame()
for (i in 1:100) {
  sampled <- data_combined %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_sample(n = top_n_combined)
  for (ps in phylostrata_levels) {
    in_top <- sum(sampled$Phylostrata == ps)
    in_bg <- sum(data_combined$Phylostrata == ps) - in_top
    not_in_ps_top <- top_n_combined - in_top
    not_in_ps_bg <- nrow(data_combined) - top_n_combined - in_bg
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
  labs(title = "Bootstrapped Phylostrata Enrichment Across All Cancers (OncoDB)",
       x = "Phylostrata", y = "Median Odds Ratio")

ggsave("plots/OncoDB/OncoDB_combined_bootstrap_enrichment.png", plot = boot_plot, width = 8, height = 6, dpi = 300)

### === Meta-analysis of Fisher P-values === ###

# List all Fisher test files
fisher_files <- list.files("stats/OncoDB/fisher", pattern = "_fisher_test\\.txt$", full.names = TRUE)

# Initialise lists to store results & wipe previous inputs
results <- list()
meta_df <- data.frame()

for (f in fisher_files) {
  # Extract project ID from filename
  project_id <- str_extract(basename(f), "(?<=OncoDB_)[A-Z]+(?=_fisher_test)")
  
  # Get cancer name from project_info
  cancer_name <- project_info %>%
    filter(project_id == !!project_id) %>%
    pull(cancer_type)

  lines <- readLines(f)
  
  # Find the p-value line
  p_val_line <- grep("p-value\\s*=\\s*", lines, value = TRUE)
  
  if (length(p_val_line) > 0) {
    # Extract the number using regex
    p_val_match <- str_match(p_val_line[1], "p-value\\s*=\\s*([0-9.eE-]+)")[,2]
    p_val <- as.numeric(p_val_match)
    
    if (!is.na(p_val)) {
      results[[length(results) + 1]] <- data.frame(
        project_id = project_id,
        cancer_type = cancer_name,
        pval = p_val,
        stringsAsFactors = FALSE
      )
    } else {
      message("Failed to convert p-value for: ", project_id)
    }
  } else {
    message("No p-value line found in: ", project_id)
  }
}

# Combine all results into one dataframe
meta_df <- bind_rows(results) %>%
  mutate(log_pval = -log10(pval)) %>%
  arrange(pval)  # sort by p-value

# Check if we have any results
if (nrow(meta_df) == 0) {
  stop("No valid results found - check your input files")
}

# Plot meta-analysis
meta_plot <- ggplot(meta_df, aes(x = reorder(cancer_type, log_pval), y = log_pval, fill = log_pval)) +
  geom_bar(stat = "identity") +
  coord_flip(ylim = c(0, 4)) +
  scale_fill_viridis_c() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = -log10(0.05), label = "p = 0.05", 
           hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Fisher Test Significance Across Cancer Types (OncoDB)",
    x = "Cancer Type",
    y = expression(-log[10](p-value))
  )
show(meta_plot)

# Save plot
ggsave("plots/OncoDB/OncoDB_fisher_meta_pvals_barplot.png", plot = meta_plot, width = 10, height = 6, dpi = 300)

# Save summary table
write.csv(meta_df, "stats/OncoDB/OncoDB_meta_fisher_pvals.csv", row.names = FALSE)

# Combine using Fisher's method (meta-analysis)
pvals_meta_clean <- na.omit(meta_df$pval)
if (length(pvals_meta_clean) > 0) {
  meta_result <- sumlog(pvals_meta_clean)
  capture.output(print(meta_result), file = "stats/OncoDB/OncoDB_meta_analysis_fisher.txt")
  print(meta_result)
} else {
  warning("No valid p-values found for meta-analysis")
}
