setwd("~/Desktop/gene_age_project/")
library(readr)
library(dplyr)
library(ggplot2)
library(fgsea)
library(viridis)
library(ggrepel)

# Load data
project_info <- read_csv("metadata/project_info.csv")

# Project ID and Cancer Name
project_id <- str_extract(basename("OncoDB_data_annotated/expression_diff_ACC.csv"), "(?<=expression_diff_)[^\\.]+")
cancer_type_name <- project_info %>%
  filter(project_id == !!project_id) %>%
  pull(cancer_type)

data <- read.csv("OncoDB_data_annotated/expression_diff_ACC.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(data) <- c("gene_name", "padj", "cancer_median", "normal_median", "log2FoldChange", "Entrez", "Phylostrata")
data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))
  
# Create phylostrata gene sets
phylostrata_sets <- split(data$gene_name, data$Phylostrata)

data <- data %>% 
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))
  
# Enhanced ranking metric with aggressive tie-breaking
ranked_genes <- data %>%
  mutate(
    ranking_metric = -log10(padj) * sign(log2FoldChange) +
        sign(log2FoldChange) * (abs(log2FoldChange)/max(abs(log2FoldChange))/1e6) +
        (cancer_median/max(cancer_median, na.rm = TRUE)/1e9) +
        (normal_median/max(normal_median, na.rm = TRUE)/1e12) +
        runif(n(), 0, 1e-12)  # Add tiny random noise to break any remaining ties
    ) %>%
    arrange(desc(ranking_metric)) %>%
    {setNames(.$ranking_metric, .$gene_name)}
  
# Report tie percentage (should be 0% after this treatment)
tie_percentage <- mean(duplicated(ranked_genes)) * 100
  cat("Tie percentage after breaking:", round(tie_percentage, 2), "%\n")
  
# Run GSEA with strict parameters
set.seed(42)
gsea_results <- fgseaMultilevel(
    pathways = phylostrata_sets,
    stats = ranked_genes,
    minSize = 15,
    maxSize = 1000,
    sampleSize = 100,
    eps = 0,  # Most precise p-value calculation
    nPermSimple = 10000  # More permutations for problematic pathways
  )

### VISUALIZATION 1: GSEA Enrichment Plot (Top Pathway)
top_pathway <- gsea_results %>% 
    arrange(padj) %>% 
    slice(1) %>% 
    pull(pathway)
show(top_pathway)

if (length(top_pathway) > 0) {
    enrichment_plot <- plotEnrichment(
      phylostrata_sets[[top_pathway]],
      ranked_genes
    ) + 
      labs(title = paste("GSEA Enrichment for Phylostrata", top_pathway, "in", cancer_type_name, "(OncoDB)"),
           subtitle = paste("NES =", round(gsea_results$NES[1], 2), 
                            "FDR =", format.pval(gsea_results$padj[1], digits = 2)),
           x = "Rank",
           y = "Enrichment Score") +
      theme(plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10))
}
show(enrichment_plot)
  
### VISUALIZATION 2: NES Barplot
gsea_plot_data <- gsea_results %>%
    mutate(pathway = factor(pathway, levels = pathway[order(NES)])) %>%
    filter(padj < 0.2) %>%
    mutate(Significance = case_when(
      padj < 0.01 ~ "FDR < 0.01",
      padj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "FDR < 0.2"
    ))
  
  if (nrow(gsea_plot_data) > 0) {
    nes_plot <- ggplot(gsea_plot_data, aes(x = pathway, y = NES, fill = Significance)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_viridis_d(name = "Significance", direction = -1) +
      coord_flip() +
      theme_bw() +
      labs(title = paste("Phylostrata Enrichment in", cancer_type_name, "(OncoDB)"),
           x = "Phylostrata",
           y = "Normalised Enrichment Score (NES)") +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 12, face = "bold"))
  }
show(nes_plot)

### VISUALIZATION 3: Leading Edge Genes
if (length(top_pathway) > 0) {
  
  # Automatically get correct index for top_pathway
  top_index <- which(gsea_results$pathway == top_pathway)
  
  if (!is.null(gsea_results$leadingEdge[[top_index]])) {
    leading_genes <- gsea_results$leadingEdge[[top_index]]
    
    leading_data <- data %>%
      filter(gene_name %in% leading_genes) %>%
      arrange(desc(log2FoldChange)) %>%
      select(gene_name, log2FoldChange, Phylostrata) %>%
      mutate(Phylostrata = factor(Phylostrata))  # no fixed levels
    
    leading_edge_plot <- ggplot(leading_data, 
                                aes(x = reorder(gene_name, log2FoldChange), 
                                    y = log2FoldChange, 
                                    fill = Phylostrata)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_viridis_d(name = "Phylostrata") +  # only shows what's present
      theme_bw() +
      labs(title = paste0("Leading Edge Genes from Phylostrata ", top_pathway),
           subtitle = paste("in", cancer_type_name, "(OncoDB)"),
           x = "Gene",
           y = "Log2 Fold Change") +
      theme(axis.text.y = element_text(size = 8),
            plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 10))
    
      show(leading_edge_plot)
    }
  }


### VISUALIZATION 4: Volcano Plot with Top Genes
volcano_data <- data %>%
    mutate(Significance = case_when(
      padj < 0.001 & abs(log2FoldChange) > 2 ~ "FDR < 0.001 & |FC| > 2",
      padj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not significant"
    ))
  
top_genes <- data %>%
    arrange(padj) %>%
    head(20) %>%
    mutate(Significance = case_when(
      padj < 0.001 & abs(log2FoldChange) > 2 ~ "FDR < 0.001 & |FC| > 2",
      padj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not significant"
    ))
  
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("red", "orange", "gray60")) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene_name),
      color = "black",
      size = 3,
      box.padding = 0.5,
      max.overlaps = 50
    ) +
    theme_bw() +
    labs(title = paste0("Volcano Plot: ", cancer_type_name, " (OncoDB)"),
         x = "Log2 Fold Change",
         y = "-log10(Adjusted p-value)") +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, face = "bold"))
show(volcano_plot)
  
### VISUALIZATION 5: Ranking Distribution
rank_df <- data.frame(
    gene = names(ranked_genes),
    rank = 1:length(ranked_genes),
    value = ranked_genes
  )
  
rank_plot <- ggplot(rank_df, aes(x = rank, y = value, color = as.factor(data$Phylostrata))) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_viridis_d(name = "Phylostrata") +
    theme_bw() +
    labs(title = paste0("Gene Ranking Distribution: ", cancer_type_name, " (OncoDB)"),
         x = "Rank",
         y = "Ranking Metric") +
    theme(plot.title = element_text(size = 12, face = "bold"))
show(rank_plot)
