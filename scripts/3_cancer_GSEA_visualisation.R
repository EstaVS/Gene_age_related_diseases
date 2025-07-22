# Set working directory and load required packages
setwd("~/Desktop/gene_age_project/")
library(readr)
library(dplyr)
library(ggplot2)
library(fgsea)
library(viridis)
library(ggrepel)
library(stringr)
library(patchwork)  # for plot composition

# Define cancer types and data sources
cancers <- c("BRCA", "KIRC", "LIHC")
sources <- c("TCGA", "OncoDB")

# Initialise plot lists
enrichment_list <- list()
nes_list <- list()
leading_list <- list()
volcano_list <- list()
rank_list <- list()

for (source in sources) {
  for (cancer in cancers) {
    
    gsea_path <- paste0("stats/", source, "/GSEA/", cancer, "_", source, "_gsea_results.rds")
    expr_path <- if (source == "TCGA") {
      paste0("TCGA_DE_results/TCGA-", cancer, "_DE.csv")
    } else {
      paste0("OncoDB_data_annotated/expression_diff_", cancer, ".csv")
    }
    
    gsea_results <- readRDS(gsea_path)
    data <- read_csv(expr_path, show_col_types = FALSE) %>%
      filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata)) %>%
      select(any_of(c(
        "gene_name", "log2FoldChange", "padj", "Phylostrata",
        "baseMean", "lfcSE", "stat", "pvalue",
        "gene_id", "Entrez", "cancer_median", "normal_median"
      )))
    
    data <- data %>%
      mutate(ranking_metric = -log10(padj + 1e-300) * sign(log2FoldChange) +
               sign(log2FoldChange) * (abs(log2FoldChange) / max(abs(log2FoldChange), na.rm = TRUE) / 1e6) +
               (if ("cancer_median" %in% names(.)) cancer_median / max(cancer_median, na.rm = TRUE) / 1e9 else 0) +
               (if ("normal_median" %in% names(.)) normal_median / max(normal_median, na.rm = TRUE) / 1e12 else 0) +
               runif(n(), 0, 1e-12))
    
    ranked_genes <- data %>%
      filter(is.finite(ranking_metric)) %>%
      arrange(desc(ranking_metric)) %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      {setNames(.$ranking_metric, .$gene_name)}
    
    phylostrata_sets <- split(data$gene_name, data$Phylostrata)
    
    top_pathway <- gsea_results %>%
      arrange(padj) %>%
      filter(!is.na(pathway)) %>%
      slice(1) %>%
      pull(pathway)
    
    if (top_pathway %in% names(phylostrata_sets)) {
      enrichment_list[[paste(cancer, source, sep = "_")]] <- 
        plotEnrichment(phylostrata_sets[[top_pathway]], ranked_genes) +
        labs(title = paste0(cancer, " (", source,")", " for PS", top_pathway),
             subtitle = paste0("NES = ", round(gsea_results$NES[1], 2), ", FDR = ", format.pval(gsea_results$padj[1], digits = 2))) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none")
    }
    
    nes_data <- gsea_results %>%
      filter(padj < 0.2) %>%
      mutate(PhyloNum = as.numeric(gsub("[^0-9]", "", pathway)),
             pathway = factor(pathway, levels = pathway[order(PhyloNum)]),
             Significance = case_when(
               padj < 0.01 ~ "FDR < 0.01",
               padj < 0.05 ~ "FDR < 0.05",
               TRUE ~ "FDR < 0.2"
             ))
    
    nes_data <- nes_data %>%
      mutate(Significance = factor(Significance, levels = c("FDR < 0.01", "FDR < 0.05", "FDR < 0.2")))
    
    if (nrow(nes_data) > 0) {
      nes_list[[paste(cancer, source, sep = "_")]] <- ggplot(nes_data, aes(x = pathway, y = NES, fill = Significance)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(
          values = c(
            "FDR < 0.01" = "#FDE725", 
            "FDR < 0.05" = "#21918C", 
            "FDR < 0.2"  = "#440154" 
          )
        ) +
        coord_flip() +
        scale_y_continuous(limits = c(-4, 4)) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none") +
        labs(title = paste0(cancer, " (", source, ")"), x = "Phylostrata", y = "NES")
    }
    
    top_index <- which(gsea_results$pathway == top_pathway)
    if (!is.null(gsea_results$leadingEdge[[top_index]])) {
      leading_data <- data %>%
        filter(gene_name %in% gsea_results$leadingEdge[[top_index]]) %>%
        arrange(desc(log2FoldChange)) %>%
        slice_head(n = 30)
      
      leading_list[[paste(cancer, source, sep = "_")]] <- ggplot(leading_data, aes(x = reorder(gene_name, log2FoldChange),
                                                                                   y = log2FoldChange, fill = factor(Phylostrata))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_viridis_d(name = "Phylostrata") +
        labs(title = paste0(cancer, " (", source, ")", " for PS", top_pathway),
             x = "Gene Name") +
        theme_bw(base_size = 10) +
        theme(legend.position = "none") 
    }
    
    volcano_data <- data %>%
      mutate(Significance = case_when(
        padj < 0.001 & abs(log2FoldChange) > 2 ~ "FDR < 0.001 & |FC| > 2",
        padj < 0.05 ~ "FDR < 0.05",
        TRUE ~ "Not significant"
      ))
    top_genes <- volcano_data %>% arrange(padj) %>% slice_head(n = 20)
    volcano_list[[paste(cancer, source, sep = "_")]] <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("red", "orange", "gray60")) +
      geom_text_repel(data = top_genes, aes(label = gene_name), size = 2.5) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none") +
      labs(title = paste0(cancer, " (", source, ")"),
           x = "Log2 Fold Change", y = "-log10(FDR)")
    
    rank_df <- data.frame(
      gene = names(ranked_genes),
      rank = seq_along(ranked_genes),
      value = ranked_genes
    ) %>% left_join(data, by = c("gene" = "gene_name"))
    rank_list[[paste(cancer, source, sep = "_")]] <- ggplot(rank_df, aes(x = rank, y = value, color = as.factor(Phylostrata))) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_viridis_d(name = "Phylostrata") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none") +
      labs(title = paste0(cancer, " (", source, ")"),
           x = "Rank", y = "Ranking Metric")
  }
}

# Assemble composite visualisations for each plot type
enrichment_combined <- wrap_plots(enrichment_list, ncol = 3) + 
    plot_annotation(title = "GSEA Enrichment Plots for Top Phylostratum",
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

nes_combined <- wrap_plots(nes_list, ncol = 3) +
    plot_annotation(title = "Normalised Enrichment Scores (NES) Barplots",
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

leading_combined <- wrap_plots(leading_list, ncol = 3) +
  plot_annotation(title = "Leading Edge Genes Barplots for Top Phylostratum",
                  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

volcano_combined <- wrap_plots(volcano_list, ncol = 3) +
    plot_annotation(title = "Volcano Plots",
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5 )))

rank_combined <- wrap_plots(rank_list, ncol = 3) + 
    plot_annotation(title = "Gene Ranking Distributions Plots",
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5 )))

# Display all combined plots
print(enrichment_combined)
print(nes_combined)
print(leading_combined)
print(volcano_combined)
print(rank_combined)

# Save plots
ggsave(paste0("plots/Final_report/enrichment_plot.png"),
       plot = enrichment_combined, width = 10, height = 7, dpi = 300)
ggsave(paste0("plots/Final_report/NES_plot.png"),
       plot = nes_combined, width = 10, height = 6, dpi = 300)
ggsave(paste0("plots/Final_report/leading_edge_plot.png"),
       plot = leading_combined, width = 10, height = 8, dpi = 300)
ggsave(paste0("plots/Final_report/volcano_plot.png"),
       plot = volcano_combined, width = 12, height = 8, dpi = 300)
ggsave(paste0("plots/Final_report/ranking_dist_plot.png"),
       plot = rank_combined, width = 10, height = 7, dpi = 300)