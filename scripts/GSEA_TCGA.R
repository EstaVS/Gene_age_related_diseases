setwd("~/Desktop/gene_age_project/")
library(readr)
library(dplyr)
library(ggplot2)
library(fgsea)
library(viridis)
library(ggrepel)
library(stringr)

# Load metadata
project_info <- read_csv("metadata/project_info.csv")
files <- list.files("TCGA_DE_results", pattern = "^TCGA-[A-Z]+_DE\\.csv$", full.names = TRUE)

# Create output directories
dir.create("stats/TCGA/GSEA", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/TCGA/GSEA", recursive = TRUE, showWarnings = FALSE)

## Main Analysis Loop
for (project_path in files) {
  project_id <- str_extract(basename(project_path), "(?<=TCGA-)[A-Z0-9]+(?=_DE\\.csv)")
  cancer_type_name <- project_info %>%
    filter(project_id == !!project_id) %>%
    pull(cancer_type)
  
  cat("\nProcessing:", project_id, "-", cancer_type_name, "\n")
  dir.create(paste0("plots/TCGA/GSEA/", project_id), recursive = TRUE, showWarnings = FALSE)
  
  # Read and prepare data
  data <- read.csv(project_path, stringsAsFactors = FALSE) %>%
    select(gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
           gene_id, Entrez, Phylostrata, cancer_median, normal_median) %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))
  
  # Create phylostrata gene sets from the file directly
  phylostrata_sets <- split(data$gene_name, data$Phylostrata)
  
  ranked_genes <- data %>%
    mutate(
      ranking_metric = -log10(padj + 1e-300) * sign(log2FoldChange) +
        sign(log2FoldChange) * (abs(log2FoldChange) / max(abs(log2FoldChange), na.rm = TRUE) / 1e6) +
        (cancer_median / max(cancer_median, na.rm = TRUE) / 1e9) +
        (normal_median / max(normal_median, na.rm = TRUE) / 1e12) +
        runif(n(), 0, 1e-12)
    ) %>%
    filter(is.finite(ranking_metric)) %>%
    arrange(desc(ranking_metric)) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    {setNames(.$ranking_metric, .$gene_name)}
  
  # Report tie percentage
  tie_percentage <- mean(duplicated(ranked_genes)) * 100
  cat("Tie percentage after breaking:", round(tie_percentage, 2), "%\n")
  
  # Run GSEA
  set.seed(42)
  gsea_results <- fgseaMultilevel(
    pathways = phylostrata_sets,
    stats = ranked_genes,
    minSize = 15,
    maxSize = 1000,
    sampleSize = 100,
    eps = 0,
    nPermSimple = 10000
  )
  
  saveRDS(gsea_results, paste0("stats/TCGA/GSEA/", project_id, "_TCGA_gsea_results.rds"))
  
  ### Visualisation 1: GSEA Enrichment Plot (Top Pathway)
  top_pathway <- gsea_results %>%
    arrange(padj) %>%
    slice(1) %>%
    pull(pathway)
  
  if (length(top_pathway) > 0) {
    enrichment_plot <- plotEnrichment(
      phylostrata_sets[[top_pathway]],
      ranked_genes
    ) + 
      labs(title = paste("GSEA Enrichment for phylostrata", top_pathway, "in", project_id, "(TCGA)"),
           subtitle = paste("NES =", round(gsea_results$NES[1], 2), 
                            "FDR =", format.pval(gsea_results$padj[1], digits = 2))) +
      theme(plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.key.size = unit(1, "cm"))
    
    ggsave(paste0("plots/TCGA/GSEA/", project_id, "/", project_id, "_TCGA_enrichment_plot.png"),
           plot = enrichment_plot, width = 8, height = 6, dpi = 300)
  }
  
  ### Visualisation 2: NES barplot
  gsea_plot_data <- gsea_results %>%
    filter(padj < 0.2) %>%
    mutate(
      PhyloNum = as.numeric(gsub("[^0-9]", "", pathway)), 
      pathway = factor(pathway, levels = pathway[order(PhyloNum)]), 
      Significance = case_when(
        padj < 0.01 ~ "FDR < 0.01",
        padj < 0.05 ~ "FDR < 0.05",
        TRUE ~ "FDR < 0.2"
      )
    )
  
  if (nrow(gsea_plot_data) > 0) {
    nes_plot <- ggplot(gsea_plot_data, aes(x = pathway, y = NES, fill = Significance)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_viridis_d(name = "Significance", direction = -1) +
      coord_flip() +
      theme_bw() +
      labs(title = paste("Phylostrata Enrichment in", project_id, "(TCGA)"),
           x = "Phylostrata",
           y = "Normalised Enrichment Score (NES)") +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.key.size = unit(1, "cm"))
    
    ggsave(paste0("plots/TCGA/GSEA/", project_id, "/", project_id, "_TCGA_nes_barplot.png"),
           plot = nes_plot, width = 8, height = 6, dpi = 300)
  }
  
  ### Visualisation 3: Leading Edge Genes
  if (length(top_pathway) > 0) {
    
    # Automatically get correct index for top_pathway
    top_index <- which(gsea_results$pathway == top_pathway)
    
    if (!is.null(gsea_results$leadingEdge[[top_index]])) {
      leading_genes <- gsea_results$leadingEdge[[top_index]]
    
      leading_data <- data %>%
        filter(gene_name %in% leading_genes) %>%
        arrange(desc(log2FoldChange)) %>%
        slice_head(n = 30) %>% # Capped at top 30 genes for visual reasons
        select(gene_name, log2FoldChange, Phylostrata) %>%
        mutate(Phylostrata = factor(Phylostrata))
      
      leading_edge_plot <- ggplot(leading_data, 
                                  aes(x = reorder(gene_name, log2FoldChange), 
                                      y = log2FoldChange, 
                                      fill = Phylostrata)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_viridis_d(name = "Phylostrata") +
        theme_bw() +
        labs(title = paste0("Leading Edge Genes from Phylostrata ", top_pathway),
             subtitle = paste("in", project_id, "(TCGA)"),
             x = "Gene",
             y = "Log2 Fold Change") +
        theme(axis.text.y = element_text(size = 8),
              plot.title = element_text(size = 16, face = "bold"),
              plot.subtitle = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              axis.text.x = element_text(size = 12),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.key.size = unit(1, "cm"))
      
      ggsave(paste0("plots/TCGA/GSEA/", project_id, "/", project_id, "_TCGA_leading_edge_genes.png"),
             plot = leading_edge_plot, width = 10, height = 8, dpi = 300)
      }
    }
  
  ### Visualisation 4: Volcano Plot with Top Genes
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
    labs(title = paste0("Volcano Plot: ", project_id, " (TCGA)"),
         x = "Log2 Fold Change",
         y = "-log10(Adjusted p-value)") +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"))
  
  ggsave(paste0("plots/TCGA/GSEA/", project_id, "/", project_id, "_TCGA_volcano_plot.png"),
         plot = volcano_plot, width = 8, height = 8, dpi = 300)
  
  ### Visualisation 5: Ranking Distribution
  rank_df <- data.frame(
    gene = names(ranked_genes),
    rank = 1:length(ranked_genes),
    value = ranked_genes
  ) %>% left_join(data, by = c("gene" = "gene_name"))
  
  rank_plot <- ggplot(rank_df, aes(x = rank, y = value, color = as.factor(Phylostrata))) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_viridis_d(name = "Phylostrata") +
    theme_bw() +
    labs(title = paste0("Gene Ranking Distribution: ", project_id, " (TCGA)"),
         x = "Rank",
         y = "Ranking Metric") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.5, "cm"))
  
  ggsave(paste0("plots/TCGA/GSEA/", project_id, "/", project_id, "_TCGA_ranking_distribution.png"),
         plot = rank_plot, width = 10, height = 6, dpi = 300)
}


### Combined Analysis Across All Cancers (with visualisation) ###

if (length(files) > 1) {
  
  # Combine all data
  combined_data <- lapply(files, function(f) {
    data <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
    colnames(data) <- c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
           "gene_id", "Entrez", "Phylostrata", "cancer_median", "normal_median")
    data$padj <- as.numeric(gsub("q<=|Q<=|Q=|q=", "", data$padj))
    data %>% filter(!is.na(padj), !is.na(log2FoldChange), !is.na(Phylostrata))
  }) %>% bind_rows()
  
  # Create combined ranked list
  combined_ranked <- combined_data %>%
    mutate(
      ranking_metric = -log10(padj + 1e-300) * sign(log2FoldChange) +
        sign(log2FoldChange) * (abs(log2FoldChange) / max(abs(log2FoldChange), na.rm = TRUE) / 1e3) +
        (cancer_median/max(cancer_median, na.rm = TRUE)/1e6) +
        (normal_median/max(normal_median, na.rm = TRUE)/1e9)
    ) %>%
    filter(is.finite(ranking_metric)) %>%
    arrange(desc(ranking_metric)) %>%
    {setNames(.$ranking_metric, .$gene_name)}
  
  # Run combined GSEA
  set.seed(42)
  combined_gsea <- fgseaMultilevel(
    pathways = phylostrata_sets,
    stats = combined_ranked,
    minSize = 15,
    maxSize = 1000,
    sampleSize = 100
  )
  
  write_csv(combined_gsea, "stats/TCGA/GSEA/TCGA_combined_gsea_results.csv")
  
  # Combined NES Plot with numeric ordering
  combined_nes_plot <- combined_gsea %>%
    filter(padj < 0.2) %>%
    mutate(
      PhyloNum = as.numeric(gsub("[^0-9]", "", pathway)),
      pathway = factor(pathway, levels = pathway[order(PhyloNum)])
    ) %>%
    ggplot(aes(x = pathway, y = NES, fill = -log10(padj))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_viridis_c(name = "-log10(FDR)") +
    coord_flip() +
    theme_bw() +
    labs(title = "Combined Phylostrata Enrichment Across All Cancers (TCGA)",
         x = "Phylostrata",
         y = "Normalised Enrichment Score (NES)") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1, "cm"))
  show(combined_nes_plot)
  
  ggsave("plots/TCGA/GSEA/TCGA_combined_nes_plot.png",
         plot = combined_nes_plot, width = 10, height = 8, dpi = 300)
}
