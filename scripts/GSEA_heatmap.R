setwd("~/Desktop/gene_age_project/")

library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)
library(gtable)
library(ComplexHeatmap)
library(circlize)
library(viridis)

### TCGA ###

# Find GSEA results files
gsea_dir <- "stats/TCGA/GSEA/"
gsea_files <- list.files(
  path = gsea_dir,
  pattern = "^[A-Z]+_TCGA_gsea_results\\.csv$",  # Exclude combined summary file
  full.names = TRUE
)

# Combine NES scores from all files
gsea_list <- lapply(gsea_files, function(file) {
  project_id <- str_extract(basename(file), "^[A-Z]+(?=_TCGA_gsea_results\\.csv)")
  df <- read_csv(file, show_col_types = FALSE) %>%
    select(pathway, NES) %>%
    mutate(project_id = project_id)
  return(df)
})

# Merge into long format then pivot to matrix
gsea_combined <- bind_rows(gsea_list)

nes_matrix <- gsea_combined %>%
  pivot_wider(names_from = pathway, values_from = NES) %>%
  column_to_rownames("project_id") %>%
  as.matrix()

# Remove any rows/columns with NA
nes_matrix_clean <- nes_matrix
nes_matrix_clean[is.na(nes_matrix_clean)] <- 0

# Change phylostrata order to accending numerically
col_order <- as.character(sort(as.numeric(colnames(nes_matrix_clean))))
nes_matrix_clean <- nes_matrix_clean[, col_order]

# Plot heatmap w/ ComplexHeatmap
ht <- Heatmap(
  nes_matrix_clean,
  name = "NES",  # Legend title
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "Phylostrata",  # x-axis title
  row_title = "Cancer Types",    # y-axis title
  column_names_side = "top",     # Move phylostrata labels to top
  row_names_side = "left",       # Move cancer type labels to left
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_rot = 0, # Fix phylostrata text rotation
  heatmap_legend_param = list(
    title = "Normalised Enrichment Scores (NES)",
    title_position = "topcenter",     # Title above the vertical legend
    legend_direction = "horizontal",    # Horizontal colour bar
    title_gp = gpar(fontsize = 12, fontface = "bold"),   # Legend title size
    labels_gp = gpar(fontsize = 10),                     # Legend text size
    legend_width = unit(8, "cm"),                        # Colour bar width
    legend_height = unit(0.7, "cm")                      # Colour bar height
  )
)

# Save plot as png
png("plots/TCGA/GSEA/TCGA_GSEA_heatmap.png", width = 12, height = 6, units = "in", res = 300)

# Draw with a top-level plot title
draw(ht,
     heatmap_legend_side = "bottom",
     column_title = "Phylostrata Enrichment (NES) Across TCGA Cancers",
     column_title_gp = gpar(fontsize = 16)
)

dev.off()

### OncoDB ###

# Find GSEA results files
gsea_dir <- "stats/OncoDB//GSEA/"
gsea_files <- list.files(
  path = gsea_dir,
  pattern = "^[A-Z]+_OncoDB_gsea_results\\.csv$",  # Exclude combined summary file
  full.names = TRUE
)

# Combine NES scores from all files
gsea_list <- lapply(gsea_files, function(file) {
  project_id <- str_extract(basename(file), "^[A-Z]+(?=_OncoDB_gsea_results\\.csv)")
  df <- read_csv(file, show_col_types = FALSE) %>%
    select(pathway, NES) %>%
    mutate(project_id = project_id)
  return(df)
})

# Merge into long format then pivot to matrix
gsea_combined <- bind_rows(gsea_list)

nes_matrix <- gsea_combined %>%
  pivot_wider(names_from = pathway, values_from = NES) %>%
  column_to_rownames("project_id") %>%
  as.matrix()

# Remove any rows/columns with NA
nes_matrix_clean <- nes_matrix
nes_matrix_clean[is.na(nes_matrix_clean)] <- 0

# Change phylostrata order to accending numerically
col_order <- as.character(sort(as.numeric(colnames(nes_matrix_clean))))
nes_matrix_clean <- nes_matrix_clean[, col_order]

# Plot heatmap w/ ComplexHeatmap
ht <- Heatmap(
  nes_matrix_clean,
  name = "NES",  # Legend title
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "Phylostrata",  # x-axis title
  row_title = "Cancer Types",    # y-axis title
  column_names_side = "top",     # Move phylostrata labels to top
  row_names_side = "left",       # Move cancer type labels to left
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_names_rot = 0, # Fix phylostrata text rotation
  heatmap_legend_param = list(
    title = "Normalised Enrichment Scores (NES)",
    title_position = "topcenter",     # Title above the vertical legend
    legend_direction = "horizontal",    # Horizontal colour bar
    title_gp = gpar(fontsize = 12, fontface = "bold"),   # Legend title size
    labels_gp = gpar(fontsize = 10),                     # Legend text size
    legend_width = unit(8, "cm"),                        # Colour bar width
    legend_height = unit(0.7, "cm")                      # Colour bar height
  )
)

# Save plot as png
png("plots/OncoDB/GSEA/OncoDB_GSEA_heatmap.png", width = 12, height = 6, units = "in", res = 300)

# Draw with a top-level plot title
draw(ht,
     heatmap_legend_side = "bottom",
     column_title = "Phylostrata Enrichment (NES) Across OncoDB Cancers",
     column_title_gp = gpar(fontsize = 16)
)

dev.off()
