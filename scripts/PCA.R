setwd("~/Desktop/gene_age_project/")

library(readr)
library(dplyr)
library(stringr)
library(FactoMineR)
library(factoextra)

### TCGA ###

# TCGA GSEA results
tcga_gsea_files <- list.files("stats/TCGA/GSEA", pattern = "_gsea_results.csv$", full.names = TRUE)
tcga_gsea_files <- tcga_gsea_files[!grepl("TCGA_combined_gsea_results.csv", tcga_gsea_files)] # Exclude combined GSEA results file

# store NES values in a list
nes_matrix <- list()

for (file in tcga_gsea_files) {
  project_id <- str_extract(basename(file), "^[^_]+")
  
  gsea_res <- read_csv(file) %>%
    select(pathway, NES) %>%
    filter(!is.na(NES))
  
  # use pathways as rownames
  rownames(gsea_res) <- gsea_res$pathway
  gsea_res$pathway <- NULL
  
  nes_matrix[[project_id]] <- gsea_res$NES
}

# Create NES matrix
nes_mat <- bind_rows(
  lapply(nes_matrix, function(x) as.data.frame(t(x))),
  .id = "project_id"
) %>%
  column_to_rownames("project_id") %>%
  as.matrix()

# Removing columns (i.e. Phylostrata) with more than 20% NAs
nes_mat_filtered <- nes_mat[, colMeans(is.na(nes_mat)) < 0.2]

# Imputing NAs
nes_mat_impute <- nes_mat # Copy of nes_mat to preserve original
nes_mat_impute[is.na(nes_mat_impute)] <- mean(nes_mat, na.rm = TRUE) # Replaced NAs with column mean

# Run PCA
pca_res <- prcomp(nes_mat_impute, scale. = TRUE)

# Plot
tcga_pca <- fviz_pca_biplot(
  pca_res,
  label = "var",               
  habillage = rownames(nes_mat),  # color by project_id
  addEllipses = FALSE,
  repel = TRUE,
  geom = "point",
  pointshape = 20   
) +
  labs(
    title = "PCA of GSEA Normalised Enrichment Score by Cancer (TCGA)",
    subtitle = "+ Additional Phylostrata contributions",
    x = paste0("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2), 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "%)")
  ) +
  theme_minimal()
show(tcga_pca)

ggsave("plots/TCGA/TCGA_PCA_nes_plot.png",
       plot = tcga_pca, width = 10, height = 8, dpi = 300)

### OncoDB ###

# OncoDB GSEA results
oncodb_gsea_files <- list.files("stats/OncoDB/GSEA", pattern = "_gsea_results.csv$", full.names = TRUE)
oncodb_gsea_files <- oncodb_gsea_files[!grepl("OncoDB_combined_gsea_results.csv", oncodb_gsea_files)] # Exclude combined GSEA results file

# store NES values in a list
nes_matrix <- list()

for (file in oncodb_gsea_files) {
  project_id <- str_extract(basename(file), "^[^_]+")
  
  gsea_res <- read_csv(file) %>%
    select(pathway, NES) %>%
    filter(!is.na(NES))
  
  # use pathways as rownames
  rownames(gsea_res) <- gsea_res$pathway
  gsea_res$pathway <- NULL
  
  nes_matrix[[project_id]] <- gsea_res$NES
}

# Create NES matrix
nes_mat <- bind_rows(
  lapply(nes_matrix, function(x) as.data.frame(t(x))),
  .id = "project_id"
) %>%
  column_to_rownames("project_id") %>%
  as.matrix()

# Removing columns (i.e. Phylostrata) with more than 20% NAs
nes_mat_filtered <- nes_mat[, colMeans(is.na(nes_mat)) < 0.2]

# Imputing NAs
nes_mat_impute <- nes_mat # Copy of nes_mat to preserve original
nes_mat_impute[is.na(nes_mat_impute)] <- mean(nes_mat, na.rm = TRUE) # Replaced NAs with column mean

# Run PCA
pca_res <- prcomp(nes_mat_impute, scale. = TRUE)

# Plot
oncodb_pca <- fviz_pca_biplot(
  pca_res,
  label = "var",               
  habillage = rownames(nes_mat),  # color by project_id
  addEllipses = FALSE,
  repel = TRUE,
  geom = "point",
  pointshape = 20   
) +
  labs(
    title = "PCA of GSEA Normalised Enrichment Score by Cancer (OncoDB)",
    subtitle = "+ Additional Phylostrata contributions",
    x = paste0("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2), 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "%)")
  ) +
  theme_minimal()
show(oncodb_pca)

ggsave("plots/OncoDB/OncoDB_PCA_nes_plot.png",
       plot = tcga_pca, width = 10, height = 8, dpi = 300)
