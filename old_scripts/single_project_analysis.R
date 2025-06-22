### Setup ###
setwd("/scratch/prj/bmb_phyl_chron_disease/esta_proj/")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(umap)
library(pheatmap)

# Mapping TCGA project IDs to cancer types (alphabetically sorted)
cancer_types <- c(
  "TCGA-ACC" = "Adrenocortical Carcinoma",
  "TCGA-BLCA" = "Bladder Urothelial Carcinoma",
  "TCGA-BRCA" = "Breast Invasive Carcinoma",
  "TCGA-CESC" = "Cervical Squamous Cell Carcinoma & Endocervical Adenocarcinoma",
  "TCGA-CHOL" = "Cholangiocarcinoma",
  "TCGA-COAD" = "Colon Adenocarcinoma",
  "TCGA-DLBC" = "Diffuse Large B-cell Lymphoma",
  "TCGA-ESCA" = "Esophageal Carcinoma",
  "TCGA-GBM" = "Glioblastoma Multiforme",
  "TCGA-HNSC" = "Head & Neck Squamous Cell Carcinoma",
  "TCGA-KICH" = "Kidney Chromophobe",
  "TCGA-KIRC" = "Kidney Renal Clear Cell Carcinoma",
  "TCGA-KIRP" = "Kidney Renal Papillary Cell Carcinoma",
  "TCGA-LAML" = "Acute Myeloid Leukemia",
  "TCGA-LGG" = "Low Grade Glioma",
  "TCGA-LIHC" = "Liver Hepatocellular Carcinoma",
  "TCGA-LUAD" = "Lung Adenocarcinoma",
  "TCGA-LUSC" = "Lung Squamous Cell Carcinoma",
  "TCGA-MESO" = "Mesothelioma",
  "TCGA-OV" = "Ovarian Serous Cystadenocarcinoma",
  "TCGA-PAAD" = "Pancreatic Adenocarcinoma",
  "TCGA-PCPG" = "Pheochromocytoma & Paraganglioma",
  "TCGA-PRAD" = "Prostate Adenocarcinoma",
  "TCGA-READ" = "Rectum Adenocarcinoma",
  "TCGA-SARC" = "Sarcoma",
  "TCGA-SKCM" = "Skin Cutaneous Melanoma",
  "TCGA-STAD" = "Stomach Adenocarcinoma",
  "TCGA-THCA" = "Thyroid Carcinoma",
  "TCGA-THYM" = "Thymoma",
  "TCGA-TGCT" = "Testicular Germ Cell Tumors",
  "TCGA-UCEC" = "Uterine Corpus Endometrial Carcinoma",
  "TCGA-UCS" = "Uterine Carcinosarcoma",
  "TCGA-UVM" = "Uveal Melanoma"
)

### TCGA DE Pipeline with Phylostrata, UMAP, Fisher's Test, Wilcoxon Heatmap ###

# Setup output directories
dir.create("plots", showWarnings = FALSE)
dir.create("stats", showWarnings = FALSE)
dir.create("DE_results", showWarnings = FALSE)

# Mapping TCGA project IDs to cancer types
cancer_types <- c("TCGA-BRCA" = "Breast Invasive Carcinoma")
project_id <- "TCGA-BRCA"
cancer_type_name <- cancer_types[project_id]

### Download & Prepare Data ###
query <- GDCquery(
  project = project_id,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", directory = "GDCdata")
data <- GDCprepare(query)

phylostrata <- read_csv("phylostrata/gene_phylostrata.txt") %>%
  rename(gene_name = GeneID)

combined_data <- data[, data$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]

gene_mapping <- as.data.frame(rowData(combined_data)) %>%
  mutate(gene_id = rownames(.)) %>%
  select(gene_id, gene_name) %>%
  filter(gene_name %in% phylostrata$gene_name)

combined_data <- combined_data[rownames(combined_data) %in% gene_mapping$gene_id, ]
gene_mapping <- gene_mapping %>% filter(gene_id %in% rownames(combined_data))

### DESeq2 Analysis ###
dds <- DESeqDataSetFromMatrix(countData = assay(combined_data),
                              colData = colData(combined_data),
                              design = ~ sample_type)
dds <- DESeq(dds)

res_df <- as.data.frame(results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))) %>%
  mutate(gene_id = rownames(.)) %>%
  left_join(gene_mapping, by = "gene_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
  select(gene_name, everything()) %>%
  left_join(phylostrata, by = "gene_name")

write_csv(res_df, paste0("DE_results/", project_id, "_DE.csv"))

### UMAP on DE statistics ###
res_stats <- res_df %>% drop_na(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
numeric_data <- res_stats %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
umap_result <- umap(scale(numeric_data))

umap_df <- as.data.frame(umap_result$layout) %>%
  mutate(gene_name = res_stats$gene_name,
         phylostrata = as.factor(res_stats$Phylostrata))

umap_plot <- ggplot(umap_df, aes(x = V1, y = V2, color = phylostrata)) +
  geom_point(alpha = 0.7, size = 1.5) +
  labs(title = paste("UMAP of DE Features for", cancer_type_name),
       x = "UMAP 1", y = "UMAP 2", color = "Phylostrata") +
  theme_minimal()

ggsave(paste0("plots/", project_id, "_umap.png"), plot = umap_plot, width = 8, height = 6)

### Fisher's Exact Test ###
res_df <- res_df %>%
  mutate(Regulation = case_when(
    padj < 0.05 & log2FoldChange >= 2 ~ "Up",
    padj < 0.05 & log2FoldChange <= -2 ~ "Down",
    TRUE ~ "NS"
  ),
  Phylostrata_Group = ifelse(Phylostrata <= 6, "Ancient", "Recent"))

deg_filtered <- res_df %>% filter(Regulation != "NS", !is.na(Phylostrata_Group))
contingency <- table(deg_filtered$Phylostrata_Group, deg_filtered$Regulation)
fisher_res <- fisher.test(contingency, simulate.p.value = TRUE, B = 10000)
write.csv(contingency, paste0("stats/", project_id, "_contingency.csv"))
capture.output(fisher_res, file = paste0("stats/", project_id, "_fisher_test.txt"))

### Pairwise Wilcoxon Test + Heatmap ###
res_df_clean <- res_df %>% filter(!is.na(log2FoldChange), !is.na(Phylostrata), padj < 0.001)
pairwise_wilcox <- pairwise.wilcox.test(
  res_df_clean$log2FoldChange,
  as.factor(res_df_clean$Phylostrata),
  p.adjust.method = "BH"
)
capture.output(pairwise_wilcox, file = paste0("stats/", project_id, "_pairwise_wilcox.txt"))

# Format to symmetric matrix
wilcox_matrix <- as.matrix(as.data.frame(pairwise_wilcox$p.value))
for (i in 1:nrow(wilcox_matrix)) {
  for (j in 1:ncol(wilcox_matrix)) {
    if (is.na(wilcox_matrix[i, j])) {
      wilcox_matrix[i, j] <- wilcox_matrix[j, i]
    }
  }
}

phylostrata_in_matrix <- colnames(wilcox_matrix)
rownames(wilcox_matrix) <- phylostrata_in_matrix
colnames(wilcox_matrix) <- phylostrata_in_matrix

rownames(wilcox_matrix) <- phylostrata_in_matrix
colnames(wilcox_matrix) <- phylostrata_in_matrix
rownames(wilcox_matrix) <- colnames(wilcox_matrix) <- levels(as.factor(res_df_clean$Phylostrata))

pheatmap(
  wilcox_matrix,
  display_numbers = TRUE,
  color = colorRampPalette(c("white", "red"))(100),
  main = paste("Pairwise Wilcoxon p-values (", cancer_type_name, ")")
)
ggsave(paste0("plots/", project_id, "_wilcox_heatmap.png"), width = 8, height = 6)

### Kruskal-Wallis Test ###
kruskal_test_result <- kruskal.test(log2FoldChange ~ as.factor(Phylostrata), data = res_df_clean)
message("Kruskal-Wallis p-value: ", kruskal_test_result$p.value)
