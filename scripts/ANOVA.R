setwd("~/Desktop/gene_age_project/")
library(stringr)
library(tibble)
library(dplyr)
library(readr)
library(broom)
library(purrr)
library(tidyverse)

### Analysis of Variance ###

# Combine all OncoDB DE results
oncodb_de_files <- list.files("OncoDB_data_annotated/", pattern = "expression_diff_.*\\.csv", full.names = TRUE)

all_oncodb_de <- lapply(oncodb_de_files, function(file) {
  project_id <- str_extract(basename(file), "(?<=expression_diff_)[^\\.]+")
  read_csv(file, show_col_types = FALSE) %>%
    mutate(
      Entrez = as.character(Entrez),
      project_id = project_id
    )
})

combined_oncodb_de <- bind_rows(all_oncodb_de)

# Global ANOVA Model - OncoDB only

global_model <- lm(log2FoldChange ~ Phylostrata + project_id, data = combined_oncodb_de)
anova(global_model) # Significant but no strong biological signal

# Integrating TCGA data

combined_tcga_de <- read_csv("TCGA_DE_results/all_TCGAprojects_combined_DE.csv", show_col_types = FALSE)
combined_tcga_de <- combined_tcga_de %>%
  mutate(project_id = str_remove(project_id, "^TCGA-")) # Removing the "TCGA-" section from the project IDs

# Global ANOVA Model - TCGA only
global_model <- lm(log2FoldChange ~ Phylostrata + project_id, data = combined_tcga_de)
anova(global_model) # Significant but also weak signal

### Combined ANOVA ###

# Add source labels
combined_oncodb_de$source <- "OncoDB"
combined_tcga_de$source   <- "TCGA"

# Combine both datasets
combined_de_all <- bind_rows(combined_oncodb_de, combined_tcga_de)

# Remove unnecessary columns
combined_de_all <- combined_de_all %>%
  select(gene_name, padj, cancer_median, normal_median, log2FoldChange,
         Entrez, Phylostrata, project_id, source)

# Inspecting the source weightings
table(combined_de_all$source)
obs_counts <- table(combined_de_all$source)
prop.table(obs_counts)

# Applying a source_weight column
combined_de_all <- combined_de_all %>%
  mutate(source_weight = 1 / as.numeric(table(source)[source]))

# Verify total weight by source
combined_de_all %>%
  group_by(source) %>%
  summarise(n_genes = n(),
            total_weight = sum(source_weight),
            prop_weight = total_weight / sum(source_weight))

## Define a weighted linear model
weighted_model <- lm(
  log2FoldChange ~ Phylostrata + project_id + source,
  data = combined_de_all,
  weights = source_weight
)

anova(weighted_model)

### Cancer Specific ANOVA ###

# Single source ANOVA and combined source ANOVA

# Get the unique cancers
cancers <- unique(combined_de_all$project_id)

for (cancer in cancers) {
  
  # Subset for this cancer
  subset_cancer <- combined_de_all %>% filter(project_id == cancer)
  
  # List the sources found
  sources <- unique(subset_cancer$source)
  
  # Display project_id and source(s)
  cat("\n------", cancer, "------\n")
  cat("Source(s):", paste(sources, collapse = ", "), "\n")
  
  # If only 1 source (either OncoDB or TCGA)
  if (length(sources) == 1) {
    cat("Running ANOVA: log2FoldChange ~ Phylostrata\n")
    
    model <- lm(log2FoldChange ~ Phylostrata, data = subset_cancer)
    anova_res <- anova(model)
    print(anova_res)
    
  } else if (length(sources) > 1) {
    cat("Running ANOVA: log2FoldChange ~ Phylostrata + source\n")
    
    model <- lm(log2FoldChange ~ Phylostrata + source, data = subset_cancer)
    anova_res <- anova(model)
    print(anova_res)
  }
}

### Visualisation ###
anova_summary <- read.csv("anova_summary.csv")
show(anova_summary)

# Pivot for stacked bar:
anova_long <- anova_summary %>%
  pivot_longer(
    cols = c("Phylostrata_pct", "Source_pct", "Residual_pct"),
    names_to = "Factor",
    values_to = "Proportion"
  ) %>%
  mutate(
    Factor = case_when(
      Factor == "Phylostrata_pct" ~ "Phylostrata",
      Factor == "Source_pct" ~ "Data_Source",
      Factor == "Residual_pct" ~ "Residual",
      TRUE ~ Factor
    )
  )%>%
  mutate(
    Factor = factor(Factor, levels = c("Residual", "Phylostrata", "Data_Source"))
  )

anova_plot <- ggplot(anova_long, aes(x = Cancer, y = Proportion, fill = Factor)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual(
    values = c(
      "Data_Source" = "royalblue",
      "Phylostrata" = "firebrick2",
      "Residual" = "palegreen"
    )) +
  labs(y = "Proportion of Variance Explained",
       title = "ANOVA Variance Partitioning by Cancer",
       fill = "Factor") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
show(anova_plot)

ggsave(paste0("plots/anova.png"),plot = anova_plot, width = 8, height = 5, dpi = 300)
