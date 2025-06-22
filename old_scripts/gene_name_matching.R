setwd("~/Desktop/gene_age_project/")

phylostrata <- read.table("phylostratum_database/gene_phylostrata.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
raw_data <- read.csv("differential_expression_data.csv")
metadata <- read.table("metadata/human.txt", sep = ",", header = TRUE, stringsAsFactors = FALSE)
# Metadata is a file obtained from Ensembl Biomart API tool

### Defining vectors to compare across dataframes ###

# Gene names
phylostrata_genes <- unique(phylostrata$GeneID) # 17318
rawdata_genes <- unique(raw_data$gene_name) # 59427
metadata_genes <- unique(metadata$Gene.name) # 41150

# Entrez ID
phylo_entrez <- unique(phylostrata$Entrez) # 16734
metadata_entrez <- unique(metadata$NCBI.gene..formerly.Entrezgene..ID) # 28331

# Ensembl ID
metadata_ensembl <- unique(metadata$Gene.stable.ID) # 86402
rawdata_ensembl <- unique(raw_data$gene_id)
head(rawdata_ensembl) # IDs have .00 at the end which shows the version number - we need to remove these
# Remove version numbers from rawdata_ensembl (everything after the first dot)
rawdata_ensembl <- sub("\\..*", "", rawdata_ensembl)
# 60660

### Comparing Gene Names ###

# Raw data vs Phylostrata
print(length(intersect(rawdata_genes, phylostrata_genes))) # Gene names in both
print(length(setdiff(rawdata_genes, phylostrata_genes))) # Gene names in Raw data but not Phylostrata
print(length(setdiff(phylostrata_genes, rawdata_genes))) # Gene names in Phylostrata but not Raw data

# Raw data vs Metadata
print(length(intersect(metadata_genes, rawdata_genes))) # Genes in both
print(length(setdiff(metadata_genes, rawdata_genes))) # Genes in Metadata but not Raw data
print(length(setdiff(rawdata_genes, metadata_genes))) # Genes in Raw data but not Metadata

# Metadata vs Phylostrata
print(length(intersect(metadata_genes, phylostrata_genes))) # Genes in both
print(length(setdiff(metadata_genes, phylostrata_genes))) # Genes in Metadata but not Phylostrata
print(length(setdiff(phylostrata_genes, metadata_genes))) # Genes in Phylostrata but not Metadata

### Comparing Entrez ID - only Metadata vs Phylostrata

print(length(intersect(metadata_entrez, phylo_entrez))) # Genes in both
print(length(setdiff(metadata_entrez, phylo_entrez))) # Genes in Metadata but not Phylostrata
print(length(setdiff(phylo_entrez, metadata_entrez))) # Genes in Phylostrata but not Metadata

### Comparing Ensembl ID - only Raw data vs Metadata

print(length(intersect(metadata_ensembl, rawdata_ensembl))) # Genes in both
print(length(setdiff(metadata_ensembl, rawdata_ensembl))) # Genes in Metadata but not Raw data
print(length(setdiff(rawdata_ensembl, metadata_ensembl))) # Genes in Raw data but not Metadata

### Applying some gene name mutations based on the HGNC Multi-symbol checking tool data

# Define variables for mismatching genes names across Raw data & Metadata so we can investigate manually
missing_genes <- setdiff(expression_genes, metadata_gene_names)
head(missing_genes)
write.csv(missing_genes, "missing_genes.csv", row.names = FALSE, quote = FALSE)

# Loading HGNC Gene Mismatch table
HGNC_genenames <- read.csv("gene_mismatch/Matched_HGNC_genes.csv")
head(HGNC_genenames)

# Applying new gene names (approx 900 gene name changes)
matched_genes <- all_genes %>%
  left_join(HGNC_genenames, by = c("gene_name" = "Input"), relationship = "many-to-many") %>%
  mutate(gene_updated = ifelse(!is.na(Approved.symbol), Approved.symbol, gene_name)) %>%
  select(gene_updated, everything())  # Puts gene_updated in front

matched_genes <- matched_genes %>%
  select(-Match.type, -Approved.symbol, -Approved.name, -HGNC.ID, -Location)

matched_genes <- matched_genes %>%
  select(-gene_name)

matched_genes <- matched_genes %>%
  mutate(gene_name = gene_updated) %>%
  select(gene_name, everything()) %>%
  select(-gene_updated)

print(head(matched_genes))

# Check gene names with new changes

fixed_gene_names <- unique(matched_genes$gene_name)
metadata_gene_names <- unique(metadata$Gene.name)

# Raw data vs Metadata

print(length(intersect(metadata_genes, fixed_gene_names))) # Genes in both
print(length(setdiff(metadata_genes, fixed_gene_names))) # Genes in Metadata but not fixed Raw data
print(length(setdiff(fixed_gene_names, metadata_genes))) # Genes in fixed Raw data but not Metadata

# Old ones to compare - not a great improvement. Only 600 more matches
print(length(intersect(metadata_genes, rawdata_genes)))
print(length(setdiff(metadata_genes, rawdata_genes)))
print(length(setdiff(rawdata_genes, metadata_genes)))

# Raw data vs Phylostrata - fewer genes matching now

print(length(intersect(phylostrata_genes, fixed_gene_names))) # Genes in both
print(length(setdiff(phylostrata_genes, fixed_gene_names))) # Genes in Phylostrata but not fixed Raw data
print(length(setdiff(fixed_gene_names, phylostrata_genes))) # Genes in fixed Raw data but not Phylostrata

print(length(intersect(rawdata_genes, phylostrata_genes))) # Gene names in both
print(length(setdiff(rawdata_genes, phylostrata_genes))) # Gene names in old Raw data but not Phylostrata
print(length(setdiff(phylostrata_genes, rawdata_genes))) # Gene names in Phylostrata but not old Raw data

### Final decision to just exclude genes not in Phylostrata
