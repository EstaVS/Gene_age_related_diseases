setwd("~/Desktop/gene_age_project/")

### Adding ensembl ID to phylostrata table ###

# Load both files
phylostrata <- read.table("phylostrata/gene_phylostrata.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)
metadata <- read.table("metadata/gene_id_metadata.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# View column names to match (adjust if needed)
head(phylostrata)
head(metadata)

# Modifying metadata column names
colnames(metadata) <- c("EnsemblID", "GeneID", "Entrez")

### Merge via GeneID ###
merged <- merge(phylostrata, metadata, by = "GeneID")

sum(merged$Entrez.x == merged$Entrez.y)
sum(merged$Entrez.x != merged$Entrez.y)

# Explore dataframe to find with Entrez column to remove
length(unique(merged$Entrez.y)) # More unique Entrez
length(unique(merged$Entrez.x))

# Remove the Entrez.x column
merged$Entrez.x <- NULL

# Rename Entrez.y to Entrez
colnames(merged)[colnames(merged) == "Entrez.y"] <- "Entrez"

### Merge based on Entrez ###
merged_2 <- merge(phylostrata, metadata, by = "Entrez")

sum(merged_2$GeneID.x == merged_2$GeneID.y)
sum(merged_2$GeneID.x != merged_2$GeneID.y)

length(unique(merged_2$GeneID.y))
length(unique(merged_2$GeneID.x))

# Explore merged_2 mismatched
library(dplyr)
merged_2 %>% filter(GeneID.x != GeneID.y) # GeneID.x includes less meaninful Gene symbols

# Remove the Entrez.x column
merged_2$GeneID.x <- NULL

# Rename Entrez.y to Entrez
colnames(merged_2)[colnames(merged_2) == "GeneID.y"] <- "GeneID"

### Check uniqueness of EnsemblID ###

sum(duplicated(merged$EnsemblID))
sum(duplicated(merged_2$EnsemblID))

### Maximising EnsemblID uniqueness
# Removing duplicates & unmatched Entrez or GeneIDs

# Remove duplicated Ensembl IDs
merged_filtered <- merged[!duplicated(merged$EnsemblID), ]
length(merged_filtered$EnsemblID)

merged2_filtered <- merged_2[!duplicated(merged_2$EnsemblID), ]
length(merged2_filtered$EnsemblID)

# Identify Ensembl IDs in merged2_filtered that are NOT in merged_filtered
missing_ensembl_ids <- setdiff(merged2_filtered$EnsemblID, merged_filtered$EnsemblID)

# Subset those rows from merged2_filtered
additional_rows <- merged2_filtered[merged2_filtered$EnsemblID %in% missing_ensembl_ids, ]

# Combine with merged_filtered
merged_combined <- rbind(merged_filtered, additional_rows)

# Save combined file
write.table(merged_combined, "phylostrata_ensembl.txt", sep = ",", row.names = FALSE, quote = FALSE)
