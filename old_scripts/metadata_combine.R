### Combining all metadata files into one and removing duplicate gene information ###

setwd("~/Desktop/gene_age_project/metadata/")
library(dplyr)

# Load the files
file1 <- read.table("human.txt", sep = ",", header = TRUE, stringsAsFactors = FALSE)
file2 <- read.table("mouse.txt", sep = ",", header = TRUE, stringsAsFactors = FALSE)
file3 <- read.table("chicken.txt", sep = ",", header = TRUE, stringsAsFactors = FALSE)
file4 <- read.table("zebrafish.txt", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Combine all files
combined_data <- rbind(file1, file2, file3, file4)

unique_data <- combined_data %>%
  distinct(NCBI.gene..formerly.Entrezgene..ID, .keep_all = TRUE)
  # Keep the first occurrence of each Entrez gene ID

# Save the combined and deduplicated data
write.table(unique_data, file = "gene_id_metadata.txt", sep = ",", row.names = FALSE, quote = FALSE)

### Didn't improve results but keeping script for workflow proof ###