setwd("~/Desktop/gene_age_project/")

# Set your directory
data_dir <- "OncoDB_data/"
an_data_dir <- "OncoDB_data_annotated/"  

# List the relevant files
files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)

# Loop through each file, count rows, and print
for (file in files) {
  data <- read.table(file, sep = "\t")
  cat(basename(file), ":", nrow(data), "rows\n")
}

files <- list.files(an_data_dir, pattern = "\\.csv$", full.names = TRUE)

# Loop through each file, count rows, and print
for (file in files) {
  data <- read.csv(file)
  cat(basename(file), ":", nrow(data), "rows\n")
}
