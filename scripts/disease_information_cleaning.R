# Load required libraries
library(tidyverse)
library(glue)

# File paths
INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/metadata/Disease_information.txt"
OUTPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/Disease_information_cleaned.csv"

# Read tab-separated file
disease_raw <- read.delim(INPUT_FILE, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)
cat(glue("Loaded: {nrow(disease_raw)} rows, {ncol(disease_raw)} columns\n\n"))

# Standardise column names
colnames(disease_raw) <- c("disease_id", "disease_name", "omim_ids", "gene_accession_id")

# Trim whitespace from all character columns
disease_clean <- disease_raw %>%
  mutate(across(everything(), str_trim))

# Remove rows with missing values
disease_clean <- disease_clean %>%
  filter(
    !is.na(disease_id) & disease_id != "",
    !is.na(disease_name) & disease_name != "",
    !is.na(omim_ids) & omim_ids != "",
    !is.na(gene_accession_id) & gene_accession_id != ""
  )

# Remove duplicate rows
disease_clean <- disease_clean[!duplicated(disease_clean), ]

# Split OMIM IDs into separate rows
# Each OMIM ID gets its own row, keeping disease_id, disease_name, and gene_accession_id
disease_clean_normalized <- disease_clean %>%
  separate_rows(omim_ids, sep = "\\|") %>%   # split by pipe "|"
  mutate(omim_ids = str_trim(omim_ids))      # remove whitespace after splitting

# Validate IDs (optional)
invalid_disease_ids <- sum(!str_detect(disease_clean_normalized$disease_id, "^DOID:\\d+$"))
invalid_gene_ids <- sum(!str_detect(disease_clean_normalized$gene_accession_id, "^MGI:\\d+$"))
cat(glue("Invalid disease IDs: {invalid_disease_ids}\n"))
cat(glue("Invalid gene IDs: {invalid_gene_ids}\n"))

# Save normalized data
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
write.csv(disease_clean_normalized, OUTPUT_FILE, row.names = FALSE)

# Verify
if (file.exists(OUTPUT_FILE)) {
  file_size_mb <- file.info(OUTPUT_FILE)$size / (1024^2)
  cat(glue(" Normalized cleaned data saved: {OUTPUT_FILE}\n"))
  cat(glue("  Rows: {nrow(disease_clean_normalized)}, Size: {round(file_size_mb, 2)} MB\n"))
}
