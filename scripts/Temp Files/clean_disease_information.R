# Load required libraries
library(tidyverse)  # For data manipulation (dplyr, stringr)
library(glue)       # For formatted string output

# File paths
INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/metadata/Disease_information.txt"
OUTPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/cleaned_disease_information.csv"

# Read tab-separated file
# read.delim() is designed for tab-separated files (unlike read.csv() which expects commas)
# sep = "\t" specifies tab as delimiter
# quote = "\"" handles disease names containing commas within quotes
# stringsAsFactors = FALSE keeps text as character strings, not factors
disease_raw <- read.delim(INPUT_FILE, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)
cat(glue("Loaded: {nrow(disease_raw)} rows, {ncol(disease_raw)} columns\n\n"))

# Standardise column names
# Replace original column names with database-friendly names (no spaces, lowercase)
colnames(disease_raw) <- c("disease_id", "disease_name", "omim_ids", "gene_accession_id")

# Trim whitespace from all character columns
# str_trim() removes leading and trailing spaces from text
# Prevents issues like "DOID:123 " being different from "DOID:123"
disease_clean <- disease_raw %>%
  mutate(
    disease_id = str_trim(disease_id),
    disease_name = str_trim(disease_name),
    omim_ids = str_trim(omim_ids),
    gene_accession_id = str_trim(gene_accession_id)
  )

# Remove rows with missing values
# Track how many rows we have before filtering
rows_before <- nrow(disease_clean)

# Filter out rows where ANY critical field is NA or empty string
# All four fields are required for database integrity
disease_clean <- disease_clean %>%
  filter(
    !is.na(disease_id) & disease_id != "",
    !is.na(disease_name) & disease_name != "",
    !is.na(omim_ids) & omim_ids != "",
    !is.na(gene_accession_id) & gene_accession_id != ""
  )

# Calculate and report how many rows were removed
rows_removed <- rows_before - nrow(disease_clean)
if (rows_removed > 0) {
  cat(glue("Rows removed due to missing values: {rows_removed}\n"))
}
cat(glue("Rows after cleaning: {nrow(disease_clean)}\n\n"))

# Remove duplicate rows
# duplicated() returns TRUE for rows that are exact duplicates (after first occurrence)
# sum() counts how many duplicates exist
duplicates <- sum(duplicated(disease_clean))

if (duplicates > 0) {
  # Keep only unique rows (first occurrence of each duplicate)
  # !duplicated() inverts the boolean (FALSE becomes TRUE)
  disease_clean <- disease_clean[!duplicated(disease_clean), ]
  cat(glue("Duplicate rows removed: {duplicates}\n"))
}
cat(glue("Rows after removing duplicates: {nrow(disease_clean)}\n\n"))

# Validate ID formats
# Check disease IDs follow format: DOID:xxxxxxx (where x = digits)
# str_detect() with regex "^DOID:\\d+$" validates format
# ^ = start of string, \\d+ = one or more digits, $ = end of string
invalid_disease_ids <- sum(!str_detect(disease_clean$disease_id, "^DOID:\\d+$"))

# Check gene IDs follow format: MGI:xxxxxxx
invalid_gene_ids <- sum(!str_detect(disease_clean$gene_accession_id, "^MGI:\\d+$"))

# Report validation results
if (invalid_disease_ids > 0) {
  cat(glue("Warning: {invalid_disease_ids} invalid disease IDs found\n"))
} else {
  cat("✓ All disease IDs valid (DOID:xxxxxxx)\n")
}

if (invalid_gene_ids > 0) {
  cat(glue("Warning: {invalid_gene_ids} invalid gene IDs found\n"))
} else {
  cat("✓ All gene IDs valid (MGI:xxxxxxx)\n")
}
cat("\n")

# Save cleaned data
# Create output directory if it doesn't exist
# recursive = TRUE creates parent directories if needed
# showWarnings = FALSE prevents error messages if directory already exists
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)

# Write cleaned data to CSV
# row.names = FALSE prevents adding an extra column of row numbers
write.csv(disease_clean, OUTPUT_FILE, row.names = FALSE)

# Verify file was saved and report details
if (file.exists(OUTPUT_FILE)) {
  # Convert file size from bytes to megabytes (divide by 1024^2)
  file_size_mb <- file.info(OUTPUT_FILE)$size / (1024^2)
  cat(glue("✓ Cleaned data saved: {OUTPUT_FILE}\n"))
  cat(glue("  Rows: {nrow(disease_clean)}, Size: {round(file_size_mb, 2)} MB\n\n"))
}