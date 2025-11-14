# Load required libraries
library(tidyverse)  # For data manipulation (dplyr, tidyr, stringr)
library(glue)       # For formatted string output

# File paths
INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/metadata/Disease_information.txt"
OUTPUT_DISEASE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/cleaned_disease_table.csv"
OUTPUT_ASSOCIATION <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/cleaned_gene_disease_association.csv"

# Read tab-separated file
# read.delim() is designed for tab-separated files (unlike read.csv() which expects commas)
# quote = "\"" handles disease names that contain commas within quotes
# stringsAsFactors = FALSE keeps text as character strings, not factors
disease_raw <- read.delim(INPUT_FILE, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)
cat(glue("Loaded: {nrow(disease_raw)} rows, {ncol(disease_raw)} columns\n\n"))

# Standardise column names
# Replace original column names with database-friendly names (no spaces, lowercase)
colnames(disease_raw) <- c("disease_id", "disease_name", "omim_ids", "gene_accession_id")

# Trim whitespace from all character columns
# str_trim() removes leading and trailing spaces from text
# This prevents issues like "DOID:123 " being different from "DOID:123"
disease_clean <- disease_raw %>%
  mutate(
    disease_id = str_trim(disease_id),
    disease_name = str_trim(disease_name),
    omim_ids = str_trim(omim_ids),
    gene_accession_id = str_trim(gene_accession_id)
  )

# Remove rows with missing values
# Track how many rows we started with
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

# Analyse multiple OMIM IDs
# Some diseases have multiple OMIM IDs separated by pipe (|) symbol
# Example: "OMIM:101800|OMIM:614613" represents genetic heterogeneity
# (same gene causing different disease variants)
multi_omim_stats <- disease_clean %>%
  mutate(has_multiple_omim = str_detect(omim_ids, "\\|")) %>%  # Check if pipe exists
  summarise(
    total_rows = n(),
    rows_with_multiple_omim = sum(has_multiple_omim),
    percentage = round(sum(has_multiple_omim) / n() * 100, 2)
  )

cat(glue("Rows with multiple OMIM IDs: {multi_omim_stats$rows_with_multiple_omim} ({multi_omim_stats$percentage}%)\n\n"))

# Create Disease table (unique diseases)
# This table stores unique disease information with disease_id as primary key
# distinct() keeps only the first occurrence of each disease_id
# This handles cases where the same disease appears multiple times with different genes
disease_table <- disease_clean %>%
  select(disease_id, disease_name, omim_ids) %>%
  distinct(disease_id, .keep_all = TRUE)

cat(glue("Unique diseases: {nrow(disease_table)}\n\n"))

# Create GeneDiseaseAssociation table (many-to-many)
# This table represents the relationship between diseases and genes
# One disease can be associated with multiple genes (genetic heterogeneity)
# One gene can be associated with multiple diseases (pleiotropy)
# distinct() removes any duplicate gene-disease pairs
gene_disease_association <- disease_clean %>%
  select(disease_id, gene_accession_id) %>%
  distinct()

cat(glue("Gene-disease associations: {nrow(gene_disease_association)}\n\n"))

# Relationship statistics
# Calculate how many genes are associated with each disease
genes_per_disease <- gene_disease_association %>%
  group_by(disease_id) %>%              # Group by disease
  summarise(gene_count = n()) %>%       # Count genes per disease
  summarise(                            # Calculate summary statistics
    min_genes = min(gene_count),
    max_genes = max(gene_count),
    avg_genes = round(mean(gene_count), 2),
    median_genes = median(gene_count)
  )

# Calculate how many diseases are associated with each gene
diseases_per_gene <- gene_disease_association %>%
  group_by(gene_accession_id) %>%       # Group by gene
  summarise(disease_count = n()) %>%    # Count diseases per gene
  summarise(                            # Calculate summary statistics
    min_diseases = min(disease_count),
    max_diseases = max(disease_count),
    avg_diseases = round(mean(disease_count), 2),
    median_diseases = median(disease_count)
  )

# Display statistics
# sep = "" concatenates all elements without adding spaces
cat("Genes per disease: min=", genes_per_disease$min_genes, 
    ", max=", genes_per_disease$max_genes, 
    ", avg=", genes_per_disease$avg_genes, 
    ", median=", genes_per_disease$median_genes, "\n", sep = "")
cat("Diseases per gene: min=", diseases_per_gene$min_diseases, 
    ", max=", diseases_per_gene$max_diseases, 
    ", avg=", diseases_per_gene$avg_diseases, 
    ", median=", diseases_per_gene$median_diseases, "\n\n", sep = "")

# Validate ID formats
# Check that disease IDs follow the correct format: DOID:xxxxxxx (where x = digits)
# str_detect() with regex pattern "^DOID:\\d+$" checks format
# ^ = start of string, \\d+ = one or more digits, $ = end of string
invalid_disease_ids <- disease_table %>%
  filter(!str_detect(disease_id, "^DOID:\\d+$"))

# Check that gene accession IDs follow the correct format: MGI:xxxxxxx
invalid_gene_ids <- gene_disease_association %>%
  filter(!str_detect(gene_accession_id, "^MGI:\\d+$"))

# validation results
if (nrow(invalid_disease_ids) > 0) {
  cat(glue("Warning: {nrow(invalid_disease_ids)} invalid disease IDs found\n"))
} else {
  cat("✓ All disease IDs valid (DOID:xxxxxxx)\n")
}

if (nrow(invalid_gene_ids) > 0) {
  cat(glue("Warning: {nrow(invalid_gene_ids)} invalid gene IDs found\n"))
} else {
  cat("✓ All gene IDs valid (MGI:xxxxxxx)\n")
}
cat("\n")

# Save cleaned data
# Create output directory if it doesn't exist
# recursive = TRUE creates parent directories if needed
# showWarnings = FALSE prevents error messages if directory already exists
dir.create(dirname(OUTPUT_DISEASE), recursive = TRUE, showWarnings = FALSE)

# Write both tables to CSV files
# row.names = FALSE prevents adding an extra column of row numbers
write.csv(disease_table, OUTPUT_DISEASE, row.names = FALSE)
write.csv(gene_disease_association, OUTPUT_ASSOCIATION, row.names = FALSE)

# Verify files were saved and report file sizes
if (file.exists(OUTPUT_DISEASE)) {
  file_size_mb <- file.info(OUTPUT_DISEASE)$size / (1024^2)  # Convert bytes to MB
  cat(glue("✓ Disease table saved: {OUTPUT_DISEASE}\n"))
  cat(glue("  Rows: {nrow(disease_table)}, Size: {round(file_size_mb, 2)} MB\n\n"))
}

if (file.exists(OUTPUT_ASSOCIATION)) {
  file_size_mb <- file.info(OUTPUT_ASSOCIATION)$size / (1024^2)  # Convert bytes to MB
  cat(glue("✓ Gene-disease association saved: {OUTPUT_ASSOCIATION}\n"))
  cat(glue("  Rows: {nrow(gene_disease_association)}, Size: {round(file_size_mb, 2)} MB\n\n"))
}