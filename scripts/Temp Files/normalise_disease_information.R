# ============================================================================
# Script: normalize_disease_information.R
# Purpose: Properly normalize disease data by splitting OMIM IDs
# This creates THREE normalized tables instead of one denormalized file
# ============================================================================

library(tidyverse)
library(glue)

# ============================================================================
# FILE PATHS - CORRECTED FOR YOUR FOLDER STRUCTURE
# ============================================================================

INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/cleaned_disease_information.csv"
OUTPUT_DISEASE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/normalised_disease_table.csv"
OUTPUT_OMIM <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/normalised_omim_table.csv"
OUTPUT_ASSOCIATION <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/normalised_gene_disease_association.csv"

cat("============================================\n")
cat("NORMALIZING DISEASE INFORMATION\n")
cat("Splitting OMIM IDs into separate table\n")
cat("============================================\n\n")

# Check if input file exists
if (!file.exists(INPUT_FILE)) {
  cat("ERROR: Input file not found at:\n")
  cat(INPUT_FILE, "\n")
  stop("Please check the file path")
}

# Read the cleaned disease file
disease_raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
cat(glue("Loaded: {nrow(disease_raw)} rows, {ncol(disease_raw)} columns\n\n"))

# Verify column structure
cat("Columns found:\n")
print(colnames(disease_raw))
cat("\n")

# Show example of the normalization issue
cat("Example of non-normalized data (multiple OMIM IDs):\n")
example <- disease_raw %>% 
  filter(str_detect(omim_ids, "\\|")) %>%
  head(3) %>%
  select(disease_name, omim_ids)
print(example)
cat("\n")

multi_omim_count <- sum(str_detect(disease_raw$omim_ids, "\\|"))
cat(glue("Rows with multiple OMIM IDs: {multi_omim_count} ({round(multi_omim_count/nrow(disease_raw)*100, 1)}%)\n\n"))

# ============================================================================
# TABLE 1: Disease (unique diseases only, no OMIM IDs)
# ============================================================================

cat("Creating Disease table (unique diseases)...\n")

disease_table <- disease_raw %>%
  select(disease_id, disease_name) %>%
  distinct(disease_id, .keep_all = TRUE)

cat(glue("  Unique diseases: {nrow(disease_table)}\n\n"))

# ============================================================================
# TABLE 2: OMIM (split pipe-separated OMIM IDs)
# ============================================================================

cat("Creating OMIM table (splitting pipe-separated values)...\n")

# Function to split OMIM IDs and create rows
split_omim_ids <- function(disease_id, omim_string) {
  # Split by pipe character
  omim_list <- unlist(strsplit(omim_string, "\\|"))
  
  # Create a data frame with one row per OMIM ID
  data.frame(
    omim_id = trimws(omim_list),
    disease_id = disease_id,
    stringsAsFactors = FALSE
  )
}

# Apply to all rows
omim_table <- disease_raw %>%
  select(disease_id, omim_ids) %>%
  distinct() %>%
  # Split each row into multiple rows (one per OMIM ID)
  rowwise() %>%
  do(split_omim_ids(.$disease_id, .$omim_ids)) %>%
  ungroup() %>%
  # Remove duplicates
  distinct()

cat(glue("  Total OMIM entries: {nrow(omim_table)}\n"))
cat(glue("  Unique OMIM IDs: {n_distinct(omim_table$omim_id)}\n\n"))

# Show example of normalized data
cat("Example of normalized OMIM data:\n")
cat("BEFORE (denormalized):\n")
example_before <- disease_raw %>%
  filter(str_detect(omim_ids, "\\|")) %>%
  head(1) %>%
  select(disease_id, disease_name, omim_ids)
print(example_before)

cat("\nAFTER (normalized):\n")
example_after <- omim_table %>%
  filter(disease_id == example_before$disease_id[1])
print(example_after)
cat("\n")

# ============================================================================
# TABLE 3: GeneDiseaseAssociation (many-to-many)
# ============================================================================

cat("Creating GeneDiseaseAssociation table...\n")

gene_disease_association <- disease_raw %>%
  select(gene_accession_id, disease_id) %>%
  distinct()

cat(glue("  Gene-disease associations: {nrow(gene_disease_association)}\n\n"))

# ============================================================================
# STATISTICS
# ============================================================================

cat("============================================\n")
cat("NORMALIZATION STATISTICS\n")
cat("============================================\n\n")

cat("BEFORE (Denormalized):\n")
cat(glue("  1 table with {nrow(disease_raw)} rows\n"))
cat(glue("  Multiple OMIM IDs in {multi_omim_count} rows ({round(multi_omim_count/nrow(disease_raw)*100, 1)}%)\n"))
cat("  Violates 1NF (First Normal Form)\n\n")

cat("AFTER (Normalized - 3NF):\n")
cat(glue("  Disease table: {nrow(disease_table)} unique diseases\n"))
cat(glue("  OMIM table: {nrow(omim_table)} OMIM-disease associations\n"))
cat(glue("  GeneDiseaseAssociation: {nrow(gene_disease_association)} gene-disease associations\n\n"))

# Calculate average OMIM IDs per disease
avg_omim <- omim_table %>%
  group_by(disease_id) %>%
  summarise(count = n()) %>%
  summarise(
    min = min(count),
    max = max(count),
    avg = round(mean(count), 2),
    median = median(count)
  )

cat("OMIM IDs per disease:\n")
cat(glue("  Min: {avg_omim$min}, Max: {avg_omim$max}, Avg: {avg_omim$avg}, Median: {avg_omim$median}\n\n"))

# ============================================================================
# VALIDATION
# ============================================================================

cat("============================================\n")
cat("DATA QUALITY VALIDATION\n")
cat("============================================\n\n")

# Check for duplicate disease IDs
dup_diseases <- disease_table %>%
  group_by(disease_id) %>%
  summarise(count = n()) %>%
  filter(count > 1)

if (nrow(dup_diseases) > 0) {
  cat("⚠ Warning: Duplicate disease IDs found\n")
} else {
  cat("✓ No duplicate disease IDs\n")
}

# Check for duplicate OMIM IDs
dup_omim <- omim_table %>%
  group_by(omim_id, disease_id) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(count > 1)

if (nrow(dup_omim) > 0) {
  cat("⚠ Warning: Duplicate OMIM-disease pairs found\n")
} else {
  cat("✓ No duplicate OMIM-disease pairs\n")
}

# Check for duplicate gene-disease associations
dup_assoc <- gene_disease_association %>%
  group_by(gene_accession_id, disease_id) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(count > 1)

if (nrow(dup_assoc) > 0) {
  cat("⚠ Warning: Duplicate gene-disease associations found\n")
} else {
  cat("✓ No duplicate gene-disease associations\n")
}

cat("\n")

# Validate OMIM ID format
invalid_omim <- omim_table %>%
  filter(!str_detect(omim_id, "^OMIM:\\d+$"))

if (nrow(invalid_omim) > 0) {
  cat(glue("⚠ Warning: {nrow(invalid_omim)} invalid OMIM IDs found\n"))
  cat("Examples:\n")
  print(head(invalid_omim, 5))
} else {
  cat("✓ All OMIM IDs valid format (OMIM:xxxxxxx)\n")
}

cat("\n")

# ============================================================================
# SAVE NORMALIZED TABLES
# ============================================================================

cat("============================================\n")
cat("SAVING NORMALIZED TABLES\n")
cat("============================================\n\n")

# Save Disease table
write.csv(disease_table, OUTPUT_DISEASE, row.names = FALSE)
cat(glue("✓ Disease table saved: {OUTPUT_DISEASE}\n"))
cat(glue("  Rows: {nrow(disease_table)}\n\n"))

# Save OMIM table
write.csv(omim_table, OUTPUT_OMIM, row.names = FALSE)
cat(glue("✓ OMIM table saved: {OUTPUT_OMIM}\n"))
cat(glue("  Rows: {nrow(omim_table)}\n\n"))

# Save GeneDiseaseAssociation table
write.csv(gene_disease_association, OUTPUT_ASSOCIATION, row.names = FALSE)
cat(glue("✓ GeneDiseaseAssociation table saved: {OUTPUT_ASSOCIATION}\n"))
cat(glue("  Rows: {nrow(gene_disease_association)}\n\n"))

# ============================================================================
# SUMMARY
# ============================================================================

cat("============================================\n")
cat("NORMALIZATION COMPLETE!\n")
cat("============================================\n\n")

cat("Three normalized CSV files created:\n")
cat(glue("  1. {basename(OUTPUT_DISEASE)} ({nrow(disease_table)} rows)\n"))
cat(glue("  2. {basename(OUTPUT_OMIM)} ({nrow(omim_table)} rows)\n"))
cat(glue("  3. {basename(OUTPUT_ASSOCIATION)} ({nrow(gene_disease_association)} rows)\n\n"))

cat("Database Schema (3NF):\n")
cat("  Disease (disease_id [PK], disease_name)\n")
cat("  OMIM (omim_id [PK], disease_id [FK])\n")
cat("  GeneDiseaseAssociation (gene_accession_id [FK], disease_id [FK])\n\n")

cat("Next steps:\n")
cat("  1. Load these 3 CSV files into your MySQL database\n")
cat("  2. Create foreign key constraints\n")
cat("  3. Mention normalization in your report\n\n")

cat("For your report:\n")
cat('  "We identified a First Normal Form violation where OMIM IDs were\n')
cat('   stored as pipe-separated values. We normalized this by creating a\n')
cat('   separate OMIM table, achieving Third Normal Form (3NF)."\n\n')