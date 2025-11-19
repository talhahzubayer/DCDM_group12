# Load required libraries
library(tidyverse)  # For data manipulation functions
library(glue)       # For string formatting in messages

#Input: CSV file that is generated from the consolidate_all_csv_files.R script
INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/consolidated_phenotype_analysis_results.csv"

# Output: the final cleaned data ready for database loading
OUTPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/cleaned_phenotype_analysis_results.csv"

# Reference files: rows that were removed during cleaning
DUPLICATES_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/removed_duplicates_from_consolidated_file.csv"
REMOVED_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/SOP_removed_rows_from_consolidated_file.csv"

# Read the consolidated CSV file
# stringsAsFactors = FALSE prevents automatic conversion to factors
data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

# The initial data dimensions
cat(glue("Loaded: {nrow(data)} rows, {ncol(data)} columns\n\n"))


# REMOVE DUPLICATES

# duplicated() returns TRUE for duplicate rows (after first occurrence)
# Example: if row 5 is identical to row 2, row 5 will be TRUE
duplicate_rows <- duplicated(data)

# Count how many duplicate rows exist
n_duplicates <- sum(duplicate_rows)

# Handle duplicates if they exist
if (n_duplicates > 0) {
  cat(glue("Duplicate rows found: {n_duplicates}\n"))
  
  # Keep only unique rows (first occurrence of each duplicate)
  # !duplicate_rows inverts the boolean vector (FALSE becomes TRUE)
  cleaned_data <- data[!duplicate_rows, ]
  
  # Saving duplicate rows for reference (to check what was removed)
  removed_duplicates <- data[duplicate_rows, ]
  write.csv(removed_duplicates, DUPLICATES_FILE, row.names = FALSE)
  
  cat(glue("✓ Saved removed duplicates to: {DUPLICATES_FILE}\n"))
  cat(glue("✓ Rows after removing duplicates: {nrow(cleaned_data)}\n\n"))
} else {
  # No duplicates found - proceed with all data
  cat("✓ No duplicate rows found - data is already unique!\n")
  cleaned_data <- data
  cat(glue("✓ Proceeding with all {nrow(cleaned_data)} rows\n\n"))
}

# REMOVE MISSING VALUES (NA and empty strings)

# Replace empty strings and text "NA" with proper NA values
# This standardises missing data representation
cleaned_data[cleaned_data == ""] <- NA
cleaned_data[cleaned_data == "NA"] <- NA

# Record row count before NA removal
rows_before_na <- nrow(cleaned_data)
cat(glue("Rows before removing NAs: {rows_before_na}\n\n"))

# na.omit() removes any row that contains at least one NA value - this ensures complete data for all analyses
cleaned_data <- na.omit(cleaned_data)


# Calculate how many rows were removed
rows_removed_after_na <- rows_before_na - nrow(cleaned_data) 

cat(glue("Rows removed due to NAs: {rows_removed_after_na}\n\n"))
cat(glue("Rows after removing NAs: {nrow(cleaned_data)}\n\n"))

# APPLY SOP VALIDATION RULES
# Following IMPC Standard Operating Procedure specifications

# Create a boolean vector to track valid rows
# Initially all TRUE (all rows assumed valid) - it will become FALSE for rows that fail any validation check
valid_rows <- rep(TRUE, nrow(cleaned_data))

cat("Validation checks:\n")

# Helper function: checks if string length is within range
check_length <- function(column, min_len, max_len) {
  
  # nchar() counts characters in each value
  # as.character() ensures we are working with strings
  nchar_val <- nchar(as.character(column))
  
  # Return TRUE if length is within range, FALSE otherwise
  return(nchar_val >= min_len & nchar_val <= max_len)
}


# SOP RULE 1: analysis_id must be exactly 15 characters
cat("Checking analysis_id (must be 15 chars)...\n")
valid_analysis_id <- nchar(as.character(cleaned_data$analysis_id)) == 15

# Update valid_rows: keep TRUE only if analysis_id is valid
# The & operator performs element-wise AND
valid_rows <- valid_rows & valid_analysis_id
cat(glue("    Invalid: {sum(!valid_analysis_id)}\n"))

# SOP RULE 2: gene_accession_id must be 9-11 characters
cat("\nChecking gene_accession_id (9-11 chars)...\n")
# Helper function is used here to check length range
valid_gene_acc <- check_length(cleaned_data$gene_accession_id, 9, 11)
valid_rows <- valid_rows & valid_gene_acc
cat(glue("    Invalid: {sum(!valid_gene_acc)}\n"))

# SOP RULE 3: gene_symbol must be 1-13 characters
cat("\nChecking gene_symbol (1-13 chars)...\n")
valid_gene_symbol <- check_length(cleaned_data$gene_symbol, 1, 13)
valid_rows <- valid_rows & valid_gene_symbol
cat(glue("    Invalid: {sum(!valid_gene_symbol)}\n"))

# SOP RULE 4: mouse_strain must be 3-5 characters
cat("\nChecking mouse_strain (3-5 chars)...\n")
valid_strain <- check_length(cleaned_data$mouse_strain, 3, 5)
valid_rows <- valid_rows & valid_strain
cat(glue("    Invalid: {sum(!valid_strain)}\n"))

# SOP RULE 5: mouse_life_stage must be 4-17 characters
cat("\nChecking mouse_life_stage (4-17 chars)...\n")
valid_life_stage <- check_length(cleaned_data$mouse_life_stage, 4, 17)
valid_rows <- valid_rows & valid_life_stage
cat(glue("    Invalid: {sum(!valid_life_stage)}\n"))

# SOP RULE 6: parameter_id must be 15-20 characters
cat("\nChecking parameter_id (15-20 chars)...\n")
valid_param_id <- check_length(cleaned_data$parameter_id, 15, 20)
valid_rows <- valid_rows & valid_param_id
cat(glue("    Invalid: {sum(!valid_param_id)}\n"))

# SOP RULE 7: parameter_name must be 2-74 characters
cat("\nChecking parameter_name (2-74 chars)...\n")
valid_param_name <- check_length(cleaned_data$parameter_name, 2, 74)
valid_rows <- valid_rows & valid_param_name
cat(glue("    Invalid: {sum(!valid_param_name)}\n"))

# SOP RULE 8: pvalue must be numeric and between 0 and 1
cat("\nChecking pvalue (numeric, 0-1)...\n")

# Convert pvalue column to numeric type
# suppressWarnings() hides warnings about non-numeric values becoming NA
# Any text values will become NA, which we'll catch in the next check
cleaned_data$pvalue <- suppressWarnings(as.numeric(cleaned_data$pvalue))

# Check 3 conditions:
# 1. Value is not NA (conversion was successful)
# 2. Value is >= 0 (minimum valid p-value)
# 3. Value is <= 1 (maximum valid p-value)
valid_pvalue <- !is.na(cleaned_data$pvalue) & 
  (cleaned_data$pvalue >= 0) & 
  (cleaned_data$pvalue <= 1)

valid_rows <- valid_rows & valid_pvalue
cat(glue("    Invalid: {sum(!valid_pvalue)}\n"))

cat("\n")


# final_cleaned: rows that passed all 8 SOP validation checks
final_cleaned <- cleaned_data[valid_rows, ]

# removed_rows: rows that failed at least one validation check
removed_rows <- cleaned_data[!valid_rows, ]

cat(glue("Rows passing SOP validation: {nrow(final_cleaned)}\n\n"))
cat(glue("Rows failing SOP validation: {nrow(removed_rows)}\n\n"))


# REMOVE OUTLIERS USING IQR METHOD

cat("Removing outliers using IQR method...\n")

# Identify numeric columns in the cleaned dataset
numeric_cols <- sapply(final_cleaned, is.numeric)

# Function to flag outliers using IQR method
# IQR (Interquartile Range) = Q3 - Q1
# Outliers are values outside: [Q1 - 1.5*IQR, Q3 + 1.5*IQR]
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  # Values within 1.5 * IQR from Q1 and Q3 are kept
  x >= (Q1 - 1.5 * IQR_val) & x <= (Q3 + 1.5 * IQR_val)
}

# Initialize a logical vector (all TRUE initially)
keep_rows <- rep(TRUE, nrow(final_cleaned))

# Apply outlier filter for each numeric column
for (col in names(final_cleaned)[numeric_cols]) {
  valid <- remove_outliers(final_cleaned[[col]])
  keep_rows <- keep_rows & valid
}

# Record counts before removing outliers
rows_before_outliers <- nrow(final_cleaned)

# Filter data to remove outliers
final_cleaned <- final_cleaned[keep_rows, ]

# Calculate and report outliers removed
outliers_removed <- rows_before_outliers - nrow(final_cleaned)
cat(glue("Outliers removed: {outliers_removed}\n"))
cat(glue("Rows after removing outliers: {nrow(final_cleaned)}\n\n"))


# CONVERT GENE_ACCESSION_ID TO UPPERCASE

cat("Converting gene_accession_id to uppercase...\n")
final_cleaned$gene_accession_id <- toupper(final_cleaned$gene_accession_id)
cat("✓ gene_accession_id converted to uppercase\n\n")


# SAVE CLEANED DATA

# Create output directory if needed
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)

# Save final cleaned data
write.csv(final_cleaned, OUTPUT_FILE, row.names = FALSE)

# Save removed rows (SOP violations) for reference
if (nrow(removed_rows) > 0) {
  write.csv(removed_rows, REMOVED_FILE, row.names = FALSE)
}

# Verify and report
if (file.exists(OUTPUT_FILE)) {
  file_size_mb <- file.info(OUTPUT_FILE)$size / (1024^2)
  cat(glue("✓ Cleaned data saved: {OUTPUT_FILE}\n"))
  cat(glue("  File size: {round(file_size_mb, 2)} MB\n\n"))
}

if (nrow(removed_rows) > 0 && file.exists(REMOVED_FILE)) {
  cat(glue("✓ Removed rows saved: {REMOVED_FILE}\n\n"))
}

if (n_duplicates > 0 && file.exists(DUPLICATES_FILE)) {
  cat(glue("✓ Duplicates saved: {DUPLICATES_FILE}\n\n"))
}