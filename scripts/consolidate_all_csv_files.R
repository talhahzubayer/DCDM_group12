# Load required libraries
library(tidyverse)  # For data manipulation functions
library(glue)       # For string formatting in messages

# Directory containing original CSV files (vertical format)
CSV_DIR <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/raw/phenotype_analysis_results"

# Directory to save transposed CSV files (horizontal format)
TRANSPOSED_CSV_DIR <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/raw/phenotype_analysis_results_transposed"

# Final output file path
OUTPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/consolidated_phenotype_analysis_results.csv"

# Create transposed directory if it doesn't exist
# recursive = TRUE creates parent directories if needed
dir.create(TRANSPOSED_CSV_DIR, recursive = TRUE, showWarnings = FALSE)


# Transposing all CSV files

# Find all CSV files in the original directory
# pattern = "\\.csv$" matches files ending in .csv
# full.names = TRUE returns complete file paths
files <- list.files(path = CSV_DIR, pattern = "\\.csv$", full.names = TRUE)
n_files <- length(files)

cat(glue("Found {n_files} CSV files to transpose\n\n"))

# Loop to transpose each CSV file
for (i in seq_along(files)) {
  file <- files[i]
  
  # Read CSV without assuming a header row
  # Original format has key-value pairs in 2 columns
  data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE)
  
  # Transpose: rotate data 90 degrees
  # Converts vertical key-value pairs to horizontal columns
  transposed <- t(data)
  
  # Convert matrix back to data frame
  transposed_df <- as.data.frame(transposed, stringsAsFactors = FALSE)
  
  # Use first row as column names
  # First row contains the keys (gene_accession_id, gene_symbol, etc.)
  colnames(transposed_df) <- transposed_df[1, ]
  
  # Remove the first row (now redundant as column names)
  transposed_df <- transposed_df[-1, ]
  
  # Reset row names to NULL (cleaner output)
  rownames(transposed_df) <- NULL
  
  # Create new filename with "transposed_" prefix
  # basename() extracts filename from full path
  # file.path() creates proper path for the new directory
  new_name <- file.path(TRANSPOSED_CSV_DIR, paste0("transposed_", basename(file)))
  
  # Save the transposed file to new directory
  write.csv(transposed_df, new_name, row.names = FALSE)
  
  # Display progress counter (updates on same line)
  # \r returns cursor to start of line
  # flush.console() forces immediate display update
  cat(glue("\rTransposed: {i}/{n_files}"))
  flush.console()
}

cat("\n\n✓ All CSV files have been transposed!\n\n")

# Verify the number of transposed files created
n_transposed <- length(list.files(path = TRANSPOSED_CSV_DIR, pattern = "^transposed_.*\\.csv$"))
cat(glue("Total transposed files: {n_transposed}\n\n"))


# Standardising column names for all transposed files

# List all transposed CSV files from transposed directory
transposed_files <- list.files(path = TRANSPOSED_CSV_DIR, pattern = "^transposed_.*\\.csv$", full.names = TRUE)

# Define the new column names
new_colnames <- c(
  "gene_accession_id",
  "gene_symbol",
  "mouse_strain",
  "mouse_life_stage",
  "parameter_id",
  "pvalue",
  "parameter_name",
  "analysis_id"
)

# Apply standardised names to each transposed file
for (i in seq_along(transposed_files)) {
  file <- transposed_files[i]
  
  # Read the transposed CSV
  # check.names = FALSE prevents R from modifying column names
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Replace column names with our standardised names
  colnames(df) <- new_colnames
  
  # Overwrite the file with standardised column names
  write.csv(df, file, row.names = FALSE)
  
  # Display progress
  cat(glue("\rUpdated: {i}/{length(transposed_files)}"))
  flush.console()
}

cat("\n\n✓ All column names standardised!\n\n")


# Combining all transposed files into one

# List all transposed CSV files from transposed directory
transposed_files <- list.files(path = TRANSPOSED_CSV_DIR, pattern = "^transposed_.*\\.csv$", full.names = TRUE)

# Combine all files into one data frame
# lapply() reads each file
# do.call(rbind, ...) stacks all rows together
# This works because all files now have identical column names
combined <- do.call(rbind, lapply(transposed_files, read.csv, stringsAsFactors = FALSE))


# Create output directory if needed
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)

# Save combined file
write.csv(combined, OUTPUT_FILE, row.names = FALSE)

if (file.exists(OUTPUT_FILE)) {
  file_size_mb <- file.info(OUTPUT_FILE)$size / (1024^2)
  cat(glue("✓ Saved to: {OUTPUT_FILE}\n"))
  cat(glue("✓ File size: {round(file_size_mb, 2)} MB\n"))
} else {
  cat("✗ ERROR: Failed to save file\n")
}