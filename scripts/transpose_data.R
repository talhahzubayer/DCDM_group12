# Load required libraries
library(tidyverse)  # For data manipulation functions
library(glue)       # For string formatting in messages

# Directory containing original CSV files (vertical format)
CSV_DIR <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/raw/phenotype_analysis_results"

# Directory to save transposed CSV files (horizontal format)
TRANSPOSED_CSV_DIR <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/raw/phenotype_analysis_results_transposed"

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