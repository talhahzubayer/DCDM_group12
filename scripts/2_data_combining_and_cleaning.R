# FILE DIRECTORY DECLARATIONS

# Folder containing Transposed CSVs 
transposed_folder <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/raw/phenotype_analysis_results_transposed"

# Output files
file_combined <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/data_table_combined.csv"
file_cleaned <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/data_table_cleaned.csv"
file_removed <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/data_table_removed.csv"
file_noNA <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/data_table_noNA.csv"
file_noNA_noDup_noOutliers <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/data_table_noNA_noDup_noOutliers.csv"
file_final <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/clean_table_final_UPPERCASE.csv"


# PART 1: COMBINE ALL TRANSPOSED CSV FILES

# List all transposed CSV files
file_list <- list.files(path = transposed_folder, pattern = "\\.csv$", full.names = TRUE)

# Read all files into a list
df_list <- lapply(file_list, function(file) {
  read.csv(file, stringsAsFactors = FALSE)
})

# Define the expected 8 columns
expected_cols <- c("gene_accession_id", "gene_symbol", "mouse_strain", 
                   "mouse_life_stage", "parameter_id", "pvalue", 
                   "parameter_name", "analysis_id")

# Align all data frames to have these exact columns
df_list_aligned <- lapply(df_list, function(df) {
  
  # Clean up column name spaces
  names(df) <- trimws(names(df))
  
  # Add missing columns as NA
  missing_cols <- setdiff(expected_cols, names(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  
  # Keep only the expected columns
  df <- df[expected_cols]
  
  return(df)
})

# Combine all aligned data frames into one table
combined_data <- do.call(rbind, df_list_aligned)

# Save the combined table
write.csv(combined_data, file_combined, row.names = FALSE)

cat("Combined table created with", nrow(combined_data), "rows and", ncol(combined_data), "columns.\n")


# PART 2: CLEAN THE COMBINED DATA

# Helper function to check column length
check_length <- function(column, min_len, max_len) {
  nchar_val <- nchar(as.character(column))
  return(nchar_val >= min_len & nchar_val <= max_len)
}

# Create a copy for cleaning
cleaned_data <- combined_data

# Initialize a logical vector to track valid rows
valid_rows <- rep(TRUE, nrow(cleaned_data))

# Apply SOP validation rules
valid_rows <- valid_rows & (nchar(as.character(cleaned_data$analysis_id)) == 15)
valid_rows <- valid_rows & check_length(cleaned_data$gene_accession_id, 9, 11)
valid_rows <- valid_rows & check_length(cleaned_data$gene_symbol, 1, 13)
valid_rows <- valid_rows & check_length(cleaned_data$mouse_strain, 3, 5)
valid_rows <- valid_rows & check_length(cleaned_data$mouse_life_stage, 4, 17)
valid_rows <- valid_rows & check_length(cleaned_data$parameter_id, 15, 20)
valid_rows <- valid_rows & check_length(cleaned_data$parameter_name, 2, 74)
cleaned_data$pvalue <- as.numeric(cleaned_data$pvalue)
valid_rows <- valid_rows & (cleaned_data$pvalue >= 0 & cleaned_data$pvalue <= 1)

# Split into cleaned and removed tables
final_cleaned <- cleaned_data[valid_rows, ]
removed_data <- cleaned_data[!valid_rows, ]

# Save the cleaned table and removed rows
write.csv(final_cleaned, file_cleaned, row.names = FALSE)
write.csv(removed_data, file_removed, row.names = FALSE)

cat("Cleaning complete.\n")
cat("Remaining rows (cleaned):", nrow(final_cleaned), "\n")
cat("Removed rows (failed SOP):", nrow(removed_data), "\n")


# PART 3: REMOVE NA VALUES

# Remove any rows containing NA
noNA_data <- na.omit(final_cleaned)

# Save the table without NAs
write.csv(noNA_data, file_noNA, row.names = FALSE)

cat("New table created without NAs.\n")
cat("Rows remaining:", nrow(noNA_data), "\n")


# PART 4: REMOVE DUPLICATES AND OUTLIERS

# Remove duplicate rows based on all columns
dedup_data <- noNA_data[!duplicated(noNA_data), ]

cat("Duplicates removed:", nrow(noNA_data) - nrow(dedup_data), "\n")

# Identify numeric columns
numeric_cols <- sapply(dedup_data, is.numeric)

# Function to flag outliers using IQR method
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  # Values within 1.5 * IQR from Q1 and Q3 are kept
  x >= (Q1 - 1.5 * IQR_val) & x <= (Q3 + 1.5 * IQR_val)
}

# Initialize a logical vector
keep_rows <- rep(TRUE, nrow(dedup_data))

# Apply outlier filter for all numeric columns
for (col in names(dedup_data)[numeric_cols]) {
  valid <- remove_outliers(dedup_data[[col]])
  keep_rows <- keep_rows & valid
}

# Filter data to remove outliers
cleaned_final <- dedup_data[keep_rows, ]

cat("Outliers removed:", nrow(dedup_data) - nrow(cleaned_final), "\n")

# Save the final cleaned table
write.csv(cleaned_final, file_noNA_noDup_noOutliers, row.names = FALSE)

cat("Further cleaning complete.\n")
cat("Rows remaining after removing duplicates and outliers:", nrow(cleaned_final), "\n")


# PART 5: CONVERT GENE_ACCESSION_ID TO UPPERCASE

# Make a copy
final_data <- cleaned_final

# Convert gene_accession_id column to uppercase
final_data$gene_accession_id <- toupper(final_data$gene_accession_id)

# Save the final table with uppercase gene IDs
write.csv(final_data, file_final, row.names = FALSE)

cat("Final table saved with uppercase gene_accession_id.\n")
