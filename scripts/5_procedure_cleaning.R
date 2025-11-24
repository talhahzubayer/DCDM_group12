# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)  # For replace_na function

# Input and output file paths
INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/metadata/IMPC_procedure.txt"
INTERMEDIATE_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_procedure.csv"
OUTPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_procedure_cleaned2.csv"

# Read the comma-separated text file
data_table <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

# View the first few rows
head(data_table)

# Save it as a CSV file
write.csv(data_table, INTERMEDIATE_FILE, row.names = FALSE)

# Read CSV
df <- read_csv(INTERMEDIATE_FILE)

# Reorder columns: move impcParameterOrigId first
df <- df %>% select(impcParameterOrigId, everything())

# Check for duplicates
duplicate_check <- df %>%
  count(impcParameterOrigId) %>%
  filter(n > 1)

cat("=== DUPLICATE CHECK ===\n")
cat("Total rows:", nrow(df), "\n")
cat("Unique origin IDs:", n_distinct(df$impcParameterOrigId), "\n")
cat("Duplicate ID groups:", nrow(duplicate_check), "\n")

# Remove duplicates if any
df_clean <- if (nrow(duplicate_check) > 0) {
  cat("\nDuplicate IDs found. Removing duplicates...\n")
  df %>% distinct(impcParameterOrigId, .keep_all = TRUE)
} else {
  cat("\nNo duplicates found.\n")
  df
}

# Clean data
df_clean <- df_clean %>%
  mutate(
    description = description %>%
      str_replace_all("&nbsp;", " ") %>%
      str_trim() %>%
      str_squish() %>%
      replace_na("No description available"),
    name = str_trim(name),
    isMandatory = str_trim(isMandatory) %>% toupper()
  )

# Final summary
cat("\n=== FINAL CLEANED DATA ===\n")
cat("Total records:", nrow(df_clean), "\n")
cat("Unique procedures:", n_distinct(df_clean$name), "\n")
cat("Mandatory procedures:", sum(df_clean$isMandatory == "TRUE"), "\n")
cat("Optional procedures:", sum(df_clean$isMandatory == "FALSE"), "\n")

# Sample data
cat("\nSample of cleaned data:\n")
print(head(df_clean, 10))

# Save cleaned data to CSV
write_csv(df_clean, OUTPUT_FILE)
cat("\nCleaned data saved to:", OUTPUT_FILE, "\n")