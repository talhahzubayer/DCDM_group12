library(dplyr)
library(readr)
library(stringr) # for str_trim

# Directory path of Where the IMPC_parameter_description.txt is stored
INPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/metadata/IMPC_parameter_description.txt"

# Final clean file directory path
OUTPUT_FILE <- "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_parameter_description_cleaned.csv"

# Read the file; blanks read as NA automatically
data <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Trim whitespace from 'impcParameterOrigId' and 'parameterId'
data <- data %>%
  mutate(
    impcParameterOrigId = str_trim(impcParameterOrigId, side = "both"),
    parameterId = str_trim(parameterId, side = "both")
  )

# Remove exact duplicate rows
data_clean <- distinct(data)

# Remove duplicates by impcParameterOrigId, treating NA as unique
data_clean <- data_clean %>%
  filter(!is.na(impcParameterOrigId)) %>%
  distinct(impcParameterOrigId, .keep_all = TRUE) %>%
  bind_rows(
    data_clean %>% filter(is.na(impcParameterOrigId))
  )

# Remove duplicates by parameterId, treating NA as unique
data_clean <- data_clean %>%
  filter(!is.na(parameterId)) %>%
  distinct(parameterId, .keep_all = TRUE) %>%
  bind_rows(
    data_clean %>% filter(is.na(parameterId))
  )

# Write out the cleaned data with blanks as "NA"
write_csv(data_clean, OUTPUT_FILE, na = "NA")