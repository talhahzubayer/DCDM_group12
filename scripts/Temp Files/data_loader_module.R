library(tidyverse)
library(RMySQL)


# CONFIGURATION

# Data source: "csv" or "database"
DATA_SOURCE <- "csv"  # ← CHANGE THIS TO "database" WHEN READY!

# CSV file paths (update these to match your folder structure)
CSV_FILES <- list(
  phenotype = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/clean_table_final_UPPERCASE.csv",
  disease = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/Disease_information_cleaned.csv",
  parameter = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_parameter_description_cleaned.csv",
  procedure = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_procedure_cleaned.csv"
)

# Database configuration (fill in when database is ready)
DB_CONFIG <- list(
  host = "localhost",           # Update when database ready
  dbname = "impc_group12",      # Update when database ready
  user = "your_username",       # Update when database ready
  password = "your_password"    # Update when database ready
)


# HELPER FUNCTION: Parameter Categorisation

categorize_parameter <- function(param_name, param_id) {
  param_lower <- tolower(paste(param_name, param_id, sep = " "))
  
  if (grepl("weight|mass|bw|body weight", param_lower)) return("Weight")
  if (grepl("image|xray|picture", param_lower)) return("Images")
  if (grepl("brain|neuro|behavior|cognitive|memory|learning|startle|inhibition|locomotor|activity|arousal|vocalization|movement|speed|distance|rearing|mobile|holepoke", param_lower)) return("Brain")
  if (grepl("eye|retina|lens|vision|optic", param_lower)) return("Vision/Eye")
  if (grepl("blood|immune|cell|hemoglobin|platelet|lymph|plasma", param_lower)) return("Blood")
  if (grepl("heart|cardio|pressure|pulse|vessel", param_lower)) return("Cardiovascular")
  if (grepl("muscle|strength|motor", param_lower)) return("Muscular")
  if (grepl("metabolism|glucose|lipid|insulin|cholesterol|fat", param_lower)) return("Metabolic")
  if (grepl("lung|respirat|breath|hypoxia", param_lower)) return("Respiratory")
  if (grepl("placenta|vagina|ovary|teste|uterus|prostate|penis|umbilic", param_lower)) return("Reproductive")
  if (grepl("skin|coat|hair|dermis", param_lower)) return("Coat/Skin")
  if (grepl("housing|cage", param_lower)) return("Housing")
  if (grepl("condition|context|cue", param_lower)) return("Conditions")
  if (grepl("cbc|biochemical", param_lower)) return("Biochemical")
  if (grepl("equipment", param_lower)) return("Equipment")
  
  return("Other")
}


# DATA LOADING FUNCTION: LOAD FROM CSV FILES

load_data_from_csv <- function(include_procedure = TRUE) {
  cat("Loading data from CSV files...\n")
  
  # Check if all required files exist
  required_files <- c("phenotype", "parameter")
  if (include_procedure) required_files <- c(required_files, "procedure")
  
  missing_files <- c()
  for (name in required_files) {
    if (!file.exists(CSV_FILES[[name]])) {
      missing_files <- c(missing_files, CSV_FILES[[name]])
    }
  }
  
  if (length(missing_files) > 0) {
    stop("Missing CSV files:\n", paste(missing_files, collapse = "\n"))
  }
  
  # Load phenotype analysis data (main data)
  cat("  Loading phenotype analysis data...\n")
  df <- read.csv(CSV_FILES$phenotype, stringsAsFactors = FALSE)
  cat(glue::glue("  Loaded {nrow(df)} phenotype records\n"))
  
  # Load parameter description file (for linking)
  cat("  Loading parameter descriptions...\n")
  param_desc <- read.csv(CSV_FILES$parameter, stringsAsFactors = FALSE)
  cat(glue::glue("  Loaded {nrow(param_desc)} parameter descriptions\n"))
  
  # Load procedure information if requested
  if (include_procedure) {
    cat("  Loading procedure information...\n")
    procedure_df <- read.csv(CSV_FILES$procedure, stringsAsFactors = FALSE)
    cat(glue::glue("  Loaded {nrow(procedure_df)} procedure records\n"))
    
    # Link procedure data to parameter descriptions via impcParameterOrigId
    # CRITICAL: After left_join with suffix, columns are renamed:
    #   - name becomes name_param and name_proc
    #   - description becomes description_param and description_proc
    #   - isMandatory stays as isMandatory (no duplicate)
    param_desc_with_procedure <- param_desc %>%
      left_join(procedure_df, by = "impcParameterOrigId", suffix = c("_param", "_proc")) %>%
      select(parameterId, 
             parameter_description = description_param,
             procedure_name = name_proc,             
             procedure_description = description_proc,
             is_mandatory = isMandatory)
    
    # Now join with main phenotype data
    df <- df %>%
      left_join(param_desc_with_procedure, by = c("parameter_id" = "parameterId"))
    
    cat(glue::glue("  Linked procedure info to {sum(!is.na(df$procedure_description))} parameters\n"))
  }
  
  # If procedure info was not loaded, add empty columns for compatibility
  if (!include_procedure || !"procedure_name" %in% colnames(df)) {
    cat("  Adding empty procedure columns for compatibility\n")
    df$procedure_name <- NA_character_
    df$procedure_description <- NA_character_
    df$is_mandatory <- NA
  }
  
  # Add parameter categories
  cat("  Adding parameter categories...\n")
  df$category <- mapply(categorize_parameter, df$parameter_name, df$parameter_id)
  
  # Create combined variables for UI
  df$gene_combined <- paste0(df$gene_symbol, " (", df$gene_accession_id, ")")
  df$parameter_combined <- paste0(df$parameter_name, " (", df$parameter_id, ")")
  
  # Handle zero p-values
  min_nonzero <- min(df$pvalue[df$pvalue > 0], na.rm = TRUE)
  if (any(df$pvalue == 0)) {
    cat("  Replacing zero p-values with minimum non-zero value\n")
    df$pvalue[df$pvalue == 0] <- min_nonzero
  }
  
  # Add transformed p-value
  df$neg_log10_pvalue <- -log10(df$pvalue)
  
  # Convert to factors for better performance
  df$gene_symbol <- as.factor(df$gene_symbol)
  df$parameter_name <- as.factor(df$parameter_name)
  df$category <- as.factor(df$category)
  df$mouse_strain <- as.factor(df$mouse_strain)
  df$mouse_life_stage <- as.factor(df$mouse_life_stage)
  
  cat("  Data loading complete!\n\n")
  cat(glue::glue("Summary:\n"))
  cat(glue::glue("  Total records: {nrow(df)}\n"))
  cat(glue::glue("  Unique genes: {n_distinct(df$gene_symbol)}\n"))
  cat(glue::glue("  Unique parameters: {n_distinct(df$parameter_id)}\n"))
  cat(glue::glue("  Mouse strains: {n_distinct(df$mouse_strain)}\n"))
  cat(glue::glue("  Life stages: {n_distinct(df$mouse_life_stage)}\n"))
  if ("procedure_name" %in% colnames(df)) {
    cat(glue::glue("  Parameters with procedure info: {sum(!is.na(df$procedure_description))}\n"))
  }
  cat("\n")
  
  return(df)
}


# DATA LOADING FUNCTION: LOAD FROM DATABASE

load_data_from_database <- function() {
  cat("Loading data from MySQL database...\n")
  
  # Connect to database
  con <- tryCatch({
    dbConnect(MySQL(),
              host = DB_CONFIG$host,
              dbname = DB_CONFIG$dbname,
              user = DB_CONFIG$user,
              password = DB_CONFIG$password)
  }, error = function(e) {
    stop("Database connection failed: ", e$message, 
         "\n\nMake sure to update DB_CONFIG at the top of this file!")
  })
  
  cat("  Connected to database\n")
  
  # Enhanced query with procedure information
  query <- "
    SELECT 
      a.analysis_id,
      a.pvalue,
      g.gene_accession_id,
      g.gene_symbol,
      p.parameter_id,
      p.parameter_name,
      p.parameter_group as category,
      p.parameter_description,
      proc.procedure_name,
      proc.procedure_description,
      proc.is_mandatory,
      a.mouse_strain,
      a.mouse_life_stage
    FROM Analysis a
    INNER JOIN Gene g ON g.gene_accession_id = a.gene_accession_id
    INNER JOIN Parameter p ON p.parameter_id = a.parameter_id
    LEFT JOIN Procedure proc ON p.impc_parameter_orig_id = proc.impc_parameter_orig_id
  "
  
  cat("  Executing query...\n")
  df <- tryCatch({
    dbGetQuery(con, query)
  }, error = function(e) {
    dbDisconnect(con)
    stop("Query failed: ", e$message)
  })
  
  dbDisconnect(con)
  cat("  Database connection closed\n")
  
  cat(glue::glue("  Loaded {nrow(df)} phenotype records\n"))
  
  # Add categories if not present or NULL
  if (!"category" %in% colnames(df) || all(is.na(df$category))) {
    cat("  Adding parameter categories...\n")
    df$category <- mapply(categorize_parameter, df$parameter_name, df$parameter_id)
  }
  
  # Create combined variables for UI
  df$gene_combined <- paste0(df$gene_symbol, " (", df$gene_accession_id, ")")
  df$parameter_combined <- paste0(df$parameter_name, " (", df$parameter_id, ")")
  
  # Handle zero p-values
  min_nonzero <- min(df$pvalue[df$pvalue > 0], na.rm = TRUE)
  if (any(df$pvalue == 0)) {
    cat("  Replacing zero p-values with minimum non-zero value\n")
    df$pvalue[df$pvalue == 0] <- min_nonzero
  }
  
  # Add transformed p-value
  df$neg_log10_pvalue <- -log10(df$pvalue)
  
  # Convert to factors
  df$gene_symbol <- as.factor(df$gene_symbol)
  df$parameter_name <- as.factor(df$parameter_name)
  df$category <- as.factor(df$category)
  df$mouse_strain <- as.factor(df$mouse_strain)
  df$mouse_life_stage <- as.factor(df$mouse_life_stage)
  
  cat("  Data loading complete!\n\n")
  cat(glue::glue("Summary:\n"))
  cat(glue::glue("  Total records: {nrow(df)}\n"))
  cat(glue::glue("  Unique genes: {n_distinct(df$gene_symbol)}\n"))
  cat(glue::glue("  Unique parameters: {n_distinct(df$parameter_id)}\n"))
  cat(glue::glue("  Mouse strains: {n_distinct(df$mouse_strain)}\n"))
  cat(glue::glue("  Life stages: {n_distinct(df$mouse_life_stage)}\n"))
  cat(glue::glue("  Parameters with procedure info: {sum(!is.na(df$procedure_description))}\n\n"))
  
  return(df)
}


# MAIN DATA LOADING FUNCTION (SWITCH BETWEEN SOURCES)

load_data <- function(source = DATA_SOURCE, include_procedure = TRUE) {
  cat("============================================\n")
  cat("IMPC DATA LOADER (ENHANCED)\n")
  cat("============================================\n")
  cat(glue::glue("Data source: {toupper(source)}\n"))
  cat(glue::glue("Include procedure info: {include_procedure}\n\n"))
  
  if (source == "csv") {
    return(load_data_from_csv(include_procedure = include_procedure))
  } else if (source == "database") {
    return(load_data_from_database())
  } else {
    stop("Invalid data source. Use 'csv' or 'database'")
  }
}


# USAGE EXAMPLES
# Example 1: Load from CSV with procedure info (default)
# data <- load_data()

# Example 2: Load from CSV without procedure info (faster)
# data <- load_data(include_procedure = FALSE)

# Example 3: Load from database (when ready)
# data <- load_data(source = "database")

cat("============================================\n")
cat("DATA LOADER MODULE LOADED\n")
cat("============================================\n")
cat("Current data source:", toupper(DATA_SOURCE), "\n")
cat("Procedure information: ENABLED\n\n")
cat("To change data source:\n")
cat("  1. Update DATA_SOURCE variable at top of file\n")
cat("  2. Or call: load_data(source = 'csv') or load_data(source = 'database')\n\n")
cat("To disable procedure info:\n")
cat("  load_data(include_procedure = FALSE)\n\n")
cat("To use in RShiny:\n")
cat("  data <- reactive({ load_data() })\n\n")