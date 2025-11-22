library(tidyverse)
library(RMySQL)


# ============================================================================
# CONFIGURATION
# ============================================================================

# Data source: "csv" or "database"
DATA_SOURCE <- "csv"  # ← CHANGE THIS TO "database" WHEN READY!

# CSV file paths (update these to match your folder structure)
CSV_FILES <- list(
  phenotype = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/clean_table_final_UPPERCASE.csv",
  disease = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/Disease_information_cleaned.csv",
  parameter = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_parameter_description_cleaned.csv",
  procedure = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_procedure_cleaned.csv"
)

# Database configuration
# UPDATE THESE VALUES to match your MySQL setup
DB_CONFIG <- list(
  host = "localhost",           # Usually "localhost" for local MySQL
  dbname = "DCDM_CW1_GROUP12",  # Your database name from SQL script
  user = "your_username",       # Your MySQL username
  password = "your_password"    # Your MySQL password
)


# ============================================================================
# HELPER FUNCTION: Parameter Categorisation (FOR CSV LOADING ONLY)
# ============================================================================
# NOTE: When loading from database, parameter groups come from Parameter_Groups table
# This function is only used for CSV files that don't have group information

categorise_parameter <- function(param_name, param_id) {
  param_lower <- tolower(paste(param_name, param_id, sep = " "))
  
  # Housing & Environment
  if (grepl("housing|cage|nutrition|diet|temperature|humidity|light|water|impc_hou", param_lower)) {
    return("Housing & Environment")
  }
  
  # Structural Phenotype
  if (grepl("weight|severity score|skin|coat|bone|skeleton|muscle|organ|anatomical|structure|length|thickness|size|dimension|morphology|composition|impc_owt|impc_his|impc_dxa|impc_bwt|impc_eye", param_lower)) {
    return("Structural Phenotype")
  }
  
  # Clinical Chemistry/Blood
  if (grepl("blood|plasma|serum|hemoglobin|lymphocyte|neutrophil|erythrocyte|platelet|leukocyte|glucose|cholesterol|protein|hormone|immunoglobulin|assay|bilirubin|alt|enzyme|chemistry|metabolite|impc_cbc|impc_cld", param_lower)) {
    return("Clinical Chemistry/Blood")
  }
  
  # Embryo & Development
  if (grepl("embryo|placenta|fetus|development|viability|dead|stillborn|gestation|umbilic|yolk sac|defect|anomaly|impc_evo", param_lower)) {
    return("Embryo & Development")
  }
  
  # Limb Function & Performance
  if (grepl("limb|forelimb|hindlimb|grip|paw|digit|impc_grs|impc_rot", param_lower)) {
    return("Limb Function & Performance")
  }
  
  # Behavioral & Neurological
  if (grepl("behavior|activity|locomotor|anxiety|sleep|startle|memory|neurological|brainwave|eeg|impc_oft|impc_aas|impc_ecs|impc_het", param_lower)) {
    return("Behavioral & Neurological")
  }
  
  # Cardiovascular & ECG
  if (grepl("ecg|heart rate|cardiac|electrocardiogram|blood pressure|impc_ecg", param_lower)) {
    return("Cardiovascular & ECG")
  }
  
  # Default to Other
  return("Other")
}


# ============================================================================
# DATA LOADING FUNCTION: LOAD FROM CSV FILES
# ============================================================================

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
    
    # Link procedure data to parameter descriptions
    # Match column names from your SQL schema
    param_desc_with_procedure <- param_desc %>%
      left_join(procedure_df, 
                by = "impcParameterOrigId", 
                suffix = c("_param", "_proc")) %>%
      select(parameterId, 
             parameter_description = description_param,
             procedure_name = name_proc,             
             procedure_description = description_proc,
             is_mandatory = isMandatory)
    
    # Join with main phenotype data
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
  
  # Add parameter categories using NEW GROUP NAMES
  cat("  Adding parameter categories...\n")
  df$category <- mapply(categorise_parameter, df$parameter_name, df$parameter_id)
  
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


# ============================================================================
# DATA LOADING FUNCTION: LOAD FROM DATABASE
# ============================================================================
# UPDATED to match your actual database schema from SQL scripts

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
         "\n\nMake sure to update DB_CONFIG at the top of this file!",
         "\nCheck: host, dbname, username, password")
  })
  
  cat("  Connected to database successfully\n")
  
  # CORRECTED QUERY based on actual database schema (v3.0)
  # Table names: Genes, Analysis, Parameters, Parameter_Groups, Procedures, Diseases
  # FIXED: Removed non-existent MouseCharacteristics table
  # FIXED: parameter_name comes from Parameters (P), not Genes (G)
  # FIXED: mouse_strain and mouse_life_stage are directly in Analysis (A)
  # ADDED: Disease associations for better dashboard integration
  query <- "
    SELECT 
      A.analysis_key_id AS analysis_id,
      A.pvalue,
      A.analysis_id AS original_analysis_id,
      G.gene_accession_id,
      G.gene_symbol,
      P.IMPC_parameter_id AS parameter_id,
      P.parameter_name,
      PG.group_name AS category,
      P.parameter_description,
      PROC.procedure_name,
      PROC.procedure_description,
      PROC.mandatory AS is_mandatory,
      A.mouse_strain,
      A.mouse_life_stage,
      D.disease_name,
      D.DO_ID,
      D.omim_ids
    FROM Analysis A
    INNER JOIN Genes G ON G.gene_id = A.gene_id
    INNER JOIN Parameters P ON P.parameter_id = A.parameter_id
    LEFT JOIN Parameter_Groups PG ON P.group_id = PG.group_id
    LEFT JOIN Procedures PROC ON P.IMPC_procedure_id = PROC.IMPC_procedure_id
    LEFT JOIN Diseases D ON G.gene_accession_id = D.gene_accession_id
  "
  
  cat("  Executing query...\n")
  df <- tryCatch({
    dbGetQuery(con, query)
  }, error = function(e) {
    dbDisconnect(con)
    stop("Query failed: ", e$message,
         "\nThis might be due to:",
         "\n  1. Table names not matching (check: Genes, Analysis, Parameters, etc.)",
         "\n  2. Column names not matching",
         "\n  3. Foreign key relationships not set up",
         "\n  4. Data not yet imported into tables")
  })
  
  dbDisconnect(con)
  cat("  Database connection closed\n")
  
  cat(glue::glue("  Loaded {nrow(df)} phenotype records\n"))
  
  # Add categories if not present or NULL
  # (Fallback to keyword-based categorization if database groups missing)
  if (!"category" %in% colnames(df) || all(is.na(df$category))) {
    cat("  Warning: No parameter groups in database, using keyword-based categorization...\n")
    df$category <- mapply(categorise_parameter, df$parameter_name, df$parameter_id)
  } else {
    cat("  Using parameter groups from database\n")
  }
  
  # Handle cases where category is still NULL after join
  df$category[is.na(df$category)] <- "Other"
  
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
  cat(glue::glue("  Parameters with procedure info: {sum(!is.na(df$procedure_description))}\n"))
  cat(glue::glue("  Parameter groups:\n"))
  cat(paste0("    - ", names(table(df$category)), ": ", table(df$category), "\n"))
  cat("\n")
  
  return(df)
}


# ============================================================================
# MAIN DATA LOADING FUNCTION (SWITCH BETWEEN SOURCES)
# ============================================================================

load_data <- function(source = DATA_SOURCE, include_procedure = TRUE) {
  cat("============================================\n")
  cat("IMPC DATA LOADER (GROUP 12 - UPDATED)\n")
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


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Example 1: Load from CSV with procedure info (default)
# data <- load_data()

# Example 2: Load from CSV without procedure info (faster)
# data <- load_data(include_procedure = FALSE)

# Example 3: Load from database (when ready)
# data <- load_data(source = "database")

# Example 4: Override default source
# data <- load_data(source = "database")

# Example 5: In RShiny dashboard
# data <- reactive({ load_data() })


# ============================================================================
# PARAMETER GROUP COLOR SCHEME (FOR DASHBOARD USE)
# ============================================================================
# Updated to match new parameter groups from your SQL

GROUP_COLORS <- c(
  'Housing & Environment' = '#deb887',      # Burlywood
  'Structural Phenotype' = '#ff1493',       # Hot pink
  'Clinical Chemistry/Blood' = '#ff0000',   # Red
  'Embryo & Development' = '#ff69b4',       # Pink
  'Limb Function & Performance' = '#4169e1',# Royal blue
  'Behavioral & Neurological' = '#00bfff',  # Deep sky blue
  'Cardiovascular & ECG' = '#8b0000',       # Dark red
  'Other' = '#696969'                       # Dim gray
)


# ============================================================================
# MODULE INITIALIZATION MESSAGE
# ============================================================================

cat("============================================\n")
cat("DATA LOADER MODULE LOADED (v3.0 - CORRECTED)\n")
cat("============================================\n")
cat("Current data source:", toupper(DATA_SOURCE), "\n")
cat("Database name:", DB_CONFIG$dbname, "\n")
cat("Procedure information: ENABLED\n\n")

cat("UPDATES IN THIS VERSION (v3.0):\n")
cat("  ✓ CRITICAL FIX: Removed non-existent MouseCharacteristics table\n")
cat("  ✓ CRITICAL FIX: parameter_name now from Parameters (P), not Genes (G)\n")
cat("  ✓ CRITICAL FIX: mouse_strain and mouse_life_stage from Analysis (A)\n")
cat("  ✓ OPTIMIZATION: Changed Procedures join to use direct FK\n")
cat("  ✓ ENHANCEMENT: Added Disease associations to query\n")
cat("  ✓ Updated parameter groups to match SQL schema:\n")
cat("    - Housing & Environment\n")
cat("    - Structural Phenotype\n")
cat("    - Clinical Chemistry/Blood\n")
cat("    - Embryo & Development\n")
cat("    - Limb Function & Performance\n")
cat("    - Behavioral & Neurological\n")
cat("    - Cardiovascular & ECG\n")
cat("    - Other\n\n")

cat("TO SWITCH TO DATABASE:\n")
cat("  1. Update DB_CONFIG at top of file with your credentials\n")
cat("  2. Set DATA_SOURCE <- 'database'\n")
cat("  3. Or call: load_data(source = 'database')\n\n")

cat("TO DISABLE PROCEDURE INFO (CSV only):\n")
cat("  load_data(include_procedure = FALSE)\n\n")

cat("TO USE IN RSHINY:\n")
cat("  source('data_loader_module.R')\n")
cat("  data <- reactive({ load_data() })\n\n")

cat("TROUBLESHOOTING:\n")
cat("  - If database connection fails, check DB_CONFIG\n")
cat("  - If query fails, verify table names in MySQL\n")
cat("  - Check that foreign keys are set up\n")
cat("  - Ensure data has been imported to all tables\n\n")