library(tidyverse)
library(RMySQL)


# CONFIGURATION

# Data source: "csv" or "database"
DATA_SOURCE <- "database"

# CSV file paths
CSV_FILES <- list(
  phenotype = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/clean_table_final_UPPERCASE.csv",
  disease = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/Disease_information_cleaned.csv",
  parameter = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_parameter_description_cleaned.csv",
  procedure = "C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_procedure_cleaned.csv"
)

# Database configuration
DB_CONFIG <- list(
  host = "localhost",
  port = 3306,
  dbname = "dcdm_cw1_group12_test",
  user = "root",
  password = "Talhah17!"
)


# HELPER FUNCTION: Parameter Categorisation (FOR CSV ONLY)

# When loading from database, parameter groups come from Parameter_Groups table
# This function is only used for CSV files without group information

categorise_parameter <- function(param_name, param_id) {
  param_lower <- tolower(paste(param_name, param_id, sep = " "))
  
  if (grepl("housing|cage|nutrition|diet|temperature|humidity|light|water|impc_hou", param_lower)) {
    return("Housing & Environment")
  }
  if (grepl("weight|severity score|skin|coat|bone|skeleton|muscle|organ|anatomical|structure|length|thickness|size|dimension|morphology|composition|impc_owt|impc_his|impc_dxa|impc_bwt|impc_eye", param_lower)) {
    return("Structural Phenotype")
  }
  if (grepl("blood|plasma|serum|hemoglobin|lymphocyte|neutrophil|erythrocyte|platelet|leukocyte|glucose|cholesterol|protein|hormone|immunoglobulin|assay|bilirubin|alt|enzyme|chemistry|metabolite|impc_cbc|impc_cld", param_lower)) {
    return("Clinical Chemistry/Blood")
  }
  if (grepl("embryo|placenta|fetus|development|viability|dead|stillborn|gestation|umbilic|yolk sac|defect|anomaly|impc_evo", param_lower)) {
    return("Embryo & Development")
  }
  if (grepl("limb|forelimb|hindlimb|grip|paw|digit|impc_grs|impc_rot", param_lower)) {
    return("Limb Function & Performance")
  }
  if (grepl("behavior|activity|locomotor|anxiety|sleep|startle|memory|neurological|brainwave|eeg|impc_oft|impc_aas|impc_ecs|impc_het", param_lower)) {
    return("Behavioral & Neurological")
  }
  if (grepl("ecg|heart rate|cardiac|electrocardiogram|blood pressure|impc_ecg", param_lower)) {
    return("Cardiovascular & ECG")
  }
  
  return("Other")
}


# DATA LOADING: CSV FILES

load_data_from_csv <- function(include_procedure = TRUE) {
  # Validate required files exist
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
  
  # Load main phenotype data
  df <- read.csv(CSV_FILES$phenotype, stringsAsFactors = FALSE)
  param_desc <- read.csv(CSV_FILES$parameter, stringsAsFactors = FALSE)
  
  # Load and link procedure information
  if (include_procedure) {
    procedure_df <- read.csv(CSV_FILES$procedure, stringsAsFactors = FALSE)
    
    param_desc_with_procedure <- param_desc %>%
      left_join(procedure_df, 
                by = "impcParameterOrigId", 
                suffix = c("_param", "_proc")) %>%
      select(parameterId, 
             parameter_description = description_param,
             procedure_name = name_proc,             
             procedure_description = description_proc,
             is_mandatory = isMandatory)
    
    df <- df %>%
      left_join(param_desc_with_procedure, by = c("parameter_id" = "parameterId"))
  }
  
  # Add empty procedure columns if not loaded
  if (!include_procedure || !"procedure_name" %in% colnames(df)) {
    df$procedure_name <- NA_character_
    df$procedure_description <- NA_character_
    df$is_mandatory <- NA
  }
  
  # Add parameter categories
  df$category <- mapply(categorise_parameter, df$parameter_name, df$parameter_id)
  
  # Create combined variables for UI
  df$gene_combined <- paste0(df$gene_symbol, " (", df$gene_accession_id, ")")
  df$parameter_combined <- paste0(df$parameter_name, " (", df$parameter_id, ")")
  
  # Handle zero p-values
  min_nonzero <- min(df$pvalue[df$pvalue > 0], na.rm = TRUE)
  if (any(df$pvalue == 0)) {
    df$pvalue[df$pvalue == 0] <- min_nonzero
  }
  
  # Add transformed p-value
  df$neg_log10_pvalue <- -log10(df$pvalue)
  
  # Convert to factors for performance
  df$gene_symbol <- as.factor(df$gene_symbol)
  df$parameter_name <- as.factor(df$parameter_name)
  df$category <- as.factor(df$category)
  df$mouse_strain <- as.factor(df$mouse_strain)
  df$mouse_life_stage <- as.factor(df$mouse_life_stage)
  
  return(df)
}


# DATA LOADING: MYSQL DATABASE

load_data_from_database <- function() {
  # Connect to database
  con <- tryCatch({
    dbConnect(MySQL(),
              host = DB_CONFIG$host,
              dbname = DB_CONFIG$dbname,
              user = DB_CONFIG$user,
              password = DB_CONFIG$password)
  }, error = function(e) {
    stop("Database connection failed: ", e$message)
  })
  
  # Query joins all required tables
  # Analysis contains core phenotype data
  # Links to Genes, Parameters, Procedures, Parameter_Groups, and Diseases
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
  
  df <- tryCatch({
    dbGetQuery(con, query)
  }, error = function(e) {
    dbDisconnect(con)
    stop("Query failed: ", e$message)
  })
  
  dbDisconnect(con)
  
  # CRITICAL: MySQL/RMySQL may return pvalue as character type
  # Must explicitly convert to numeric for mathematical operations
  df$pvalue <- as.numeric(as.character(df$pvalue))
  
  # Remove rows with invalid p-values
  df <- df[!is.na(df$pvalue) & df$pvalue >= 0 & df$pvalue <= 1, ]
  
  # Use parameter groups from database, fallback to keyword-based if missing
  if (!"category" %in% colnames(df) || all(is.na(df$category))) {
    df$category <- mapply(categorise_parameter, df$parameter_name, df$parameter_id)
  }
  
  # Handle NULL categories
  df$category[is.na(df$category)] <- "Other"
  
  # Create combined variables for UI
  df$gene_combined <- paste0(df$gene_symbol, " (", df$gene_accession_id, ")")
  df$parameter_combined <- paste0(df$parameter_name, " (", df$parameter_id, ")")
  
  # Handle zero p-values
  min_nonzero <- min(df$pvalue[df$pvalue > 0], na.rm = TRUE)
  if (any(df$pvalue == 0, na.rm = TRUE)) {
    df$pvalue[df$pvalue == 0] <- min_nonzero
  }
  
  # Add transformed p-value
  df$neg_log10_pvalue <- -log10(df$pvalue)
  
  # Convert to factors for performance
  df$gene_symbol <- as.factor(df$gene_symbol)
  df$parameter_name <- as.factor(df$parameter_name)
  df$category <- as.factor(df$category)
  df$mouse_strain <- as.factor(df$mouse_strain)
  df$mouse_life_stage <- as.factor(df$mouse_life_stage)
  
  return(df)
}


# MAIN DATA LOADING FUNCTION

# Switches between CSV and database sources based on configuration
# Returns standardized data frame for use in RShiny dashboard

load_data <- function(source = DATA_SOURCE, include_procedure = TRUE) {
  if (source == "csv") {
    return(load_data_from_csv(include_procedure = include_procedure))
  } else if (source == "database") {
    return(load_data_from_database())
  } else {
    stop("Invalid data source. Use 'csv' or 'database'")
  }
}


# PARAMETER GROUP COLOR SCHEME
# Colors used in dashboard visualisations

GROUP_COLORS <- c(
  'Housing & Environment' = '#deb887',
  'Structural Phenotype' = '#ff1493',
  'Clinical Chemistry/Blood' = '#ff0000',
  'Embryo & Development' = '#ff69b4',
  'Limb Function & Performance' = '#4169e1',
  'Behavioral & Neurological' = '#00bfff',
  'Cardiovascular & ECG' = '#8b0000',
  'Other' = '#696969'
)