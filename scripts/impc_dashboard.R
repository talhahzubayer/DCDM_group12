library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)
library(umap)
library(metap)
library(pheatmap)
library(viridis)

# DATA LOADING
# Make sure the "data_loader_module.R" is in the same directory as the impc_dashboard script
source("data_loader_module.R")

# CONFIGURATION

QUERY_GENES <- c("Smarcd3", "Ppp3cc", "Rab12", "Klhl33")
DEFAULT_THRESHOLD <- 0.05

# Parameter group colors matching database Parameter_Groups table
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


# HELPER FUNCTIONS

# Combine multiple p-values using Fisher's method (sumlog)
# Used when multiple tests exist for same gene-parameter combination
combine_pvalue <- function(pvals) {
  pvals[pvals == 0] <- 1e-10  # Replace exact zeros to avoid log(0)
  if (length(pvals) < 2) return(pvals[1])
  else return(sumlog(pvals)$p)
}


# DASHBOARD USER INTERFACE

ui <- dashboardPage(
  dashboardHeader(title = "IMPC Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("home")),
      menuItem("Gene Analysis", tabName = "gene", icon = icon("dna")),
      menuItem("Phenotype Analysis", tabName = "phenotype", icon = icon("chart-bar")),
      menuItem("Gene Clustering", tabName = "clustering", icon = icon("project-diagram")),
      menuItem("Four Query Genes", tabName = "query_genes", icon = icon("star"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # OVERVIEW TAB
      tabItem(
        tabName = "overview",
        h2("Overview: Parameter Space Reduction"),
        p(em("Parameter groups based on IMPC biological classifications")),
        
        fluidRow(
          valueBoxOutput("total_genes"),
          valueBoxOutput("total_parameters"),
          valueBoxOutput("total_significant")
        ),
        
        fluidRow(
          box(width = 6, title = "Parameters per Group",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("overview_params")),
          box(width = 6, title = "Significant Associations per Group",
              status = "success", solidHeader = TRUE,
              plotlyOutput("overview_sig"))
        ),
        
        fluidRow(
          box(width = 12, title = "Parameter Groups Summary",
              status = "info", solidHeader = TRUE,
              p(strong("Note:"), "Groups correspond to major biological systems and experimental domains."),
              DTOutput("overview_table"))
        )
      ),
      # GENE ANALYSIS
      tabItem(
        tabName = "gene",
        h2("Gene-Centric Visualisation"),
        p(strong("Purpose:"), " Select a knockout mouse and visualise statistical scores of all phenotypes tested. 
          Shows which phenotypes are significantly affected by the gene knockout."),
        
        fluidRow(
          column(3,
                 wellPanel(
                   h4("Gene Selection"),
                   selectInput("gene_select", "Select Gene:", choices = NULL),
                   textAreaInput("gene_custom", "Or paste gene list (one per line):",
                                 value = paste(QUERY_GENES, collapse = "\n"),
                                 rows = 4, placeholder = "Gene symbols"),
                   hr(),
                   h4("Filters"),
                   selectInput("gene_category", "Parameter Group:", choices = NULL),
                   selectInput("gene_strain", "Strain:", choices = NULL),
                   selectInput("gene_lifestage", "Life Stage:", choices = NULL),
                   hr(),
                   sliderInput("gene_pvalue", "P-value Threshold:",
                               min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                   checkboxInput("gene_fdr", "Use FDR correction", value = FALSE)
                 )
          ),
          
          column(9,
                 box(width = 12, title = "Phenotype Significance Plot",
                     status = "primary", solidHeader = TRUE,
                     plotlyOutput("gene_plot", height = "500px")),
                 
                 box(width = 12, title = "Significant Phenotypes Table",
                     status = "info", solidHeader = TRUE,
                     p(strong("Columns show:"), " Parameter name, Parameter group, P-value, Disease link, Procedure"),
                     DTOutput("gene_table"))
          )
        )
      ),
      # PHENOTYPE ANALYSIS
      tabItem(
        tabName = "phenotype",
        h2("Phenotype-Centric Visualisation"),
        p(strong("Purpose:"), " Select a phenotype and visualise statistical scores of all knockout mice. 
          Shows which genes are significantly associated with that phenotype."),
        
        fluidRow(
          column(3,
                 wellPanel(
                   h4("Phenotype Selection"),
                   selectInput("pheno_category", "Parameter Group:", choices = NULL),
                   selectInput("pheno_parameter", "Select Parameter:", choices = NULL),
                   hr(),
                   h4("Parameter Details"),
                   uiOutput("pheno_details"),
                   hr(),
                   sliderInput("pheno_pvalue", "P-value Threshold:",
                               min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                   checkboxInput("pheno_fdr", "Use FDR correction", value = FALSE),
                   checkboxInput("pheno_highlight_disease", "Highlight disease-linked genes", value = TRUE)
                 )
          ),
          
          column(9,
                 box(width = 12, title = "Gene Significance Plot",
                     status = "primary", solidHeader = TRUE,
                     plotlyOutput("pheno_plot", height = "500px")),
                 
                 box(width = 12, title = "Significant Genes Table",
                     status = "info", solidHeader = TRUE,
                     p(strong("Columns show:"), " Gene symbol, MGI ID, P-value, Number of diseases linked, Top disease"),
                     DTOutput("pheno_table"))
          )
        )
      ),
      # CLUSTERING
      tabItem(
        tabName = "clustering",
        h2("Gene Clustering Based on Phenotype Profiles"),
        p(strong("Purpose:"), " Identify clusters of genes with similar phenotype scores. 
          Genes in the same cluster have similar biological effects."),
        
        fluidRow(
          column(3,
                 wellPanel(
                   h4("Clustering Settings"),
                   sliderInput("cluster_k", "Number of Clusters:",
                               min = 2, max = 10, value = 4, step = 1),
                   hr(),
                   h4("Filters"),
                   selectInput("cluster_category", "Parameter Group:", choices = NULL),
                   selectInput("cluster_strain", "Strain:", choices = NULL),
                   selectInput("cluster_lifestage", "Life Stage:", choices = NULL),
                   hr(),
                   sliderInput("cluster_pvalue", "P-value Threshold:",
                               min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                   checkboxInput("cluster_fdr", "Use FDR correction", value = FALSE),
                   hr(),
                   checkboxInput("cluster_highlight_query", "Highlight query genes", value = TRUE)
                 )
          ),
          
          column(9,
                 box(width = 12, title = "UMAP: Gene Clustering",
                     status = "primary", solidHeader = TRUE,
                     plotlyOutput("cluster_umap", height = "500px")),
                 
                 box(width = 12, title = "Correlation Heatmap",
                     status = "warning", solidHeader = TRUE,
                     plotOutput("cluster_heatmap", height = "500px")),
                 
                 box(width = 12, title = "Cluster Assignments",
                     status = "info", solidHeader = TRUE,
                     DTOutput("cluster_table"))
          )
        )
      ),
      # FOUR QUERY GENES
      tabItem(
        tabName = "query_genes",
        h2("The Four Genotypes of Interest"),
        p(strong("Purpose:"), " Summary of significant phenotypes for the four query genes: 
          Smarcd3, Ppp3cc, Rab12, Klhl33. Shows parameter groups, p-values, disease links, and procedures."),
        
        fluidRow(
          column(12,
                 wellPanel(
                   sliderInput("query_pvalue", "P-value Threshold:",
                               min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                   checkboxInput("query_fdr", "Use FDR correction", value = FALSE)
                 )
          )
        ),
        
        fluidRow(
          valueBoxOutput("query_smarcd3"),
          valueBoxOutput("query_ppp3cc"),
          valueBoxOutput("query_rab12"),
          valueBoxOutput("query_klhl33")
        ),
        
        fluidRow(
          box(width = 12, title = "Comparison of Query Genes",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("query_comparison", height = "400px"))
        ),
        
        fluidRow(
          box(width = 12, title = "Detailed Results for All Four Genes",
              status = "info", solidHeader = TRUE,
              p(strong("Shows:"), " Gene, Parameter group, Parameter name, P-value, Disease association, Procedure"),
              DTOutput("query_table"))
        )
      )
    )
  )
)


# SERVER LOGIC

server <- function(input, output, session) {
  
  # DATA LOADING & PREPARATION
  
  # Load main phenotype data and add FDR correction
  data <- reactive({
    df <- load_data()
    
    df <- df %>%
      mutate(p_adj = p.adjust(pvalue, method = "BH"),
             is_query = gene_symbol %in% QUERY_GENES)
    
    return(df)
  })
  
  # Load disease information with robust column name normalisation
  # Handles various naming conventions from different data sources
  disease_data <- reactive({
    tryCatch({
      disease <- read.csv(CSV_FILES$disease, header = TRUE, stringsAsFactors = FALSE)
      
      # Normalise gene_accession_id column name
      if (!"gene_accession_id" %in% names(disease)) {
        if ("GENE_ACCESSION_ID" %in% names(disease)) {
          disease <- dplyr::rename(disease, gene_accession_id = GENE_ACCESSION_ID)
        }
      }
      
      # Normalise disease term column name (multiple possible variants)
      if (!"disease_term" %in% names(disease)) {
        if ("DISEASE_TERM" %in% names(disease)) {
          disease <- dplyr::rename(disease, disease_term = DISEASE_TERM)
        } else if ("disease_name" %in% names(disease)) {
          disease <- dplyr::rename(disease, disease_term = disease_name)
        } else if ("DISEASE_LABEL" %in% names(disease)) {
          disease <- dplyr::rename(disease, disease_term = DISEASE_LABEL)
        } else {
          disease$disease_term <- NA_character_
        }
      }
      
      disease
    }, error = function(e) {
      data.frame(
        gene_accession_id = character(),
        disease_id        = character(),
        disease_term      = character(),
        stringsAsFactors  = FALSE
      )
    })
  })
  
  # Summarise disease associations per gene
  gene_disease_summary <- reactive({
    disease <- disease_data()
    if (nrow(disease) == 0 || !"gene_accession_id" %in% names(disease)) {
      return(data.frame(
        gene_accession_id = character(),
        n_diseases        = integer(),
        top_disease       = character(),
        stringsAsFactors  = FALSE
      ))
    }
    
    disease %>%
      group_by(gene_accession_id) %>%
      summarise(
        n_diseases = sum(!is.na(disease_term)),
        top_disease = {
          dt <- disease_term[!is.na(disease_term)]
          if (length(dt) == 0) NA_character_ else dt[1]
        },
        .groups = "drop"
      )
  })
  
  # Initialise dropdown menus with actual data from loaded dataset
  observe({
    df <- data()
    
    # Use actual categories from loaded data (works with both CSV and database)
    actual_categories <- sort(unique(as.character(df$category)))
    cats <- c("All" = "all", setNames(actual_categories, actual_categories))
    
    updateSelectInput(session, "gene_category", choices = cats)
    updateSelectInput(session, "pheno_category", choices = cats)
    updateSelectInput(session, "cluster_category", choices = cats)
    
    strains <- c("All" = "all", sort(unique(as.character(df$mouse_strain))))
    updateSelectInput(session, "gene_strain", choices = strains)
    updateSelectInput(session, "cluster_strain", choices = strains)
    
    lifestages <- c("All" = "all", sort(unique(as.character(df$mouse_life_stage))))
    updateSelectInput(session, "gene_lifestage", choices = lifestages)
    updateSelectInput(session, "cluster_lifestage", choices = lifestages)
    
    genes <- c("All" = "all", sort(unique(df$gene_combined)))
    updateSelectInput(session, "gene_select", choices = genes)
  })
  
  # Update parameter dropdown when category selection changes
  observeEvent(input$pheno_category, {
    df <- data()
    if (input$pheno_category == "all") {
      params <- sort(unique(df$parameter_combined))
    } else {
      params <- sort(unique(df[df$category == input$pheno_category, ]$parameter_combined))
    }
    updateSelectInput(session, "pheno_parameter", choices = params)
  })
  
  
  # OVERVIEW TAB
  
  output$total_genes <- renderValueBox({
    valueBox(n_distinct(data()$gene_symbol), "Total Genes", icon = icon("dna"), color = "blue")
  })
  
  output$total_parameters <- renderValueBox({
    valueBox(n_distinct(data()$parameter_id), "Total Parameters", icon = icon("list"), color = "green")
  })
  
  output$total_significant <- renderValueBox({
    n_sig <- sum(data()$pvalue <= DEFAULT_THRESHOLD, na.rm = TRUE)
    valueBox(format(n_sig, big.mark = ","), "Significant Associations (p≤0.05)", 
             icon = icon("check"), color = "orange")
  })
  
  output$overview_params <- renderPlotly({
    df <- data() %>%
      group_by(category) %>%
      summarise(n_params = n_distinct(parameter_id), .groups = 'drop') %>%
      arrange(desc(n_params))
    
    plot_ly(df, x = ~reorder(category, n_params), y = ~n_params, type = "bar",
            marker = list(color = ~category, colors = GROUP_COLORS)) %>%
      layout(xaxis = list(title = "", tickangle = -45),
             yaxis = list(title = "Number of Parameters"),
             showlegend = FALSE)
  })
  
  output$overview_sig <- renderPlotly({
    df <- data() %>%
      filter(pvalue <= DEFAULT_THRESHOLD) %>%
      group_by(category) %>%
      summarise(n_sig = n(), .groups = 'drop') %>%
      arrange(desc(n_sig))
    
    plot_ly(df, x = ~reorder(category, n_sig), y = ~n_sig, type = "bar",
            marker = list(color = ~category, colors = GROUP_COLORS)) %>%
      layout(xaxis = list(title = "", tickangle = -45),
             yaxis = list(title = "Significant Associations (p≤0.05)"),
             showlegend = FALSE)
  })
  
  output$overview_table <- renderDT({
    df <- data() %>%
      group_by(category) %>%
      summarise(`Parameters` = n_distinct(parameter_id),
                `Significant (p≤0.05)` = sum(pvalue <= 0.05, na.rm = TRUE),
                `% Significant` = round(sum(pvalue <= 0.05, na.rm = TRUE) / n() * 100, 1),
                .groups = 'drop') %>%
      arrange(desc(`Significant (p≤0.05)`))
    
    datatable(df, options = list(pageLength = 20), rownames = FALSE)
  })
  
  
  # GENE ANALYSIS TAB
  
  # Parse custom gene list from textarea input
  gene_custom_list <- reactive({
    if (is.null(input$gene_custom) || trimws(input$gene_custom) == "") return(NULL)
    genes <- unlist(strsplit(input$gene_custom, "\n"))
    genes <- trimws(genes)
    genes[genes != ""]
  })
  
  # Filter and prepare data for gene-centric analysis
  gene_data <- reactive({
    df <- data()
    
    # Apply gene filter (custom list takes precedence over dropdown)
    custom <- gene_custom_list()
    if (!is.null(custom) && length(custom) > 0) {
      df <- df[df$gene_symbol %in% custom, ]
    } else if (input$gene_select != "all") {
      df <- df[df$gene_combined == input$gene_select, ]
    }
    
    # Apply additional filters
    if (input$gene_category != "all") df <- df[df$category == input$gene_category, ]
    if (input$gene_strain != "all") df <- df[df$mouse_strain == input$gene_strain, ]
    if (input$gene_lifestage != "all") df <- df[df$mouse_life_stage == input$gene_lifestage, ]
    
    if (nrow(df) == 0) return(NULL)
    
    # Combine p-values for same gene-parameter pairs using Fisher's method
    df_combined <- df %>%
      group_by(gene_symbol, parameter_id, parameter_name, category) %>%
      summarise(pvalue = combine_pvalue(pvalue),
                p_adj = combine_pvalue(p_adj),
                procedure_name = first(procedure_name),
                gene_accession_id = first(gene_accession_id),
                .groups = 'drop')
    
    # Add disease information via left join
    disease_summary <- gene_disease_summary()
    df_combined <- df_combined %>%
      left_join(disease_summary, by = "gene_accession_id")
    
    # Select p-value based on FDR correction checkbox
    df_combined$pval_use <- if (input$gene_fdr) df_combined$p_adj else df_combined$pvalue
    df_combined$neg_log10 <- -log10(df_combined$pval_use)
    
    return(df_combined)
  })
  
  output$gene_plot <- renderPlotly({
    df <- gene_data()
    if (is.null(df) || nrow(df) == 0) {
      return(plot_ly() %>% add_annotations(text = "No data. Check filters.", showarrow = FALSE))
    }
    
    df$sig <- ifelse(df$pval_use <= input$gene_pvalue, "Significant", "Not Significant")
    
    plot_ly(df, x = ~seq_along(parameter_name), y = ~neg_log10,
            type = 'scatter', mode = 'markers',
            marker = list(size = 8, color = ~category, colors = GROUP_COLORS,
                          line = list(color = ~ifelse(sig == "Significant", "red", "white"), width = 2)),
            text = ~paste0("<b>Parameter:</b> ", parameter_name,
                           "<br><b>Group:</b> ", category,
                           "<br><b>P-value:</b> ", signif(pval_use, 3),
                           "<br><b>Procedure:</b> ", procedure_name,
                           "<br><b>Diseases:</b> ", ifelse(is.na(n_diseases), 0, n_diseases)),
            hoverinfo = 'text') %>%
      layout(xaxis = list(title = "Parameters (ordered by group)", showticklabels = FALSE),
             yaxis = list(title = paste("-log10(", if (input$gene_fdr) "FDR-adjusted " else "", "p-value)")),
             shapes = list(list(type = "line",
                                y0 = -log10(input$gene_pvalue), y1 = -log10(input$gene_pvalue),
                                x0 = 0, x1 = nrow(df),
                                line = list(color = "red", dash = "dash", width = 2))),
             showlegend = TRUE)
  })
  
  output$gene_table <- renderDT({
    df <- gene_data()
    if (is.null(df)) return(NULL)
    
    df_sig <- df[df$pval_use <= input$gene_pvalue, ]
    if (nrow(df_sig) == 0) {
      return(datatable(data.frame(Message = "No significant phenotypes"),
                       options = list(dom = 't'), rownames = FALSE))
    }
    
    summary <- df_sig %>%
      arrange(pval_use) %>%
      select(`Parameter Name` = parameter_name,
             `Parameter Group` = category,
             `P-value` = pval_use,
             `Disease Link` = n_diseases,
             `Procedure` = procedure_name) %>%
      mutate(`P-value` = format(`P-value`, scientific = TRUE, digits = 3),
             `Disease Link` = ifelse(is.na(`Disease Link`), "No",
                                     paste(`Disease Link`, "diseases")))
    
    datatable(summary, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  
  # PHENOTYPE ANALYSIS TAB
  
  # Display parameter details in sidebar
  output$pheno_details <- renderUI({
    req(input$pheno_parameter)
    df <- data()
    param_info <- df %>%
      filter(parameter_combined == input$pheno_parameter) %>%
      select(parameter_name, parameter_id, category, procedure_name) %>%
      distinct() %>%
      slice(1)
    
    if (nrow(param_info) == 0) return(p("Select a parameter"))
    
    tags$div(
      tags$p(strong("Parameter:"), param_info$parameter_name),
      tags$p(strong("ID:"), param_info$parameter_id),
      tags$p(strong("Group:"), param_info$category),
      tags$p(strong("Procedure:"), param_info$procedure_name)
    )
  })
  
  # Filter and prepare data for phenotype-centric analysis
  pheno_data <- reactive({
    req(input$pheno_parameter)
    df <- data()
    
    df <- df[df$parameter_combined == input$pheno_parameter, ]
    if (nrow(df) == 0) return(NULL)
    
    # Combine p-values per gene using Fisher's method
    df_combined <- df %>%
      group_by(gene_symbol, gene_accession_id) %>%
      summarise(pvalue = combine_pvalue(pvalue),
                p_adj = combine_pvalue(p_adj),
                .groups = 'drop')
    
    # Add disease information
    disease_summary <- gene_disease_summary()
    df_combined <- df_combined %>%
      left_join(disease_summary, by = "gene_accession_id")
    
    df_combined$pval_use <- if (input$pheno_fdr) df_combined$p_adj else df_combined$pvalue
    df_combined$neg_log10 <- -log10(df_combined$pval_use)
    df_combined$has_disease <- !is.na(df_combined$n_diseases) & df_combined$n_diseases > 0
    
    return(df_combined)
  })
  
  output$pheno_plot <- renderPlotly({
    df <- pheno_data()
    if (is.null(df)) {
      return(plot_ly() %>% add_annotations(text = "Select a parameter", showarrow = FALSE))
    }
    
    df$sig <- ifelse(df$pval_use <= input$pheno_pvalue, "Significant", "Not Significant")
    
    # Color points based on disease link or significance
    if (input$pheno_highlight_disease) {
      df$color_group <- ifelse(df$has_disease, "Disease-linked", "No disease link")
      colors <- c("Disease-linked" = "red", "No disease link" = "grey")
    } else {
      df$color_group <- df$sig
      colors <- c("Significant" = "red", "Not Significant" = "grey")
    }
    
    plot_ly(df, x = ~seq_along(gene_symbol), y = ~neg_log10,
            type = 'scatter', mode = 'markers',
            marker = list(size = 8, color = ~color_group, colors = colors,
                          line = list(color = "black", width = 1)),
            text = ~paste0("<b>Gene:</b> ", gene_symbol,
                           "<br><b>MGI ID:</b> ", gene_accession_id,
                           "<br><b>P-value:</b> ", signif(pval_use, 3),
                           "<br><b>Diseases:</b> ", ifelse(is.na(n_diseases), 0, n_diseases),
                           "<br><b>Top disease:</b> ", ifelse(is.na(top_disease), "None", top_disease)),
            hoverinfo = 'text') %>%
      layout(xaxis = list(title = "Genes", showticklabels = FALSE),
             yaxis = list(title = paste("-log10(", if (input$pheno_fdr) "FDR-adjusted " else "", "p-value)")),
             shapes = list(list(type = "line",
                                y0 = -log10(input$pheno_pvalue), y1 = -log10(input$pheno_pvalue),
                                x0 = 0, x1 = nrow(df),
                                line = list(color = "red", dash = "dash", width = 2))),
             showlegend = TRUE)
  })
  
  output$pheno_table <- renderDT({
    df <- pheno_data()
    if (is.null(df)) return(NULL)
    
    df_sig <- df[df$pval_use <= input$pheno_pvalue, ]
    if (nrow(df_sig) == 0) {
      return(datatable(data.frame(Message = "No significant genes"),
                       options = list(dom = 't'), rownames = FALSE))
    }
    
    summary <- df_sig %>%
      arrange(pval_use) %>%
      select(`Gene Symbol` = gene_symbol,
             `MGI ID` = gene_accession_id,
             `P-value` = pval_use,
             `N Diseases` = n_diseases,
             `Top Disease` = top_disease) %>%
      mutate(`P-value` = format(`P-value`, scientific = TRUE, digits = 3),
             `N Diseases` = ifelse(is.na(`N Diseases`), 0, `N Diseases`),
             `Top Disease` = ifelse(is.na(`Top Disease`), "None", `Top Disease`))
    
    datatable(summary, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  
  # CLUSTERING TAB
  
  # Prepare data matrix for UMAP and k-means clustering
  cluster_data <- reactive({
    df <- data()
    
    # Apply filters
    if (input$cluster_category != "all") df <- df[df$category == input$cluster_category, ]
    if (input$cluster_strain != "all") df <- df[df$mouse_strain == input$cluster_strain, ]
    if (input$cluster_lifestage != "all") df <- df[df$mouse_life_stage == input$cluster_lifestage, ]
    
    df$pval_use <- if (input$cluster_fdr) df$p_adj else df$pvalue
    df <- df[df$pval_use <= input$cluster_pvalue, ]
    
    if (nrow(df) == 0) return(NULL)
    
    # Combine p-values and create binary matrix (significant = 1, not = 0)
    df_combined <- df %>%
      group_by(gene_symbol, parameter_id) %>%
      summarise(pvalue = combine_pvalue(pval_use), .groups = 'drop')
    
    mat <- df_combined %>%
      mutate(sig = 1) %>%
      pivot_wider(names_from = parameter_id, values_from = sig, values_fill = 0) %>%
      as.data.frame()
    
    mat <- mat[!duplicated(mat$gene_symbol), ]
    rownames(mat) <- mat$gene_symbol
    mat <- mat[, -1]
    
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    
    # Remove constant columns (no clustering information)
    mat <- mat[, apply(mat, 2, function(x) length(unique(x)) > 1)]
    if (ncol(mat) < 2) return(NULL)
    
    # Configure UMAP with safe n_neighbors (must be < number of samples)
    n_genes <- nrow(mat)
    n_neighbors <- max(2, min(15, n_genes - 1))
    umap_config <- umap::umap.defaults
    umap_config$n_neighbors <- n_neighbors
    umap_config$random_state <- 123
    
    umap_res <- tryCatch(
      umap::umap(mat, config = umap_config),
      error = function(e) NULL
    )
    if (is.null(umap_res)) return(NULL)
    
    # Perform k-means clustering
    k <- min(input$cluster_k, nrow(mat) - 1)
    km_res <- kmeans(mat, centers = k, nstart = 25)
    
    # Combine results into dataframe
    result_df <- data.frame(
      gene = rownames(mat),
      x = umap_res$layout[,1],
      y = umap_res$layout[,2],
      cluster = as.factor(km_res$cluster),
      n_sig = rowSums(mat),
      is_query = rownames(mat) %in% QUERY_GENES
    )
    
    list(df = result_df, matrix = mat)
  })
  
  output$cluster_umap <- renderPlotly({
    obj <- cluster_data()
    if (is.null(obj)) {
      return(
        plot_ly() %>%
          add_annotations(text = "No data. Try relaxing filters.",
                          showarrow = FALSE)
      )
    }
    
    df <- obj$df
    
    # Base scatter plot with genes colored by cluster
    p <- plot_ly(
      data  = df,
      x     = ~x,
      y     = ~y,
      type  = "scatter",
      mode  = "markers",
      color = ~cluster,
      colors = "Set1",
      size  = ~pmin(n_sig, 10),
      sizes = c(4, 14),
      marker = list(opacity = 0.8),
      text = ~paste0(
        "<b>Gene:</b> ", gene,
        "<br><b>Cluster:</b> ", cluster,
        "<br><b>Sig phenotypes:</b> ", n_sig
      ),
      hoverinfo = "text"
    )
    
    # Add highlighting for query genes if checkbox selected
    if (input$cluster_highlight_query) {
      query_df <- df[df$is_query, ]
      if (nrow(query_df) > 0) {
        p <- p %>%
          add_markers(
            data = query_df,
            x    = ~x,
            y    = ~y,
            marker = list(
              size = 10,
              color = "rgba(0,0,0,0)",
              line  = list(color = "red", width = 2)
            ),
            hoverinfo  = "none",
            showlegend = FALSE
          ) %>%
          add_text(
            data = query_df,
            x    = ~x,
            y    = ~y,
            text = ~gene,
            textposition = "top center",
            textfont = list(size = 10, color = "red"),
            showlegend = FALSE
          )
      }
    }
    
    p %>%
      layout(
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2")
      )
  })
  
  output$cluster_heatmap <- renderPlot({
    obj <- cluster_data()
    if (is.null(obj)) {
      plot.new()
      text(0.5, 0.5, "No data", cex = 1.5)
      return()
    }
    
    # Calculate Pearson correlation between genes
    cor_mat <- cor(t(obj$matrix), method = "pearson")
    pheatmap(cor_mat, clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             color = colorRampPalette(c("white", "red"))(100),
             fontsize_row = 8, fontsize_col = 8,
             main = "Gene Correlation Heatmap")
  })
  
  output$cluster_table <- renderDT({
    obj <- cluster_data()
    if (is.null(obj)) return(NULL)
    
    summary <- obj$df %>%
      select(Gene = gene, Cluster = cluster, `Significant Phenotypes` = n_sig,
             `Query Gene` = is_query) %>%
      arrange(Cluster, desc(`Significant Phenotypes`))
    
    datatable(summary, options = list(pageLength = 15), rownames = FALSE) %>%
      formatStyle('Query Gene', backgroundColor = styleEqual(c(TRUE, FALSE), c('yellow', 'white')))
  })
  
  
  # FOUR QUERY GENES TAB
  
  # Summarise significant phenotypes for each query gene
  query_summary <- reactive({
    df <- data()
    df$pval_use <- if (input$query_fdr) df$p_adj else df$pvalue
    
    disease_summary <- gene_disease_summary()
    
    summary <- df %>%
      filter(gene_symbol %in% QUERY_GENES, pval_use <= input$query_pvalue) %>%
      group_by(gene_symbol, gene_accession_id) %>%
      summarise(n_sig = n(),
                top_category = names(sort(table(category), decreasing = TRUE))[1],
                min_pvalue = min(pval_use),
                .groups = 'drop') %>%
      left_join(disease_summary, by = "gene_accession_id")
    
    return(summary)
  })
  
  output$query_smarcd3 <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Smarcd3")
    if (nrow(s) == 0) {
      valueBox(0, "Smarcd3", subtitle = "No significant phenotypes", icon = icon("dna"), color = "red")
    } else {
      valueBox(s$n_sig, "Smarcd3", 
               subtitle = paste("Top:", s$top_category, "| Diseases:", ifelse(is.na(s$n_diseases), 0, s$n_diseases)),
               icon = icon("dna"), color = "red")
    }
  })
  
  output$query_ppp3cc <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Ppp3cc")
    if (nrow(s) == 0) {
      valueBox(0, "Ppp3cc", subtitle = "No significant phenotypes", icon = icon("dna"), color = "blue")
    } else {
      valueBox(s$n_sig, "Ppp3cc",
               subtitle = paste("Top:", s$top_category, "| Diseases:", ifelse(is.na(s$n_diseases), 0, s$n_diseases)),
               icon = icon("dna"), color = "blue")
    }
  })
  
  output$query_rab12 <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Rab12")
    if (nrow(s) == 0) {
      valueBox(0, "Rab12", subtitle = "No significant phenotypes", icon = icon("dna"), color = "green")
    } else {
      valueBox(s$n_sig, "Rab12",
               subtitle = paste("Top:", s$top_category, "| Diseases:", ifelse(is.na(s$n_diseases), 0, s$n_diseases)),
               icon = icon("dna"), color = "green")
    }
  })
  
  output$query_klhl33 <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Klhl33")
    if (nrow(s) == 0) {
      valueBox(0, "Klhl33", subtitle = "No significant phenotypes", icon = icon("dna"), color = "purple")
    } else {
      valueBox(s$n_sig, "Klhl33",
               subtitle = paste("Top:", s$top_category, "| Diseases:", ifelse(is.na(s$n_diseases), 0, s$n_diseases)),
               icon = icon("dna"), color = "purple")
    }
  })
  
  output$query_comparison <- renderPlotly({
    df <- data()
    df$pval_use <- if (input$query_fdr) df$p_adj else df$pvalue
    
    comp <- df %>%
      filter(gene_symbol %in% QUERY_GENES, pval_use <= input$query_pvalue) %>%
      group_by(gene_symbol, category) %>%
      summarise(n = n(), .groups = 'drop')
    
    if (nrow(comp) == 0) {
      return(plot_ly() %>% add_annotations(text = "No significant data", showarrow = FALSE))
    }
    
    plot_ly(comp, x = ~category, y = ~n, color = ~gene_symbol, type = "bar") %>%
      layout(xaxis = list(title = "Parameter Group", tickangle = -45),
             yaxis = list(title = "Number of Significant Phenotypes"),
             barmode = "group")
  })
  
  output$query_table <- renderDT({
    df <- data()
    df$pval_use <- if (input$query_fdr) df$p_adj else df$pvalue
    
    disease_summary <- gene_disease_summary()
    
    result <- df %>%
      filter(gene_symbol %in% QUERY_GENES, pval_use <= input$query_pvalue) %>%
      left_join(disease_summary, by = "gene_accession_id") %>%
      select(Gene = gene_symbol,
             `Parameter Group` = category,
             `Parameter Name` = parameter_name,
             `P-value` = pval_use,
             `Disease Association` = n_diseases,
             `Procedure` = procedure_name) %>%
      mutate(`P-value` = format(`P-value`, scientific = TRUE, digits = 3),
             `Disease Association` = ifelse(is.na(`Disease Association`), "No", 
                                            paste(`Disease Association`, "diseases"))) %>%
      arrange(Gene, `P-value`)
    
    datatable(result, options = list(pageLength = 20, scrollX = TRUE), 
              rownames = FALSE, filter = 'top')
  })
}


# RUN APPLICATION

shinyApp(ui, server)