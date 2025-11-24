library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)
library(umap)
library(metap)

# DATA LOADING
# Make sure the "data_loader_module.R" is in the same directory as the impc_dashboard script
source("data_loader_module.R")

# CONFIGURATION
# Four genes of interest for coursework analysis
QUERY_GENES <- c("Smarcd3", "Ppp3cc", "Rab12", "Klhl33")
DEFAULT_THRESHOLD <- 0.05

# Colours are consistent across all visualisations
GROUP_COLOURS <- c(
  'Housing & Environment' = '#deb887',
  'Structural Phenotype' = '#ff1493',
  'Clinical Chemistry/Blood' = '#ff0000',
  'Embryo & Development' = '#ff69b4',
  'Limb Function & Performance' = '#4169e1',
  'Behavioral & Neurological' = '#00bfff',
  'Cardiovascular & ECG' = '#8b0000',
  'Other' = '#696969'
)

# HELPER FUNCTION - fisher's method for combining p-values
# IMPC tests genes multiple times across batches/conditions, creating multiple p-values for 
# the same gene-parameter pair. Fisher's method properly combines these independent tests into a single summary p-value
combine_pvalue <- function(pvals) {
  pvals[pvals == 0] <- 1e-10
  if (length(pvals) < 2) return(pvals[1])
  else return(sumlog(pvals)$p)
}

# DASHBOARD USER INTERFACE

ui <- dashboardPage(
  dashboardHeader(title = "IMPC Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Analysis", tabName = "gene", icon = icon("dna")),
      menuItem("Phenotype Analysis", tabName = "phenotype", icon = icon("chart-bar")),
      menuItem("Gene Clustering", tabName = "clustering", icon = icon("project-diagram")),
      menuItem("Four Query Genes", tabName = "query_genes", icon = icon("star"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # GENE ANALYSIS TAB
      tabItem(
        tabName = "gene",
        h2("Gene-Centric Analysis"),
        p("Select a gene to visualise all phenotype associations. The plot shows -log10(p-value) for each parameter, 
          with points coloured by parameter group. Significant phenotypes (above the red threshold line) are listed in the table."),
        fluidRow(
          column(3,
                 selectInput("gene_select", "Select Gene:", choices = NULL),
                 sliderInput("gene_pvalue", "P-value Threshold:",
                             min = 0.001, max = 0.1, value = 0.05, step = 0.001)
          ),
          column(9,
                 plotlyOutput("gene_plot", height = "400px"),
                 DTOutput("gene_table"))
        )
      ),
      
      # PHENOTYPE ANALYSIS TAB
      tabItem(
        tabName = "phenotype",
        h2("Phenotype-Centric Analysis"),
        p("Select a phenotype parameter to visualise all gene associations. The plot shows which genes are significantly 
          associated with the selected parameter."),
        fluidRow(
          column(3,
                 selectInput("pheno_category", "Parameter Group:", choices = NULL),
                 selectInput("pheno_parameter", "Parameter:", choices = NULL),
                 hr(),
                 h4("Parameter Details"),
                 uiOutput("pheno_details"),
                 hr(),
                 sliderInput("pheno_pvalue", "P-value Threshold:",
                             min = 0.001, max = 0.1, value = 0.05, step = 0.001)
          ),
          column(9,
                 plotlyOutput("pheno_plot", height = "400px"),
                 DTOutput("pheno_table"))
        )
      ),
      
      # GENE CLUSTERING TAB
      tabItem(
        tabName = "clustering",
        h2("Gene Clustering"),
        p(strong("Purpose:"), " to identify genes with similar phenotype profiles using UMAP dimensionality reduction and k-means clustering. 
          Genes in the same cluster have similar biological effects across tested parameters."),
        fluidRow(
          column(3,
                 sliderInput("cluster_k", "Number of Clusters:",
                             min = 2, max = 10, value = 4, step = 1),
                 sliderInput("cluster_pvalue", "P-value Threshold:",
                             min = 0.001, max = 0.1, value = 0.05, step = 0.001)
          ),
          column(9,
                 plotlyOutput("cluster_umap", height = "400px"),
                 DTOutput("cluster_table"))
        )
      ),
      
      # FOUR QUERY GENES TAB
      tabItem(
        tabName = "query_genes",
        h2("Four Query Genes Summary"),
        p(strong("Purpose:")," to compare significant phenotypes across the four query genes - Smarcd3, 
          Ppp3cc, Rab12, and Klhl33. Shows counts, top biological systems, and parameter group breakdown."),
        fluidRow(
          column(12,
                 sliderInput("query_pvalue", "P-value Threshold:",
                             min = 0.001, max = 0.1, value = 0.05, step = 0.001)
          )
        ),
        fluidRow(
          valueBoxOutput("query_smarcd3"),
          valueBoxOutput("query_ppp3cc"),
          valueBoxOutput("query_rab12"),
          valueBoxOutput("query_klhl33")
        ),
        fluidRow(
          column(12,
                 plotlyOutput("query_comparison", height = "400px"),
                 DTOutput("query_table"))
        )
      )
    )
  )
)

# SERVER
server <- function(input, output, session) {
  
  # Load phenotype data (works with both CSV and MySQL via data_loader_module.R)
  data <- reactive({
    df <- load_data()
    df %>% mutate(is_query = gene_symbol %in% QUERY_GENES)
  })
  
  # Initialise dropdowns
  observe({
    df <- data()
    # Populate gene dropdown
    genes <- c("All" = "all", sort(unique(df$gene_combined)))
    updateSelectInput(session, "gene_select", choices = genes)
    
    # Get actual category names
    actual_categories <- sort(unique(as.character(df$category)))
    cats <- c("All" = "all", setNames(actual_categories, actual_categories))
    updateSelectInput(session, "pheno_category", choices = cats)
  })
  
  # Update parameter dropdown when category changes
  observeEvent(input$pheno_category, {
    df <- data()
    if (is.null(input$pheno_category) || input$pheno_category == "all") {
      params <- sort(unique(df$parameter_name))
    } else {
      params <- sort(unique(df[df$category == input$pheno_category, ]$parameter_name))
    }
    updateSelectInput(session, "pheno_parameter", choices = params)
  }, ignoreNULL = FALSE)
  
  # GENE ANALYSIS
  gene_data <- reactive({
    df <- data()
    if (input$gene_select != "all") {
      df <- df[df$gene_combined == input$gene_select, ]
    }
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      group_by(gene_symbol, parameter_id, parameter_name, category) %>%
      # When the same gene-parameter pair is tested multiple times (different batches, conditions), 
      # combine those p-values using Fisher's method
      summarise(pvalue = combine_pvalue(pvalue),
                procedure_name = first(procedure_name),
                .groups = 'drop') %>%
      mutate(neg_log10 = -log10(pvalue))
  })
  
  output$gene_plot <- renderPlotly({
    df <- gene_data()
    if (is.null(df)) return(plot_ly())
    
    df$sig <- df$pvalue <= input$gene_pvalue
    
    plot_ly(df, x = ~seq_along(parameter_name), y = ~neg_log10,
            type = 'scatter', mode = 'markers',
            color = ~category, colors = GROUP_COLOURS,
            marker = list(size = 6),
            text = ~paste0(parameter_name, "\n", category, "\nP: ", signif(pvalue, 3)),
            hoverinfo = 'text') %>%
      layout(xaxis = list(title = "Parameters"),
             yaxis = list(title = "-log10(p-value)"),
             shapes = list(list(type = "line",
                                y0 = -log10(input$gene_pvalue), 
                                y1 = -log10(input$gene_pvalue),
                                x0 = 0, x1 = nrow(df),
                                line = list(color = "red", dash = "dash"))))
  })
  
  output$gene_table <- renderDT({
    df <- gene_data()
    if (is.null(df)) return(NULL)
    
    df %>%
      filter(pvalue <= input$gene_pvalue) %>%
      arrange(pvalue) %>%
      select(Parameter = parameter_name, Group = category, 
             `P-value` = pvalue, Procedure = procedure_name) %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })
  
  # PHENOTYPE ANALYSIS
  
  # Display parameter details in sidebar
  output$pheno_details <- renderUI({
    req(input$pheno_parameter)
    df <- data()
    param_info <- df %>%
      filter(parameter_name == input$pheno_parameter) %>%
      select(parameter_name, parameter_id, category) %>%
      distinct() %>%
      slice(1)
    
    if (nrow(param_info) == 0) return(p("Select a parameter"))
    
    tags$div(
      tags$p(strong("Name:"), param_info$parameter_name),
      tags$p(strong("ID:"), param_info$parameter_id),
      tags$p(strong("Group:"), param_info$category)
    )
  })
  
  pheno_data <- reactive({
    req(input$pheno_parameter)
    df <- data()
    df <- df[df$parameter_name == input$pheno_parameter, ]
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      group_by(gene_symbol) %>%
      # combine multiple p-values for the same gene when analysing a specific parameter
      summarise(pvalue = combine_pvalue(pvalue), .groups = 'drop') %>%
      mutate(neg_log10 = -log10(pvalue))
  })
  
  output$pheno_plot <- renderPlotly({
    df <- pheno_data()
    if (is.null(df)) return(plot_ly())
    
    df$sig <- df$pvalue <= input$pheno_pvalue
    
    plot_ly(df, x = ~seq_along(gene_symbol), y = ~neg_log10,
            type = 'scatter', mode = 'markers',
            marker = list(size = 6, color = ~ifelse(sig, "red", "grey")),
            text = ~paste0(gene_symbol, "\nP: ", signif(pvalue, 3)),
            hoverinfo = 'text') %>%
      layout(xaxis = list(title = "Genes"),
             yaxis = list(title = "-log10(p-value)"),
             shapes = list(list(type = "line",
                                y0 = -log10(input$pheno_pvalue),
                                y1 = -log10(input$pheno_pvalue),
                                x0 = 0, x1 = nrow(df),
                                line = list(color = "red", dash = "dash"))))
  })
  
  output$pheno_table <- renderDT({
    df <- pheno_data()
    if (is.null(df)) return(NULL)
    
    df %>%
      filter(pvalue <= input$pheno_pvalue) %>%
      arrange(pvalue) %>%
      select(Gene = gene_symbol, `P-value` = pvalue) %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })
  
  # CLUSTERING
  cluster_data <- reactive({
    df <- data()
    df <- df[df$pvalue <= input$cluster_pvalue, ]
    if (nrow(df) == 0) return(NULL)
    
    df_combined <- df %>%
      group_by(gene_symbol, parameter_id) %>%
      # Before creating the clustering matrix, combine p-values for any gene-parameter pairs tested multiple times
      summarise(pvalue = combine_pvalue(pvalue), .groups = 'drop')
    
    # Create binary matrix (1 = significant, 0 = not significant)
    mat <- df_combined %>%
      mutate(sig = 1) %>%
      pivot_wider(names_from = parameter_id, values_from = sig, values_fill = 0) %>%
      as.data.frame()
    
    # Prepare matrix for clustering
    mat <- mat[!duplicated(mat$gene_symbol), ]
    rownames(mat) <- mat$gene_symbol
    mat <- mat[, -1]
    
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    # Remove constant columns (no clustering information)
    mat <- mat[, apply(mat, 2, function(x) length(unique(x)) > 1)]
    if (ncol(mat) < 2) return(NULL)
    
    # Configure UMAP (adaptive n_neighbors for different sample sizes)
    n_neighbors <- max(2, min(15, nrow(mat) - 1))
    umap_config <- umap::umap.defaults
    umap_config$n_neighbors <- n_neighbors
    umap_config$random_state <- 123
    
    # Run UMAP dimensionality reduction
    umap_res <- tryCatch(umap::umap(mat, config = umap_config), error = function(e) NULL)
    if (is.null(umap_res)) return(NULL)
    
    # K-means clustering on original binary matrix
    k <- min(input$cluster_k, nrow(mat) - 1)
    km_res <- kmeans(mat, centers = k, nstart = 25)
    
    # Combine results into data frame for visualisation
    data.frame(
      gene = rownames(mat),
      x = umap_res$layout[,1],
      y = umap_res$layout[,2],
      cluster = as.factor(km_res$cluster),
      n_sig = as.integer(rowSums(mat)),
      is_query = rownames(mat) %in% QUERY_GENES
    )
  })
  
  output$cluster_umap <- renderPlotly({
    df <- cluster_data()
    if (is.null(df)) return(plot_ly())
    
    p <- plot_ly(df, x = ~x, y = ~y, color = ~cluster, colors = "Set1",
                 type = "scatter", mode = "markers",
                 size = ~n_sig, sizes = c(10, 25),
                 text = ~paste0(gene, "\nCluster: ", cluster, "\nSig Phenotypes: ", n_sig),
                 hoverinfo = 'text') %>%
      layout(xaxis = list(title = "UMAP 1"),
             yaxis = list(title = "UMAP 2"))
    
    query_df <- df[df$is_query, ]
    if (nrow(query_df) > 0) {
      p <- p %>%
        add_text(data = query_df, x = ~x, y = ~y, text = ~gene,
                 textfont = list(color = "red"), showlegend = FALSE)
    }
    p
  })
  
  output$cluster_table <- renderDT({
    df <- cluster_data()
    if (is.null(df)) return(NULL)
    
    df %>%
      select(Gene = gene, Cluster = cluster, `Sig Phenotypes` = n_sig) %>%
      arrange(Cluster) %>%
      datatable(options = list(pageLength = 15), rownames = FALSE)
  })
  
  # QUERY GENES
  query_summary <- reactive({
    df <- data()
    df %>%
      filter(gene_symbol %in% QUERY_GENES, pvalue <= input$query_pvalue) %>%
      group_by(gene_symbol) %>%
      summarise(
        n_sig = n(), # Count of significant phenotypes
        # Find most common parameter group for this gene
                top_category = names(sort(table(category), decreasing = TRUE))[1],
                .groups = 'drop')
  })
  
  output$query_smarcd3 <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Smarcd3")
    if (nrow(s) == 0) {
      valueBox(0, "Smarcd3", subtitle = "No significant phenotypes", 
               icon = icon("dna"), color = "red")
    } else {
      valueBox(s$n_sig, "Smarcd3", 
               subtitle = paste("Top:", s$top_category),
               icon = icon("dna"), color = "red")
    }
  })
  
  output$query_ppp3cc <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Ppp3cc")
    if (nrow(s) == 0) {
      valueBox(0, "Ppp3cc", subtitle = "No significant phenotypes", 
               icon = icon("dna"), color = "blue")
    } else {
      valueBox(s$n_sig, "Ppp3cc",
               subtitle = paste("Top:", s$top_category),
               icon = icon("dna"), color = "blue")
    }
  })
  
  output$query_rab12 <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Rab12")
    if (nrow(s) == 0) {
      valueBox(0, "Rab12", subtitle = "No significant phenotypes", 
               icon = icon("dna"), color = "green")
    } else {
      valueBox(s$n_sig, "Rab12",
               subtitle = paste("Top:", s$top_category),
               icon = icon("dna"), color = "green")
    }
  })
  
  output$query_klhl33 <- renderValueBox({
    s <- query_summary() %>% filter(gene_symbol == "Klhl33")
    if (nrow(s) == 0) {
      valueBox(0, "Klhl33", subtitle = "No significant phenotypes", 
               icon = icon("dna"), color = "purple")
    } else {
      valueBox(s$n_sig, "Klhl33",
               subtitle = paste("Top:", s$top_category),
               icon = icon("dna"), color = "purple")
    }
  })
  
  output$query_comparison <- renderPlotly({
    df <- data()
    comp <- df %>%
      filter(gene_symbol %in% QUERY_GENES, pvalue <= input$query_pvalue) %>%
      group_by(gene_symbol, category) %>%
      summarise(n = n(), .groups = 'drop')
    
    if (nrow(comp) == 0) return(plot_ly())
    
    plot_ly(comp, x = ~category, y = ~n, color = ~gene_symbol, type = "bar") %>%
      layout(xaxis = list(title = "Parameter Group", tickangle = -45),
             yaxis = list(title = "Significant Phenotypes"),
             barmode = "group")
  })
  
  output$query_table <- renderDT({
    df <- data()
    df %>%
      filter(gene_symbol %in% QUERY_GENES, pvalue <= input$query_pvalue) %>%
      select(Gene = gene_symbol, Group = category, 
             Parameter = parameter_name, `P-value` = pvalue) %>%
      arrange(Gene, `P-value`) %>%
      datatable(options = list(pageLength = 15), rownames = FALSE)
  })
}

shinyApp(ui, server)