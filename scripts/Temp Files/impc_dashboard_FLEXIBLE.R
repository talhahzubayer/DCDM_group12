# ============================================================================
# IMPC Phenotype Analysis Dashboard - FLEXIBLE VERSION
# Works with BOTH CSV files and MySQL database
# Just change DATA_SOURCE in data_loader_module.R!
# ============================================================================

library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)
library(umap)
library(metap)
library(pheatmap)
library(viridis)
library(RColorBrewer)

# ============================================================================
# LOAD DATA MODULE
# ============================================================================

# Source the flexible data loader
# Make sure data_loader_module.R is in the same directory!
source("data_loader_module.R")

# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_PVALUE <- 0.05

# Color palette for parameter groups
GROUP_COLORS <- c(
  'Weight' = '#ff1493', 'Images' = '#dda0dd', 'Brain' = '#00bfff',
  'Blood' = '#ff0000', 'Vision/Eye' = '#ffff00', 'Cardiovascular' = '#8b0000',
  'Metabolic' = '#ff4500', 'Respiratory' = '#00fa9a', 'Muscular' = '#4169e1',
  'Reproductive' = '#ff00ff', 'Coat/Skin' = '#ffa500', 'Biochemical' = '#00ffff',
  'Equipment' = '#f0e68c', 'Housing' = '#deb887', 'Conditions' = '#e9967a',
  'Other' = '#696969', 'Non-significant' = 'gray86'
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Fisher's method for p-value combination
combine_pvalue <- function(pvals) {
  pvals[pvals == 0] <- 1e-10
  if (length(pvals) < 2) return(pvals[1])
  else return(sumlog(pvals)$p)
}

# ============================================================================
# USER INTERFACE
# ============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "IMPC Dashboard - Group 12"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Analysis", tabName = "gene_view", icon = icon("dna")),
      menuItem("Phenotype Analysis", tabName = "phenotype_view", icon = icon("chart-bar")),
      menuItem("Gene Clustering", tabName = "clustering_view", icon = icon("project-diagram")),
      br(),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    # Data source indicator
    tags$div(
      style = "position: fixed; bottom: 10px; left: 10px; background: #222; padding: 10px; border-radius: 5px; color: white; font-size: 10px;",
      tags$strong("Data Source:"),
      tags$br(),
      textOutput("data_source_indicator", inline = TRUE)
    )
  ),
  
  dashboardBody(
    tabItems(
      # Gene Analysis Tab
      tabItem(
        tabName = "gene_view",
        h2("Visualization 1: Gene-Centric Analysis"),
        p("Select genes to view all phenotypes tested and their statistical significance."),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Gene Selection"),
              selectInput("gene_select", "Select Gene:", choices = NULL),
              textAreaInput("custom_genes",
                           "Or paste custom gene list (one per line):",
                           rows = 4,
                           placeholder = "Smarcd3\nPpp3cc\nRab12\nKlhl33"),
              hr(),
              h4("Filters"),
              selectInput("gene_category_filter", "Parameter Category:", choices = NULL),
              selectInput("gene_strain_filter", "Mouse Strain:", choices = NULL),
              selectInput("gene_lifestage_filter", "Life Stage:", choices = NULL),
              numericInput("gene_pvalue_threshold", "P-value Threshold:",
                          value = DEFAULT_PVALUE, min = 0, max = 1, step = 0.01)
            )
          ),
          column(9,
            box(width = 12, title = "Phenotype Significance Plot", 
                status = "primary", solidHeader = TRUE,
                plotlyOutput("gene_plot", height = "500px")),
            box(width = 12, title = "Significant Phenotypes Table",
                status = "info", solidHeader = TRUE,
                DTOutput("gene_table"))
          )
        )
      ),
      
      # Phenotype Analysis Tab
      tabItem(
        tabName = "phenotype_view",
        h2("Visualization 2: Phenotype-Centric Analysis"),
        p("Select a phenotype to view all genes tested and their statistical significance."),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Phenotype Selection"),
              selectInput("phenotype_category_filter", "Parameter Category:", choices = NULL),
              selectInput("parameter_select", "Select Parameter:", choices = NULL),
              hr(),
              h4("Filters"),
              numericInput("phenotype_pvalue_threshold", "P-value Threshold:",
                          value = DEFAULT_PVALUE, min = 0, max = 1, step = 0.01)
            )
          ),
          column(9,
            box(width = 12, title = "Gene Significance Plot",
                status = "primary", solidHeader = TRUE,
                plotlyOutput("phenotype_plot", height = "500px")),
            box(width = 12, title = "Significant Genes Table",
                status = "info", solidHeader = TRUE,
                DTOutput("phenotype_table"))
          )
        )
      ),
      
      # Clustering Tab
      tabItem(
        tabName = "clustering_view",
        h2("Visualization 3: Gene Clustering Analysis"),
        p("Identify genes with similar phenotype profiles using UMAP dimensionality reduction."),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Clustering Options"),
              selectInput("cluster_category_filter", "Parameter Category:", choices = NULL),
              selectInput("cluster_strain_filter", "Mouse Strain:", choices = NULL),
              selectInput("cluster_lifestage_filter", "Life Stage:", choices = NULL),
              numericInput("cluster_pvalue_threshold", "P-value Threshold:",
                          value = DEFAULT_PVALUE, min = 0, max = 1, step = 0.01),
              sliderInput("k_clusters", "Number of K-means Clusters:",
                         min = 2, max = 10, value = 5, step = 1),
              hr(),
              actionButton("reset_cluster", "Reset Filters",
                          icon = icon("refresh"), class = "btn-primary"),
              br(), br(),
              downloadButton("download_cluster", "Download Results", class = "btn-success")
            )
          ),
          column(9,
            box(width = 12, title = "UMAP Clustering Plot",
                status = "primary", solidHeader = TRUE,
                plotlyOutput("cluster_plot", height = "500px")),
            box(width = 12, title = "Gene Correlation Heatmap",
                status = "primary", solidHeader = TRUE,
                plotOutput("heatmap_plot", height = "600px"))
          )
        )
      ),
      
      # About Tab
      tabItem(
        tabName = "about",
        h2("About This Dashboard"),
        box(width = 12,
            h3("IMPC Phenotype Analysis Dashboard"),
            p("Interactive exploration of IMPC mouse phenotype data."),
            
            h4("Flexible Data Source:"),
            p("This dashboard can load data from:"),
            tags$ul(
              tags$li(strong("CSV Files"), " - For development and testing"),
              tags$li(strong("MySQL Database"), " - For production deployment")
            ),
            p("Change data source by updating DATA_SOURCE in data_loader_module.R"),
            
            h4("Your Query Genes:"),
            tags$ul(
              tags$li(strong("Smarcd3"), " (MGI:1914243) - SWI/SNF related matrix associated actin dependent regulator of chromatin"),
              tags$li(strong("Ppp3cc"), " (MGI:107162) - Protein phosphatase 3 catalytic subunit gamma"),
              tags$li(strong("Rab12"), " (MGI:894284) - RAB12 member RAS oncogene family"),
              tags$li(strong("Klhl33"), " (MGI:3644593) - Kelch like family member 33")
            ),
            
            h4("Features:"),
            tags$ul(
              tags$li("Gene-centric visualization"),
              tags$li("Phenotype-centric visualization"),
              tags$li("UMAP clustering analysis"),
              tags$li("Fisher's method for p-value combination"),
              tags$li("Automated parameter categorization")
            ),
            
            h4("Data Statistics:"),
            uiOutput("data_stats")
        )
      )
    )
  )
)

# ============================================================================
# SERVER LOGIC
# ============================================================================

server <- function(input, output, session) {
  
  # Load data on startup using flexible loader
  data <- reactive({
    load_data()  # Uses DATA_SOURCE from data_loader_module.R
  })
  
  # Display data source in UI
  output$data_source_indicator <- renderText({
    toupper(DATA_SOURCE)
  })
  
  # Display data statistics in About tab
  output$data_stats <- renderUI({
    df <- data()
    tags$div(
      tags$p(strong("Total Records: "), format(nrow(df), big.mark = ",")),
      tags$p(strong("Unique Genes: "), n_distinct(df$gene_symbol)),
      tags$p(strong("Unique Parameters: "), n_distinct(df$parameter_id)),
      tags$p(strong("Mouse Strains: "), paste(levels(df$mouse_strain), collapse = ", ")),
      tags$p(strong("Life Stages: "), paste(levels(df$mouse_life_stage), collapse = ", "))
    )
  })
  
  # Initialize dropdown choices
  observe({
    df <- data()
    
    gene_choices <- c("All" = "all", sort(unique(df$gene_combined)))
    updateSelectInput(session, "gene_select", choices = gene_choices)
    
    category_choices <- c("All" = "all", sort(unique(as.character(df$category))))
    updateSelectInput(session, "gene_category_filter", choices = category_choices)
    updateSelectInput(session, "phenotype_category_filter", choices = category_choices)
    updateSelectInput(session, "cluster_category_filter", choices = category_choices)
    
    strain_choices <- c("All" = "all", sort(unique(as.character(df$mouse_strain))))
    updateSelectInput(session, "gene_strain_filter", choices = strain_choices)
    updateSelectInput(session, "cluster_strain_filter", choices = strain_choices)
    
    lifestage_choices <- c("All" = "all", sort(unique(as.character(df$mouse_life_stage))))
    updateSelectInput(session, "gene_lifestage_filter", choices = lifestage_choices)
    updateSelectInput(session, "cluster_lifestage_filter", choices = lifestage_choices)
    
    param_choices <- c("All" = "all", sort(unique(df$parameter_combined)))
    updateSelectInput(session, "parameter_select", choices = param_choices)
  })
  
  # Update parameter choices when category filter changes
  observeEvent(input$phenotype_category_filter, {
    df <- data()
    if (input$phenotype_category_filter == "all") {
      param_choices <- sort(unique(df$parameter_combined))
    } else {
      param_choices <- sort(unique(df[df$category == input$phenotype_category_filter, ]$parameter_combined))
    }
    updateSelectInput(session, "parameter_select", choices = c("All" = "all", param_choices))
  })
  
  # Process custom gene list
  custom_genes <- reactive({
    req(input$custom_genes)
    if (input$custom_genes == "") return(NULL)
    genes <- unlist(strsplit(input$custom_genes, "\n"))
    genes <- trimws(genes)
    genes[genes != ""]
  })
  
  # Filter data for gene view
  gene_data <- reactive({
    df <- data()
    
    custom <- custom_genes()
    if (!is.null(custom) && length(custom) > 0) {
      df <- df[df$gene_symbol %in% custom, ]
    } else if (input$gene_select != "all") {
      df <- df[df$gene_combined == input$gene_select, ]
    }
    
    if (input$gene_category_filter != "all") {
      df <- df[df$category == input$gene_category_filter, ]
    }
    if (input$gene_strain_filter != "all") {
      df <- df[df$mouse_strain == input$gene_strain_filter, ]
    }
    if (input$gene_lifestage_filter != "all") {
      df <- df[df$mouse_life_stage == input$gene_lifestage_filter, ]
    }
    
    df %>%
      group_by(gene_symbol, parameter_name, category, parameter_id) %>%
      summarise(
        pvalue = combine_pvalue(pvalue),
        neg_log10_pvalue = -log10(pvalue),
        .groups = 'drop'
      )
  })
  
  # Gene view plot
  output$gene_plot <- renderPlotly({
    df <- gene_data()
    if (nrow(df) == 0) {
      return(plot_ly() %>%
        add_annotations(text = "No data available for selected criteria", 
                       showarrow = FALSE, font = list(size = 16)))
    }
    
    df$significant <- ifelse(df$pvalue <= input$gene_pvalue_threshold, "Significant", "Not Significant")
    
    plot_ly(df, x = ~seq_along(parameter_name), y = ~neg_log10_pvalue,
            type = 'scatter', mode = 'markers',
            marker = list(size = 8, color = ~pvalue, colorscale = 'Viridis',
                         colorbar = list(title = "P-value"),
                         line = list(color = ~ifelse(significant == "Significant", "red", "white"), width = 2)),
            text = ~paste("Parameter:", parameter_name, "<br>Category:", category,
                         "<br>P-value:", format(pvalue, scientific = TRUE, digits = 3)),
            hoverinfo = 'text') %>%
      layout(title = "Phenotype Significance for Selected Gene(s)",
             xaxis = list(title = "Parameters", showticklabels = FALSE),
             yaxis = list(title = "-log10(p-value)"),
             shapes = list(list(type = "line",
                               y0 = -log10(input$gene_pvalue_threshold),
                               y1 = -log10(input$gene_pvalue_threshold),
                               x0 = 0, x1 = nrow(df),
                               line = list(color = "red", dash = "dash", width = 2))))
  })
  
  # Gene view table
  output$gene_table <- renderDT({
    df <- gene_data()
    df_sig <- df[df$pvalue <= input$gene_pvalue_threshold, ]
    
    if (nrow(df_sig) == 0) {
      return(datatable(data.frame(Message = "No significant phenotypes found. Try adjusting filters."),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    summary_df <- df_sig %>%
      arrange(pvalue) %>%
      select(Parameter = parameter_name, Category = category, `P-value` = pvalue) %>%
      mutate(`P-value` = format(`P-value`, scientific = TRUE, digits = 3))
    
    datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE,
                                         order = list(list(2, 'asc'))), rownames = FALSE)
  })
  
  # Phenotype view data
  phenotype_data <- reactive({
    df <- data()
    if (input$parameter_select != "all") {
      df <- df[df$parameter_combined == input$parameter_select, ]
    } else if (input$phenotype_category_filter != "all") {
      df <- df[df$category == input$phenotype_category_filter, ]
    }
    
    df %>%
      group_by(gene_symbol, gene_accession_id, parameter_name, mouse_strain) %>%
      summarise(pvalue = combine_pvalue(pvalue), neg_log10_pvalue = -log10(pvalue), .groups = 'drop')
  })
  
  # Phenotype view plot
  output$phenotype_plot <- renderPlotly({
    df <- phenotype_data()
    if (nrow(df) == 0) {
      return(plot_ly() %>%
        add_annotations(text = "No data available", showarrow = FALSE, font = list(size = 16)))
    }
    
    df$significant <- ifelse(df$pvalue <= input$phenotype_pvalue_threshold, "Significant", "Not Significant")
    
    plot_ly(df, x = ~seq_along(gene_symbol), y = ~neg_log10_pvalue,
            type = 'scatter', mode = 'markers',
            marker = list(size = 8, color = ~pvalue, colorscale = 'Plasma',
                         colorbar = list(title = "P-value"),
                         line = list(color = ~ifelse(significant == "Significant", "red", "white"), width = 2)),
            text = ~paste("Gene:", gene_symbol, "<br>MGI ID:", gene_accession_id,
                         "<br>P-value:", format(pvalue, scientific = TRUE, digits = 3)),
            hoverinfo = 'text') %>%
      layout(title = "Gene Significance for Selected Phenotype(s)",
             xaxis = list(title = "Genes", showticklabels = FALSE),
             yaxis = list(title = "-log10(p-value)"),
             shapes = list(list(type = "line",
                               y0 = -log10(input$phenotype_pvalue_threshold),
                               y1 = -log10(input$phenotype_pvalue_threshold),
                               x0 = 0, x1 = nrow(df),
                               line = list(color = "red", dash = "dash", width = 2))))
  })
  
  # Phenotype view table
  output$phenotype_table <- renderDT({
    df <- phenotype_data()
    df_sig <- df[df$pvalue <= input$phenotype_pvalue_threshold, ]
    
    if (nrow(df_sig) == 0) {
      return(datatable(data.frame(Message = "No significant genes found"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    summary_df <- df_sig %>%
      arrange(pvalue) %>%
      select(Gene = gene_symbol, `MGI ID` = gene_accession_id,
             Strain = mouse_strain, `P-value` = pvalue) %>%
      mutate(`P-value` = format(`P-value`, scientific = TRUE, digits = 3))
    
    datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE,
                                         order = list(list(3, 'asc'))), rownames = FALSE)
  })
  
  # Clustering data
  cluster_data <- reactive({
    df <- data()
    
    if (input$cluster_category_filter != "all") {
      df <- df[df$category == input$cluster_category_filter, ]
    }
    if (input$cluster_strain_filter != "all") {
      df <- df[df$mouse_strain == input$cluster_strain_filter, ]
    }
    if (input$cluster_lifestage_filter != "all") {
      df <- df[df$mouse_life_stage == input$cluster_lifestage_filter, ]
    }
    
    df <- df[df$pvalue <= input$cluster_pvalue_threshold, ]
    
    if (nrow(df) == 0) return(NULL)
    
    df_combined <- df %>%
      group_by(gene_symbol, parameter_id) %>%
      summarise(pvalue = combine_pvalue(pvalue), .groups = 'drop')
    
    mat <- df_combined %>%
      mutate(significant = 1) %>%
      pivot_wider(names_from = parameter_id, values_from = significant, values_fill = 0) %>%
      as.data.frame()
    
    rownames(mat) <- mat$gene_symbol
    mat <- mat[, -1]
    
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    
    list(matrix = mat, genes = rownames(mat), raw_data = df_combined)
  })
  
  # Clustering plot
  output$cluster_plot <- renderPlotly({
    clust <- cluster_data()
    if (is.null(clust)) {
      return(plot_ly() %>%
        add_annotations(text = "No significant phenotypes found.\nTry adjusting filters.", 
                       showarrow = FALSE, font = list(size = 16)))
    }
    
    umap_result <- umap(clust$matrix, random_state = 123)
    k <- min(input$k_clusters, nrow(clust$matrix) - 1)
    kmeans_result <- kmeans(clust$matrix, centers = k, nstart = 25)
    sig_counts <- rowSums(clust$matrix)
    
    plot_df <- data.frame(x = umap_result$layout[, 1], y = umap_result$layout[, 2],
                         gene = clust$genes, cluster = as.factor(kmeans_result$cluster),
                         sig_count = sig_counts)
    
    plot_ly(plot_df, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
            marker = list(size = ~sig_count * 2, color = ~as.numeric(cluster),
                         colorscale = 'Viridis', colorbar = list(title = "Cluster"),
                         line = list(color = "white", width = 1)),
            text = ~paste("Gene:", gene, "<br>Cluster:", cluster, "<br>Significant Phenotypes:", sig_count),
            hoverinfo = 'text') %>%
      layout(title = "UMAP Clustering of Genes", xaxis = list(title = "UMAP 1"),
             yaxis = list(title = "UMAP 2"))
  })
  
  # Heatmap
  output$heatmap_plot <- renderPlot({
    clust <- cluster_data()
    if (is.null(clust)) {
      plot.new()
      text(0.5, 0.5, "No significant phenotypes found.\nTry adjusting filters.", cex = 1.2)
    } else {
      cor_mat <- cor(t(clust$matrix), method = "pearson")
      pheatmap(cor_mat, clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean", clustering_method = "complete",
               color = colorRampPalette(c("white", "red"))(100), show_rownames = TRUE,
               show_colnames = TRUE, fontsize_row = 8, fontsize_col = 8,
               main = "Gene Correlation Heatmap")
    }
  })
  
  # Reset filters
  observeEvent(input$reset_cluster, {
    updateSelectInput(session, "cluster_category_filter", selected = "all")
    updateSelectInput(session, "cluster_strain_filter", selected = "all")
    updateSelectInput(session, "cluster_lifestage_filter", selected = "all")
    updateNumericInput(session, "cluster_pvalue_threshold", value = DEFAULT_PVALUE)
    updateSliderInput(session, "k_clusters", value = 5)
  })
  
  # Download results
  output$download_cluster <- downloadHandler(
    filename = function() paste0("clustering_results_", Sys.Date(), ".csv"),
    content = function(file) {
      clust <- cluster_data()
      if (!is.null(clust)) {
        results <- data.frame(Gene = clust$genes, Significant_Phenotypes = rowSums(clust$matrix))
        write.csv(results, file, row.names = FALSE)
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
