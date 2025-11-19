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


# LOAD DATA MODULE

# Source the flexible data loader
source("data_loader_module.R")


# CONFIGURATION

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


# HELPER FUNCTIONS

# Fisher's method for p-value combination
combine_pvalue <- function(pvals) {
  pvals[pvals == 0] <- 1e-10
  if (length(pvals) < 2) return(pvals[1])
  else return(sumlog(pvals)$p)
}


# USER INTERFACE

ui <- dashboardPage(
  dashboardHeader(title = "IMPC Dashboard - Group 12"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Analysis", tabName = "gene_view", icon = icon("dna")),
      menuItem("Bubble Chart", tabName = "bubble_view", icon = icon("circle")),
      menuItem("Phenotype Analysis", tabName = "phenotype_view", icon = icon("chart-bar")),
      menuItem("Gene Clustering", tabName = "clustering_view", icon = icon("project-diagram")),
      menuItem("Procedure Info", tabName = "procedure_view", icon = icon("table")),
      menuItem("Disease Associations", tabName = "disease_view", icon = icon("notes-medical")),
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
        h2("Visualisation 1: Gene-Centric Analysis"),
        p("Select genes to view all phenotypes tested and their statistical significance."),
        p(tags$em("Hover over points to see procedure information!")),
        
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
              checkboxInput("gene_only_mandatory", "Only Mandatory Procedures", value = FALSE),
              hr(),
              h4("P-value Options"),
              sliderInput("gene_pvalue_threshold", "P-value Threshold:",
                         min = 0.001, max = 0.1, value = 0.05, step = 0.001),
              checkboxInput("gene_use_fdr", "Use FDR-adjusted p-values", value = FALSE)
            )
          ),
          column(9,
            box(width = 12, title = "Gene Phenotype Plot",
                status = "primary", solidHeader = TRUE,
                plotlyOutput("gene_plot", height = "500px")),
            box(width = 12, title = "Significant Phenotypes",
                status = "info", solidHeader = TRUE,
                DTOutput("gene_table"))
          )
        )
      ),
      
      # Bubble Chart Tab
      tabItem(
        tabName = "bubble_view",
        h2("Bubble Chart: Gene Phenotype Significance"),
        p("Visualize the significance of multiple parameters for a selected gene."),
        p(tags$em("Bubble size is proportional to significance. Larger bubbles = more significant.")),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Gene Selection"),
              selectInput("bubble_gene_select", "Select Gene:", choices = NULL),
              hr(),
              h4("Display Options"),
              sliderInput("bubble_pvalue_threshold", "P-value Threshold:",
                         min = 0.001, max = 0.1, value = 0.05, step = 0.001),
              checkboxInput("bubble_use_fdr", "Use FDR-adjusted p-values", value = FALSE),
              selectInput("bubble_category_filter", "Parameter Category:", choices = NULL)
            )
          ),
          column(9,
            box(width = 12, title = "Bubble Chart",
                status = "primary", solidHeader = TRUE,
                plotlyOutput("bubble_chart", height = "600px")),
            box(width = 12, title = "Chart Explanation",
                status = "info", solidHeader = TRUE,
                p("This bubble chart shows:"),
                tags$ul(
                  tags$li(strong("Y-axis:"), " Transformed p-value (higher = more significant)"),
                  tags$li(strong("Bubble size:"), " Proportional to significance"),
                  tags$li(strong("Color:"), " Blue = significant, Grey = non-significant"),
                  tags$li(strong("Transformation:"), " (1 - p-value)^8 for better Visualisation")
                ))
          )
        )
      ),
      
      # Phenotype Analysis Tab
      tabItem(
        tabName = "phenotype_view",
        h2("Visualisation 2: Phenotype-Centric Analysis"),
        p("Select phenotypes to see which genes are significantly associated."),
        p(tags$em("Hover over points to see procedure information!")),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Phenotype Selection"),
              selectInput("phenotype_category_filter", "Parameter Category:", choices = NULL),
              selectInput("parameter_select", "Select Parameter:", choices = NULL),
              hr(),
              h4("P-value Options"),
              sliderInput("phenotype_pvalue_threshold", "P-value Threshold:",
                         min = 0.001, max = 0.1, value = 0.05, step = 0.001),
              checkboxInput("phenotype_use_fdr", "Use FDR-adjusted p-values", value = FALSE)
            )
          ),
          column(9,
            box(width = 12, title = "Gene Significance Plot",
                status = "primary", solidHeader = TRUE,
                plotlyOutput("phenotype_plot", height = "500px")),
            box(width = 12, title = "Significant Genes",
                status = "info", solidHeader = TRUE,
                DTOutput("phenotype_table"))
          )
        )
      ),
      
      # Gene Clustering Tab
      tabItem(
        tabName = "clustering_view",
        h2("Visualisation 3: Gene Clustering Analysis"),
        p("UMAP clustering of genes based on phenotype similarity profiles."),
        p(tags$em("Enhanced with most significant parameter information!")),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Clustering Parameters"),
              sliderInput("k_clusters", "Number of Clusters:",
                         min = 2, max = 10, value = 3, step = 1),
              hr(),
              h4("Filters"),
              selectInput("cluster_category_filter", "Parameter Category:", choices = NULL),
              selectInput("cluster_strain_filter", "Mouse Strain:", choices = NULL),
              selectInput("cluster_lifestage_filter", "Life Stage:", choices = NULL),
              checkboxInput("cluster_only_mandatory", "Only Mandatory Procedures", value = FALSE),
              hr(),
              h4("P-value Options"),
              sliderInput("cluster_pvalue_threshold", "P-value Threshold:",
                         min = 0.001, max = 0.1, value = 0.05, step = 0.001),
              checkboxInput("cluster_use_fdr", "Use FDR-adjusted p-values", value = FALSE)
            )
          ),
          column(9,
            box(width = 12, title = "UMAP Clustering",
                status = "primary", solidHeader = TRUE,
                plotlyOutput("cluster_plot", height = "500px")),
            box(width = 12, title = "Correlation Heatmap",
                status = "warning", solidHeader = TRUE,
                plotOutput("heatmap_plot", height = "500px"))
          )
        )
      ),
      
      # Procedure Info Tab
      tabItem(
        tabName = "procedure_view",
        h2("Procedure Information Reference"),
        p("Browse the complete mapping of parameters to procedures."),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Filters"),
              selectInput("proc_category_filter", "Parameter Category:", choices = NULL),
              checkboxInput("proc_only_mandatory", "Only Mandatory", value = FALSE)
            )
          ),
          column(9,
            box(width = 12, title = "Parameter-Procedure Mapping",
                status = "primary", solidHeader = TRUE,
                DTOutput("procedure_table")),
            box(width = 12, title = "Statistics",
                status = "info", solidHeader = TRUE,
                verbatimTextOutput("procedure_stats"))
          )
        )
      ),
      
      # Disease Associations Tab
      tabItem(
        tabName = "disease_view",
        h2("Disease Associations"),
        p("Explore disease associations for your query genes."),
        
        fluidRow(
          column(3,
            wellPanel(
              h4("Gene Selection"),
              selectInput("disease_gene_select", "Select Gene:", choices = NULL),
              hr(),
              h4("Filter Options"),
              checkboxInput("disease_show_all", "Show all associations", value = TRUE)
            )
          ),
          column(9,
            box(width = 12, title = "Disease Information",
                status = "primary", solidHeader = TRUE,
                DTOutput("disease_table")),
            box(width = 12, title = "Summary Statistics",
                status = "info", solidHeader = TRUE,
                uiOutput("disease_stats"))
          )
        )
      ),
      
      # About Tab 
      tabItem(
        tabName = "about",
        h2("About This Dashboard"),
        box(width = 12,
            h3("IMPC Phenotype Analysis Dashboard - MASTER VERSION"),
            p("Interactive exploration of IMPC mouse phenotype data with comprehensive features."),
            
            h4("Features:"),
            tags$ul(
              tags$li(strong("FDR Correction"), " - Benjamini-Hochberg multiple testing correction"),
              tags$li(strong("Bubble Chart"), " - Novel Visualisation with transformed p-values"),
              tags$li(strong("Disease Associations"), " - Link genes to disease information"),
              tags$li(strong("Enhanced Clustering"), " - Shows most significant parameter per gene"),
              tags$li(strong("Procedure Information"), " - Tooltips show procedure names and descriptions"),
              tags$li(strong("Mandatory Filter"), " - Filter by mandatory/optional procedures"),
              tags$li(strong("Flexible Data Source"), " - Works with CSV files or MySQL database"),
              tags$li(strong("Fisher's Method"), " - P-value combination across strains/life stages"),
              tags$li(strong("Custom Gene Lists"), " - Paste your own gene symbols"),
              tags$li(strong("Comprehensive Filtering"), " - Category, strain, life stage filters")
            ),
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
            
            h4("Data Statistics:"),
            uiOutput("data_stats")
        )
      )
    )
  )
)


# SERVER LOGIC

server <- function(input, output, session) {
  
  # Load data on startup using flexible loader
  data <- reactive({
    df <- load_data()
    
    # Add FDR correction
    cat("  Applying FDR correction (Benjamini-Hochberg)...\n")
    df <- df %>%
      mutate(p_adj = p.adjust(pvalue, method = "BH"),
             sig_fdr = ifelse(!is.na(p_adj) & p_adj < 0.05, TRUE, FALSE))
    
    cat("  FDR correction complete!\n")
    return(df)
  })
  
  # Load disease information
  disease_data <- reactive({
    tryCatch({
      disease_info <- read.csv(CSV_FILES$disease, header = TRUE, stringsAsFactors = FALSE)
      cat("  Loaded disease information\n")
      return(disease_info)
    }, error = function(e) {
      cat("  Warning: Could not load disease information:", e$message, "\n")
      return(NULL)
    })
  })
  
  # Display data source in UI
  output$data_source_indicator <- renderText({
    toupper(DATA_SOURCE)
  })
  
  # Display data statistics in About tab
  output$data_stats <- renderUI({
    df <- data()
    proc_count <- sum(!is.na(df$procedure_name))
    proc_pct <- round(proc_count / nrow(df) * 100, 1)
    fdr_sig_count <- sum(df$sig_fdr, na.rm = TRUE)
    fdr_sig_pct <- round(fdr_sig_count / nrow(df) * 100, 1)
    
    tags$div(
      tags$p(strong("Total Records: "), format(nrow(df), big.mark = ",")),
      tags$p(strong("Unique Genes: "), n_distinct(df$gene_symbol)),
      tags$p(strong("Unique Parameters: "), n_distinct(df$parameter_id)),
      tags$p(strong("Significant (FDR < 0.05): "), 
             format(fdr_sig_count, big.mark = ","), " (", fdr_sig_pct, "%)"),
      tags$p(strong("Parameters with Procedure Info: "), 
             format(proc_count, big.mark = ","), " (", proc_pct, "%)"),
      tags$p(strong("Mouse Strains: "), paste(levels(df$mouse_strain), collapse = ", ")),
      tags$p(strong("Life Stages: "), paste(levels(df$mouse_life_stage), collapse = ", "))
    )
  })
  
  # Initialise dropdown choices
  observe({
    df <- data()
    
    gene_choices <- c("All" = "all", sort(unique(df$gene_combined)))
    updateSelectInput(session, "gene_select", choices = gene_choices)
    updateSelectInput(session, "bubble_gene_select", choices = gene_choices)
    updateSelectInput(session, "disease_gene_select", choices = sort(unique(df$gene_symbol)))
    
    category_choices <- c("All" = "all", sort(unique(as.character(df$category))))
    updateSelectInput(session, "gene_category_filter", choices = category_choices)
    updateSelectInput(session, "bubble_category_filter", choices = category_choices)
    updateSelectInput(session, "phenotype_category_filter", choices = category_choices)
    updateSelectInput(session, "cluster_category_filter", choices = category_choices)
    updateSelectInput(session, "proc_category_filter", choices = category_choices)
    
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
    # Return NULL if textbox is empty or only whitespace
    if (is.null(input$custom_genes) || trimws(input$custom_genes) == "") {
      return(NULL)
    }
    # Process the gene list
    genes <- unlist(strsplit(input$custom_genes, "\n"))
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    
    # Return NULL if no valid genes after processing
    if (length(genes) == 0) {
      return(NULL)
    }
    return(genes)
  })
  
  # Filter data for gene view
  gene_data <- reactive({
    df <- data()
    
    # Apply filters first
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
    if (input$gene_only_mandatory) {
      df <- df[!is.na(df$is_mandatory) & df$is_mandatory == TRUE, ]
    }
    
    # Create procedure lookup BEFORE aggregation
    proc_lookup <- df %>%
      select(parameter_id, parameter_name, procedure_name, procedure_description, is_mandatory) %>%
      distinct()
    
    # Combine p-values using Fisher's method
    df_combined <- df %>%
      group_by(gene_symbol, parameter_id, parameter_name, category) %>%
      summarise(pvalue = combine_pvalue(pvalue),
               p_adj = combine_pvalue(p_adj), .groups = 'drop')
    
    # Add back procedure information
    df_combined <- df_combined %>%
      left_join(proc_lookup, by = c("parameter_id", "parameter_name"))
    
    # Use appropriate p-value
    if (input$gene_use_fdr) {
      df_combined$plot_pvalue <- df_combined$p_adj
    } else {
      df_combined$plot_pvalue <- df_combined$pvalue
    }
    
    df_combined$neg_log10_pvalue <- -log10(df_combined$plot_pvalue)
    
    return(df_combined)
  })
  
  # Gene view plot
  output$gene_plot <- renderPlotly({
    df <- gene_data()
    
    if (nrow(df) == 0) {
      return(plot_ly() %>%
        add_annotations(text = "No data available for selected filters.", 
                       showarrow = FALSE, font = list(size = 16)))
    }
    
    df$significant <- ifelse(df$plot_pvalue <= input$gene_pvalue_threshold, "Significant", "Not Significant")
    
    # Create enhanced tooltips with procedure information
    df$tooltip <- apply(df, 1, function(row) {
      pval_label <- if (input$gene_use_fdr) "FDR-adjusted p-value" else "P-value"
      tooltip_text <- paste0(
        "<b>Parameter:</b> ", row["parameter_name"],
        "<br><b>Category:</b> ", row["category"],
        "<br><b>", pval_label, ":</b> ", format(as.numeric(row["plot_pvalue"]), scientific = TRUE, digits = 3)
      )
      
      if (!is.na(row["procedure_name"]) && row["procedure_name"] != "") {
        tooltip_text <- paste0(tooltip_text,
                              "<br><b>Procedure:</b> ", row["procedure_name"])
      }
      
      if (!is.na(row["is_mandatory"]) && row["is_mandatory"] != "") {
        mandatory_text <- ifelse(row["is_mandatory"] == "TRUE", "Yes", "No")
        tooltip_text <- paste0(tooltip_text,
                              "<br><b>Mandatory:</b> ", mandatory_text)
      }
      
      tooltip_text
    })
    
    plot_ly(df, x = ~seq_along(parameter_name), y = ~neg_log10_pvalue,
            type = 'scatter', mode = 'markers',
            marker = list(size = 8, color = ~plot_pvalue, colorscale = 'Plasma',
                         colorbar = list(title = if (input$gene_use_fdr) "FDR p-value" else "P-value"),
                         line = list(color = ~ifelse(significant == "Significant", "red", "white"), width = 2)),
            text = ~tooltip,
            hoverinfo = 'text') %>%
      layout(title = paste("Phenotypes for", if (is.null(custom_genes())) input$gene_select else "Custom Gene List"),
             xaxis = list(title = "Parameters (ordered by category)", showticklabels = FALSE),
             yaxis = list(title = paste("-log10(", if (input$gene_use_fdr) "FDR-adjusted " else "", "p-value)", sep = "")),
             shapes = list(list(type = "line",
                               y0 = -log10(input$gene_pvalue_threshold),
                               y1 = -log10(input$gene_pvalue_threshold),
                               x0 = 0, x1 = nrow(df),
                               line = list(color = "red", dash = "dash", width = 2))))
  })
  
  # Gene view table
  output$gene_table <- renderDT({
    df <- gene_data()
    df_sig <- df[df$plot_pvalue <= input$gene_pvalue_threshold, ]
    
    if (nrow(df_sig) == 0) {
      return(datatable(data.frame(Message = "No significant phenotypes found"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    pval_col_name <- if (input$gene_use_fdr) "FDR p-value" else "P-value"
    summary_df <- df_sig %>%
      arrange(plot_pvalue) %>%
      select(Parameter = parameter_name, Category = category,
             Procedure = procedure_name, Mandatory = is_mandatory) %>%
      mutate(!!pval_col_name := format(df_sig$plot_pvalue[order(df_sig$plot_pvalue)], 
                                       scientific = TRUE, digits = 3))
    
    datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE,
                                         order = list(list(4, 'asc'))), rownames = FALSE)
  })
  
  # Bubble Chart
  output$bubble_chart <- renderPlotly({
    df <- data()
    
    # Apply gene filter
    if (input$bubble_gene_select != "all") {
      df <- df[df$gene_combined == input$bubble_gene_select, ]
    }
    
    # Apply category filter
    if (input$bubble_category_filter != "all") {
      df <- df[df$category == input$bubble_category_filter, ]
    }
    
    if (nrow(df) == 0) {
      return(plot_ly() %>%
        add_annotations(text = "No data available for selected gene.", 
                       showarrow = FALSE, font = list(size = 16)))
    }
    
    # Use appropriate p-value
    if (input$bubble_use_fdr) {
      df$plot_pvalue <- df$p_adj
    } else {
      df$plot_pvalue <- df$pvalue
    }
    
    # Aggregate by parameter (take minimum p-value)
    df_agg <- df %>%
      group_by(parameter_name, category) %>%
      summarise(p_value = min(plot_pvalue, na.rm = TRUE), .groups = "drop") %>%
      mutate(ptrans = (1 - p_value)^8,  # Transform p-value for Visualisation
             sig_group = ifelse(p_value < input$bubble_pvalue_threshold, "Significant", "Non-significant"))
    
    if (nrow(df_agg) == 0) return(NULL)
    
    plot_ly(df_agg,
            x = ~parameter_name,
            y = ~ptrans,
            type = "scatter",
            mode = "markers",
            marker = list(
              size = ~ptrans * 50,  # Bubble size proportional to transformed p-value
              color = ~sig_group,
              colors = c("Significant" = "steelblue", "Non-significant" = "grey"),
              line = list(color = "black", width = 1)
            ),
            text = ~paste0(
              "<b>Parameter:</b> ", parameter_name,
              "<br><b>Category:</b> ", category,
              "<br><b>", if (input$bubble_use_fdr) "FDR p-value" else "P-value", ":</b> ", 
              signif(p_value, 3),
              "<br><b>Significance:</b> ", sig_group
            ),
            hoverinfo = "text") %>%
      layout(
        title = paste("Bubble Chart for", input$bubble_gene_select),
        xaxis = list(title = "Parameters", tickangle = -45, automargin = TRUE),
        yaxis = list(title = "Transformed p-value [(1-p)^8]"),
        showlegend = TRUE
      )
  })
  
  # Filter data for phenotype view
  phenotype_data <- reactive({
    df <- data()
    
    if (input$parameter_select == "all") {
      return(df)
    }
    
    df <- df[df$parameter_combined == input$parameter_select, ]
    
    # Use appropriate p-value
    if (input$phenotype_use_fdr) {
      df$plot_pvalue <- df$p_adj
    } else {
      df$plot_pvalue <- df$pvalue
    }
    
    df$neg_log10_pvalue <- -log10(df$plot_pvalue)
    
    return(df)
  })
  
  # Phenotype view plot
  output$phenotype_plot <- renderPlotly({
    df <- phenotype_data()
    
    if (nrow(df) == 0 || input$parameter_select == "all") {
      return(plot_ly() %>%
        add_annotations(text = "Please select a specific parameter.", 
                       showarrow = FALSE, font = list(size = 16)))
    }
    
    df$significant <- ifelse(df$plot_pvalue <= input$phenotype_pvalue_threshold, "Significant", "Not Significant")
    
    # Create enhanced tooltips
    df$tooltip <- apply(df, 1, function(row) {
      pval_label <- if (input$phenotype_use_fdr) "FDR-adjusted p-value" else "P-value"
      tooltip_text <- paste0(
        "<b>Gene:</b> ", row["gene_symbol"],
        "<br><b>MGI ID:</b> ", row["gene_accession_id"],
        "<br><b>", pval_label, ":</b> ", format(as.numeric(row["plot_pvalue"]), scientific = TRUE, digits = 3)
      )
      
      if (!is.na(row["procedure_name"]) && row["procedure_name"] != "") {
        tooltip_text <- paste0(tooltip_text,
                              "<br><b>Procedure:</b> ", row["procedure_name"])
      }
      
      tooltip_text
    })
    
    plot_ly(df, x = ~seq_along(gene_symbol), y = ~neg_log10_pvalue,
            type = 'scatter', mode = 'markers',
            marker = list(size = 8, color = ~plot_pvalue, colorscale = 'Plasma',
                         colorbar = list(title = if (input$phenotype_use_fdr) "FDR p-value" else "P-value"),
                         line = list(color = ~ifelse(significant == "Significant", "red", "white"), width = 2)),
            text = ~tooltip,
            hoverinfo = 'text') %>%
      layout(title = "Gene Significance for Selected Phenotype(s)",
             xaxis = list(title = "Genes", showticklabels = FALSE),
             yaxis = list(title = paste("-log10(", if (input$phenotype_use_fdr) "FDR-adjusted " else "", "p-value)", sep = "")),
             shapes = list(list(type = "line",
                               y0 = -log10(input$phenotype_pvalue_threshold),
                               y1 = -log10(input$phenotype_pvalue_threshold),
                               x0 = 0, x1 = nrow(df),
                               line = list(color = "red", dash = "dash", width = 2))))
  })
  
  # Phenotype view table
  output$phenotype_table <- renderDT({
    df <- phenotype_data()
    df_sig <- df[df$plot_pvalue <= input$phenotype_pvalue_threshold, ]
    
    if (nrow(df_sig) == 0) {
      return(datatable(data.frame(Message = "No significant genes found"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    pval_col_name <- if (input$phenotype_use_fdr) "FDR p-value" else "P-value"
    summary_df <- df_sig %>%
      arrange(plot_pvalue) %>%
      select(Gene = gene_symbol, `MGI ID` = gene_accession_id,
             Strain = mouse_strain, Procedure = procedure_name) %>%
      mutate(!!pval_col_name := format(df_sig$plot_pvalue[order(df_sig$plot_pvalue)], 
                                       scientific = TRUE, digits = 3))
    
    datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE,
                                         order = list(list(4, 'asc'))), rownames = FALSE)
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
    if (input$cluster_only_mandatory) {
      df <- df[!is.na(df$is_mandatory) & df$is_mandatory == TRUE, ]
    }
    
    # Use appropriate p-value
    if (input$cluster_use_fdr) {
      df$plot_pvalue <- df$p_adj
    } else {
      df$plot_pvalue <- df$pvalue
    }
    
    df <- df[df$plot_pvalue <= input$cluster_pvalue_threshold, ]
    
    if (nrow(df) == 0) return(NULL)
    
    df_combined <- df %>%
      group_by(gene_symbol, parameter_id, parameter_name) %>%
      summarise(pvalue = combine_pvalue(plot_pvalue), .groups = 'drop')
    
    mat <- df_combined %>%
      mutate(significant = 1) %>%
      pivot_wider(names_from = parameter_id, values_from = significant, values_fill = 0) %>%
      as.data.frame()
    
    rownames(mat) <- mat$gene_symbol
    mat <- mat[, -1]
    
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    
    # Get most significant parameter for each gene
    most_sig_param <- df_combined %>%
      group_by(gene_symbol) %>%
      slice_min(order_by = pvalue, n = 1) %>%
      select(gene_symbol, most_sig_parameter = parameter_name)
    
    list(matrix = mat, genes = rownames(mat), raw_data = df_combined, 
         most_sig = most_sig_param)
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
    
    # Add most significant parameter info
    plot_df <- plot_df %>%
      left_join(clust$most_sig, by = c("gene" = "gene_symbol"))
    
    # Enhanced tooltip with most significant parameter
    plot_ly(plot_df, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
            marker = list(size = ~sig_count * 2, color = ~as.numeric(cluster),
                         colorscale = 'Viridis', colorbar = list(title = "Cluster"),
                         line = list(color = "white", width = 1)),
            text = ~paste("<b>Gene:</b>", gene, 
                         "<br><b>Cluster:</b>", cluster, 
                         "<br><b>Significant Phenotypes:</b>", sig_count,
                         "<br><b>Most Significant Parameter:</b>", most_sig_parameter),
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
  
  # Procedure Info Tab - Table
  output$procedure_table <- renderDT({
    df <- data()
    
    # Get unique parameter-procedure mappings
    proc_df <- df %>%
      select(parameter_id, parameter_name, category, procedure_name, 
             procedure_description, is_mandatory) %>%
      distinct()
    
    # Apply filters
    if (input$proc_category_filter != "all") {
      proc_df <- proc_df[proc_df$category == input$proc_category_filter, ]
    }
    if (input$proc_only_mandatory) {
      proc_df <- proc_df[!is.na(proc_df$is_mandatory) & proc_df$is_mandatory == TRUE, ]
    }
    
    # Clean up for display
    proc_df <- proc_df %>%
      arrange(category, parameter_name) %>%
      select(Category = category, Parameter = parameter_name, 
             `Procedure Name` = procedure_name, 
             `Procedure Description` = procedure_description,
             Mandatory = is_mandatory)
    
    datatable(proc_df, options = list(pageLength = 15, scrollX = TRUE), 
              rownames = FALSE, filter = 'top')
  })
  
  # Procedure Info Tab - Statistics
  output$procedure_stats <- renderText({
    df <- data()
    
    total_params <- n_distinct(df$parameter_id)
    params_with_proc <- n_distinct(df$parameter_id[!is.na(df$procedure_name)])
    mandatory_count <- sum(!is.na(df$is_mandatory) & df$is_mandatory == TRUE)
    
    paste0(
      "Total unique parameters: ", total_params, "\n",
      "Parameters with procedure info: ", params_with_proc, 
      " (", round(params_with_proc/total_params*100, 1), "%)\n",
      "Mandatory procedure entries: ", mandatory_count, "\n"
    )
  })
  
  # Disease Associations
  output$disease_table <- renderDT({
    disease_info <- disease_data()
    
    if (is.null(disease_info)) {
      return(datatable(data.frame(Message = "Disease information not available"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    if (input$disease_gene_select == "") {
      return(datatable(data.frame(Message = "Please select a gene"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    # Get gene accession ID
    df <- data()
    gene_info <- df %>%
      filter(gene_symbol == input$disease_gene_select) %>%
      select(gene_accession_id) %>%
      distinct()
    
    if (nrow(gene_info) == 0) {
      return(datatable(data.frame(Message = "Gene not found in dataset"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    # Filter disease info
    disease_filtered <- disease_info %>%
      filter(gene_accession_id %in% gene_info$gene_accession_id)
    
    if (nrow(disease_filtered) == 0) {
      return(datatable(data.frame(Message = "No disease associations found for this gene"),
                      options = list(dom = 't'), rownames = FALSE))
    }
    
    # Format for display
    disease_display <- disease_filtered %>%
      select(`Gene` = gene_accession_id, everything())
    
    datatable(disease_display, options = list(pageLength = 10, scrollX = TRUE), 
              rownames = FALSE)
  })
  
  # Disease statistics
  output$disease_stats <- renderUI({
    disease_info <- disease_data()
    
    if (is.null(disease_info)) {
      return(tags$p("Disease information not available"))
    }
    
    if (input$disease_gene_select == "") {
      return(tags$p("Select a gene to view statistics"))
    }
    
    df <- data()
    gene_info <- df %>%
      filter(gene_symbol == input$disease_gene_select) %>%
      select(gene_accession_id) %>%
      distinct()
    
    disease_filtered <- disease_info %>%
      filter(gene_accession_id %in% gene_info$gene_accession_id)
    
    tags$div(
      tags$p(strong("Gene: "), input$disease_gene_select),
      tags$p(strong("Number of disease associations: "), nrow(disease_filtered))
    )
  })
}


# RUN APPLICATION

shinyApp(ui, server)