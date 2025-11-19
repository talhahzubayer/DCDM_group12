library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(umap)

# Load your datasets
phenotype <- read.csv("C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/clean_table_final_UPPERCASE.csv", header = TRUE, stringsAsFactors = FALSE)
disease_info <- read.csv("C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/Disease_information_cleaned.csv", header = TRUE, stringsAsFactors = FALSE)
parameter_description <- read.csv("C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_parameter_description_cleaned.csv", header = TRUE, stringsAsFactors = FALSE)
procedure <- read.csv("C:/Users/Talhah Zubayer/Documents/DCDM_group12/data/processed/IMPC_procedure_cleaned.csv", header = TRUE, stringsAsFactors = FALSE)

# Merge datasets
# cuz the information from multiple tables therefore need to combine it 
#inner_join: with paramter_description: ensure only keep the records that have valid parameter information.
# > it is essential to match both paramter_id and parameter_name for accuracy 
#left_join: with procedure: add procedure metadata to each parameter, but keeps all existing phenotype records even if there is no matching procedure 
#left_join: with diease_info: links genes to disease associations, many-to-many is used because a gene can be linked to mulitple disease, and a disease can be linked to multiple genes.
# > this step enriches the dataset with biological context
dataset <- phenotype %>%
  inner_join(parameter_description, by = c("parameter_id" = "parameterId", "parameter_name" = "name")) %>%
  left_join(procedure, by = "impcParameterOrigId") %>%
  left_join(disease_info, by = c("gene_accession_id"), relationship = "many-to-many")

# Ensure pvalue is numeric
#sometime p-values are stpred as characters whrn read from csv file 
#converting them to numberic ensures that can filter, compare, and adjust p-value later 
dataset <- dataset %>% mutate(pvalue = as.numeric(pvalue))

# Set alpha threshold
#define the cutoff for statistical significance 
#p-value below 0.05 are consifered significant 
#it is standard practice in hypothesis testing and helps identify wchih associations are likely meaningful 
alpha <- 0.05

# Add FDR
#when testing many paramters for many genes, some p-value will appear significant just by chance 
#FDR adjustment using Benjamini-Hochberg(BH corrects for multiple comparisons 
#FDR = 0.05 -> among all results called significant, on average 5% might be false positives 
#BH: ranks all p-value and adjust them so that the chance of false discoveries stays below a desired threshold 
# it is widely ised on genomics and high-throughput experiments where thousands of tests are performed 
#sig_fdr is a logical column marking which association remain significant after correction 
#this step ensure that downstream plots and summaries are statistically robust 
dataset <- dataset %>%
  mutate(p_adj = p.adjust(pvalue, method = "BH"),
         sig_fdr = ifelse(!is.na(p_adj) & p_adj < alpha, TRUE, FALSE))

# Only keep your four genes
#these combine columns are particularly useful for dropdown menus and hover infor in plot 
dataset <- dataset %>%
  filter(gene_symbol %in% c("Smarcd3", "Ppp3cc", "Rab12", "Klhl33")) %>%
  mutate(gene_combined = paste(gene_symbol, gene_accession_id, sep = "_"),
         parameter_combined = paste(parameter_name, parameter_id, sep = "_"))

# ----------------- Shiny UI -----------------
ui <- dashboardPage(
  dashboardHeader(title = "Simplified IMPC Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Bubble Chart", tabName = "figure1", icon = icon("circle")),
      menuItem("Scatter Plot", tabName = "figure2", icon = icon("line-chart")),
      menuItem("Clustering Plot", tabName = "figure3", icon = icon("cogs"))
    )
  ),
  dashboardBody(
    tabItems(
      # --- Bubble Chart ---"What does this gene do?"
      #reason for it 
      #visualise the significance of multiple parameters for a single gene (or knockout mouse)
      #this can see at a glance which parameter are most affected by a selected gene
      #bubble size proportional to significance (small p-values > bigger bubbles)
      #colour indicates significance vs non-significance 
      #Identifies which phenotypes are strongly associated with a specific gene knockout 
      tabItem(tabName = "figure1",
              fluidRow(
                column(4,
                       selectInput("gene_id", "Select Gene", choices = unique(dataset$gene_combined))),
                column(4,
                       sliderInput("p_value_threshold1", "p-value threshold", min = 0, max = 1, value = 0.05, step = 0.01))
              ),
              fluidRow(
                box(title = "Bubble Chart", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("bubble_chart"))
              )
      ),
      # --- Scatter Plot --- "Which genes affect this phenotype?"
      #reason for it 
      #one specific phenotype and want to see which gene significantly affect it 
      #Highlightes which genes drive changes in a specific phenotype, making it easy to spot significant associations across multiple genes 
      tabItem(tabName = "figure2",
              fluidRow(
                column(6,
                       selectInput("parameter_id", "Select Parameter", choices = unique(dataset$parameter_combined))),
                column(6,
                       sliderInput("p_value_threshold2", "p-value threshold", min = 0, max = 1, value = 0.05, step = 0.01))
              ),
              fluidRow(
                box(title = "Scatter Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("scatter_plot"))
              )
      ),
      # --- Clustering Plot --- "Which genes behave similarly across many phenotypes?"
      #reason for it
      #instead of looking at genes one by one, can cluster genes with similar phenotypes profiles 
      #UMAP reduce high-dimensional p-value data into 2D for visualisation
      #k-means clustering shows natural grouping of genes
      #Reveals patterns and relationships across multiple genes and phenotypes, which can suggest shared pathway or biological functions 
      tabItem(tabName = "figure3",
              fluidRow(
                column(4,
                       sliderInput("p_value_threshold3", "p-value threshold", min = 0, max = 1, value = 0.05, step = 0.01)),
                column(4,
                       sliderInput("k_number", "K-number (clusters)", min = 1, max = 4, value = 2, step = 1))
              ),
              fluidRow(
                box(title = "Clustering Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("clustering_plot"))
              )
      )
    )
  )
)
#Together can explore the data from multiple perspectivces: per-gene, per-phenotype, and globally across all genes 
# ----------------- Shiny Server -----------------
server <- function(input, output, session) {
  
  # --- Bubble Chart ---
  output$bubble_chart <- renderPlotly({
    df <- dataset %>%
      filter(gene_combined == input$gene_id) %>%
      group_by(parameter_name) %>%
      summarise(p_value = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
      mutate(ptrans = (1 - p_value)^8,
             sig_group = ifelse(p_value < input$p_value_threshold1, "Significant", "Non-significant"))
    
    if(nrow(df) == 0) return(NULL)
    
    plot_ly(df,
            x = ~parameter_name,
            y = ~ptrans,  # vertical placement proportional to transformed p-value
            type = "scatter",
            mode = "markers",
            marker = list(
              size = ~ptrans * 50,  # scale bubble size
              color = ~sig_group,
              colors = c("Significant" = "steelblue", "Non-significant" = "grey"),
              line = list(color = "black", width = 1)
            ),
            text = ~paste0(
              "<b>Parameter:</b> ", parameter_name,
              "<br><b>p-value:</b> ", signif(p_value, 3),
              "<br><b>Significance:</b> ", sig_group
            ),
            hoverinfo = "text") %>%
      layout(
        title = paste("Bubble Chart for", input$gene_id),
        xaxis = list(title = "Parameters", tickangle = -45, automargin = TRUE),
        yaxis = list(title = "Transformed p-value"),
        showlegend = TRUE
      )
  })
  
  # --- Scatter Plot ---
  output$scatter_plot <- renderPlotly({
    df <- dataset %>%
      filter(parameter_combined == input$parameter_id) %>%
      mutate(neg_log10_pvalue = -log10(pvalue),
             significance = ifelse(pvalue < input$p_value_threshold2, "Significant", "Not Significant"))
    
    if(nrow(df) == 0) return(NULL)
    
    plot_ly(df,
            x = ~gene_combined,
            y = ~neg_log10_pvalue,
            type = 'scatter',
            mode = 'markers+lines',
            marker = list(color = ~pvalue, colorscale = "Plasma", size = 10),
            text = ~paste("Gene:", gene_symbol, "<br>p-value:", signif(pvalue,3))) %>%
      layout(title = paste("Scatter Plot for", input$parameter_id),
             xaxis = list(title = "Genes", showticklabels = TRUE, tickangle = -45),
             yaxis = list(title = "-log10(p-value)"))
  })
  
  # --- Clustering Plot ---
  output$clustering_plot <- renderPlotly({
    mat <- dataset %>%
      group_by(gene_combined, parameter_id) %>%
      summarise(p_value = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = parameter_id, values_from = p_value, values_fill = 1)
    
    if(nrow(mat) < 2) return(NULL)
    
    mat_df <- mat %>% column_to_rownames("gene_combined")
    
    # UMAP: safe n_neighbors
    umap_res <- umap(as.matrix(mat_df), n_neighbors = min(2, nrow(mat_df)-1), random_state = 123)
    
    # k-means
    k_number <- min(input$k_number, nrow(mat_df))
    kmeans_res <- kmeans(as.matrix(mat_df), centers = k_number, nstart = 25)
    
    # Most significant parameter
    Min_Column <- apply(mat_df, 1, function(row){
      if(all(row >= input$p_value_threshold3)) "Non-significant" else names(row)[which.min(row)]
    })
    
    umap_plot <- data.frame(
      x = umap_res$layout[,1],
      y = umap_res$layout[,2],
      gene_name = rownames(mat_df),
      cluster = as.factor(kmeans_res$cluster),
      Min_Column = Min_Column
    )
    
    plot_ly(umap_plot,
            x = ~x,
            y = ~y,
            color = ~cluster,
            type = "scatter",
            mode = "markers",
            marker = list(size = 15, line = list(color = 'black', width = 1)),
            text = ~paste0(
              "<b>Gene:</b> ", gene_name,
              "<br><b>Cluster:</b> ", cluster,
              "<br><b>Most Significant Parameter:</b> ", Min_Column
            ),
            hoverinfo = "text") %>%
      layout(
        title = "UMAP Clustering of Genes",
        xaxis = list(title = "UMAP1"),
        yaxis = list(title = "UMAP2"),
        showlegend = TRUE
      )
  })
}

# Run the Shiny app
shinyApp(ui, server)
