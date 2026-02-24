library(shiny)
library(Seurat)
library(tibble)
library(dplyr)
library(tidyr)
library(shinydashboard)
library(ggthemes)
library(gt)
library(DT)
library(shinyjs)
library(ggprism)
library(shinycssloaders)
library(ggplot2)

## Load data ----

message("Loading data...")
sc <- readRDS("data/EMP85_HematoEC_APP.rds")
message("Data loaded successfully")

message("Loading gene list, ensembl2symbol, subcellular location, cell percentage and counts data...")

ensembl_symbol <- readRDS("data/symbol2ensembl_2026_02_16.RDS")
gene_list <- readRDS("data/gene_list.RDS") 
subcellular_locations <- readRDS("data/Subcellular_locations.RDS")
cell_pct_genes <- readRDS("data/Cell_PCT_Genes_LONG.rds")
cell_count_genes <- readRDS("data/Cell_Counts_Genes.rds")

cluster_list <- readRDS("data/cluster_list.RDS")
allMarkers <- readRDS("data/all_markers.RDS")


## Define functions ----

getFeaturePlot <- function(dataset, ensembl_symbol_df, symbol_in, color, ...){
  
  ## Debugging
  # dataset = sc
  # symbol_in = "Kdr"
  # ensembl_symbol_df = ensembl_symbol
  # color = "firebrick2"
  
  # Ensembl to symbol
  
  ensembl_id <- ensembl_symbol_df %>% 
    dplyr::filter(symbol == symbol_in) %>% 
    dplyr::pull(ensembl) %>%
    as.character() %>% 
    unique()
  
  plot <- FeaturePlot(object = dataset, features = ensembl_id, slot = "data")
  plot$labels$title <- symbol_in
  plot$labels$x <- "UMAP_1"
  plot$labels$y <- "UMAP_2"
  
  if(missing(color) == FALSE){
    plot <- plot + 
      scale_colour_gradientn(colours = c("grey", color))
  }
  
  if(all(plot[[1]][["data"]][[4]] == 0)){
    plot <- plot + scale_colour_gradientn(colours = c("grey", "grey"))
  }
  
  return(plot)
}



getVln <- function(dataset, ensembl, symbol, celltype_levels_in){
  
  VlnPlot(object = dataset, features = ensembl) + 
    ggtitle(symbol) + 
    NoLegend() +
    theme(axis.title.x = element_blank()) +
    theme(plot.margin = unit(c(0, 0, 0, 1), "inches")) 
}



# Pct cells expressing gene per cluster
getPct <- function(symbol_in, ensembl_symbol_in, cell_count_genes_in, pct_df_in){
  
  # Fetch ensemblid
  ensid <- ensembl_symbol_in %>% 
    dplyr::filter(symbol == symbol_in) %>% 
    dplyr::pull(ensembl)
  
  # Get counts for query gene
  geneCellCounts <- cell_count_genes_in %>% 
    dplyr::select(Cluster, ensid, n)
  
  # Get percentages for query gene
  pct_df_counts <- pct_df_in %>% 
    dplyr::filter(symbol == symbol_in) %>% 
    dplyr::left_join(geneCellCounts, by = "Cluster") %>% 
    rename_at(vars(starts_with("ENSMUSG")),  ~"counts") %>% 
    dplyr::mutate(Cluster = factor(Cluster, levels = levels(Idents(sc))))
  
  # Plot
  ggplot(pct_df_counts, aes(x = Cluster, y = pct, fill = Cluster)) +
    geom_bar(stat="identity") +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 120)) +
    geom_text(aes(label = paste0(counts, " / ", n)), angle = 90, hjust = -0.1, fontface = "bold", family = "sans") +
    ggtitle(symbol_in) +
    ylab(paste0("Percentage of cells expressing ", symbol_in)) +
    NoLegend() +
    theme_light() +
    theme(panel.grid.major.y = element_line(color = "grey", size = 0.5, linetype = 2)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12, colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(family = "sans", size = 14),
          plot.margin = unit(c(0, 0, 0, 1), "inches")) +
    NoLegend()
}


# Custom markers

getCustomMarkers <- function(dataset, groupA, groupB, ensembl_symbol_in){
  
  # Create metadata column with custom groups
  newMetadata <- dataset@meta.data %>% 
    dplyr::mutate(custom_cluster = case_when(Cluster %in% groupA ~ "group_a",
                                             Cluster %in% groupB ~ "group_b",
                                             .default = Cluster))
  dataset <- Seurat::AddMetaData(dataset, newMetadata)
  
  # Set identities and calculate markers
  Idents(dataset) <- dataset$custom_cluster
  
  markers <- FindMarkers(dataset, ident.1 = "group_a", ident.2 = "group_b")
  
  markers <- markers %>% 
    rownames_to_column("ensembl") %>% 
    dplyr::left_join(ensembl_symbol_in, by = "ensembl") %>% 
    dplyr::relocate(symbol) %>% 
    dplyr::mutate(p_val = format(p_val, digits = 3)) %>% 
    dplyr::mutate(avg_log2FC = round(avg_log2FC, 2)) %>% 
    dplyr::mutate(p_val_adj = format(p_val_adj, digits = 3))
  
  return(markers)
}



get1vsAllMarkers <- function(df, cluster_in){
  
  df_subset <- df %>% 
    dplyr::filter(cluster == cluster_in) %>% 
    dplyr::select(-cluster)
  
  return(df_subset)
}



getSubcellularLocation <- function(df, ensembl_in){
  df <- df %>% 
    dplyr::filter(ENSEMBL_ID == ensembl_in)
}




################################################################################

ui <- dashboardPage(
  
  dashboardHeader(title = "E8.5 HaematoEC"),
  
  # Sidebar
  dashboardSidebar(
    
    sidebarMenu( id = "sidebarMenu",
                 
                 ### User inputs
                 menuItem("Inputs", tabName = "inputs", icon = icon("keyboard")),
                 
                 selectizeInput(inputId = "geneInput",
                                label = "Select a gene",
                                choices = NULL,
                                selected = character(0),
                                multiple = FALSE,
                                options = list(placeholder="None selected...")),
                 actionButton(inputId = "geneInput_button", "Run gene"),
                 
                 hr(id = "hr1"),
                 
                 selectizeInput(inputId = "clusterInputA",
                                label = "Compare these cluster/s:",
                                choices = NULL,
                                multiple = TRUE,
                                options = list(placeholder="None selected...")),
                 
                 selectizeInput(inputId = "clusterInputB",
                                label = "Versus these cluster/s",
                                choices = NULL,
                                multiple = TRUE,
                                options = list(placeholder="None selected...")),
                 
                 actionButton(inputId = "clusterInput_buton", "Calculate markers"),
                 
                 hr(id = "hr2"),
                 
                 selectizeInput(inputId = "clusterInput1vsAll",
                                label = "Select a cluster to compare against the rest of clusters",
                                choices = NULL,
                                multiple = FALSE,
                                options = list(placeholder="None selected...")),
                 
                 actionButton(inputId = "clusterInput1vsAll_buton", "Calculate markers"),
                 
                 # Output tabs
                 menuItem(h2("Outputs:")),
                 menuItem("Gene expression", tabName = "expr", icon = icon("image"),
                          badgeLabel = "Gene"),
                 menuItem("Cluster markers", tabName = "custom_markers", icon = icon("table"),
                          badgeLabel = "Custom"),
                 menuItem("Cluster markers", tabName = "1vsall_markers", icon = icon("table"),
                          badgeLabel = "1vs.All")
    )
  ),
  
  ### Body
  dashboardBody(
    
    # ### style.css needs to be in /www folder; this command needs to be in body
    # tags$head(
    #   tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    # ),
    
    # Expression UI tab
    tabItems(
      tabItem(tabName = "expr",
              fluidPage(
                fluidRow(h2(htmlOutput(outputId = "geneSymbol_text"))),
                
                # Box with reference UMAPs
                fluidRow(
                  box(title = 'Gene expression', width = 12, solidHeader = T, status = 'primary',
                      column(width = 7, img(src = "umap_NUMBERS.png", style="width: 100%")),
                      column(width = 5, plotOutput(outputId = "featurePlot", height = "auto") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                  )
                ),
                
                # Box with violin plots
                fluidRow(
                  box(title = "Gene expression: Violin plot", width = 8, solidHeader = T, status = 'primary',
                      column(width = 12, plotOutput(outputId = "vlnPlot", height = "auto") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                  )
                ),
                
                # Box with cell cluster percentages
                fluidRow(
                  box(title = '% Cells expressing gene per cluster', width = 8, solidHeader = T, status = 'primary',
                      column(width = 12, plotOutput(outputId = "pctPlot", height = "auto") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                  )
                ),
                
                # Box with Protein subcellular location
                fluidRow(
                  box(title = "Subcellular localisation of encoded protein", width = 8, solidHeader = T, status = "primary",
                      column(width = 12, DT::dataTableOutput(outputId = "subcellularLocation", width = "60%") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                      
                  )
                )
                
                
              ),
              
      ),
      
      ### Custom markers UI tab
      tabItem(tabName = "custom_markers",
              fluidPage(
                # Title and reminder of compared clusters
                fluidRow(
                  h2("Marker results"),
                  htmlOutput(outputId = "clusterInputA_text"),
                  htmlOutput(outputId = "clusterInputB_text"),
                  HTML("<br>")
                ),
                # Results
                fluidRow(
                  box(width = "100%",
                      column(width = 12, DT::dataTableOutput(outputId = "customClustersTable", width = "60%") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                  )
                )
              )
      ),
      
      ### 1vsAll markers UI tab
      tabItem(tabName = "1vsall_markers",
              fluidPage(
                # Title and reminder of compared clusters
                fluidRow(
                  h2("Marker results"),
                  htmlOutput(outputId = "clusterInput1vsAll_text"),
                  HTML("<br>")
                ),
                # Results
                fluidRow(
                  box(width = "100%",
                      column(width = 12, DT::dataTableOutput(outputId = "clusters1vsAllTable", width = "60%") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
                  )
                )
              )
      )
      
    )
    
  )
) # clusters1vsAllTable

################################################################################

server <- function(input, output, session) { 
  
  ### Automatic tab switching
  observeEvent(input$geneInput_button, {
    req(input$geneInput)
    updateTabsetPanel(session, "sidebarMenu",
                      selected = "expr")
  })
  
  observeEvent(input$clusterInput_buton, {
    req(input$clusterInputB)
    updateTabsetPanel(session, "sidebarMenu",
                      selected = "custom_markers")
  })
  
  observeEvent(input$clusterInput1vsAll_buton, {
    req(input$clusterInput1vsAll)
    updateTabsetPanel(session, "sidebarMenu",
                      selected = "1vsall_markers")
  })
  
  ### Display possible inputs in dropdown menus
  updateSelectizeInput(session, "geneInput", choices = gene_list, server = TRUE)
  updateSelectizeInput(session, "clusterInputA", choices = cluster_list, server = TRUE)
  updateSelectizeInput(session, "clusterInputB", choices = cluster_list, server = TRUE)
  updateSelectizeInput(session, "clusterInput1vsAll", choices = cluster_list, server = TRUE)
  
  ### Change tabs
  
  ### Process inputs
  # Gene symbol
  geneSymbol_query <- eventReactive(input$geneInput_button, {
    input_symbol <- input$geneInput
  })
  output$geneSymbol_text <- eventReactive(input$geneInput_button, {
    paste0("Gene expression results: <b>", geneSymbol_query(), "<b>")
  })
  # Gene ensembl
  geneEnsembl_query <- reactive({
    input_ensembl <- ensembl_symbol %>% 
      dplyr::filter(symbol == geneSymbol_query()) %>% 
      dplyr::pull(ensembl)
  })
  
  # Clusters
  clusterInputA_query <- eventReactive(input$clusterInput_buton, {
    input$clusterInputA
  })
  output$clusterInputA_text <- eventReactive(input$clusterInput_buton, {
    paste(c("Comparing clusters: <b>", paste(input$clusterInputA, collapse = ", ")), "<b>", collapse = " ")
  })
  
  clusterInputB_query <- eventReactive(input$clusterInput_buton, {
    input$clusterInputB
  })
  output$clusterInputB_text <- eventReactive(input$clusterInput_buton, {
    paste(c("Versus clusters: <b>", paste(input$clusterInputB, collapse = ", ")), "<b>", collapse = " ")
  })
  
  clusterInput1vsAll_query <- eventReactive(input$clusterInput1vsAll_buton, {
    input$clusterInput1vsAll
  })
  output$clusterInput1vsAll_text <- eventReactive(input$clusterInput1vsAll_buton, {
    paste0("Markers of cluster: <b>", input$clusterInput1vsAll, "<b> vs. rest of clusters")
  })
  
  ### Gene expression FeaturePlots
  
  # Need to set height = "auto" in plotOutput() and then set height like below
  # to keep img aspect ratio despite window resizing
  output$featurePlot <- renderPlot({
    getFeaturePlot(dataset = sc, ensembl = geneEnsembl_query(), symbol_in = geneSymbol_query(), ensembl_symbol_df = ensembl_symbol, color = "firebrick2")
  }, height = function() {session$clientData$output_featurePlot_width * 1.0})
  
  ### Violin plot
  output$vlnPlot <- renderPlot({
    getVln(dataset = sc, ensembl = geneEnsembl_query(), symbol = geneSymbol_query())
  }, height = function() {session$clientData$output_pctPlot_width * 0.35})
  
  ### Gene expression pct cells by celltype
  output$pctPlot <- renderPlot({
    getPct(symbol_in = geneSymbol_query(), 
           ensembl_symbol_in = ensembl_symbol, 
           cell_count_genes_in = cell_count_genes,
           pct_df_in = cell_pct_genes)
  }, height = function() {session$clientData$output_pctPlot_width * 0.35})
  
  ### Subcellular location
  output$subcellularLocation <- DT::renderDataTable({
    getSubcellularLocation(df = subcellular_locations, ensembl_in = geneEnsembl_query())
  }, rownames = TRUE, options = list(pageLength = 50))
  
  ### Custom clusters
  output$customClustersTable <- DT::renderDataTable({
    getCustomMarkers(dataset = sc, groupA = clusterInputA_query(), groupB = clusterInputB_query(), ensembl_symbol_in = ensembl_symbol)
  }, rownames = TRUE, options = list(pageLength = 50))
  
  ### 1vsAll markers
  output$clusters1vsAllTable <- DT::renderDataTable({
    get1vsAllMarkers(df = allMarkers, cluster_in = clusterInput1vsAll_query())
  }, rownames = TRUE, options = list(pageLength = 50))
}

################################################################################

shinyApp(ui, server)