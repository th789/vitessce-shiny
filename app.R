
# load data & packages ----------------------------------------------------

library(shiny)
library(vitessce)
library(Seurat)
source("tailored-demo-helpers.R")

#####basic demo
#load datasets
data_tcellcd8_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_results.rds")
data_pbmc_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_results.rds")
#dataset list for selection
data_list <- list(tcell_cd8="data_tcellcd8_results", pbmc="data_pbmc_results")

data_descrip_list <- list(tcell_cd8="Dataset: CD8 T-cell \n Source: Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 doi: 10.1038/ncomms14049 (2017)", 
                          pbmc="Dataset: peripheral blood mononuclear cells (PBMC) \n Source: 10x Genomics sample dataset")

data_names_list <- list(tcell_cd8="tcellcd8", pbmc="pbmc")

# descrip_tcell_cd8 = "Dataset: CD8 T-cell \n Source: Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 doi: 10.1038/ncomms14049 (2017)"
# descrip_pbmc = "Dataset: peripheral blood mononuclear cells (PBMC) \n Source: 10x Genomics sample dataset"
# data_list2 <- list(tcell_cd8=c("data_tcellcd8_results", "descrip_tcell_cd8"), 
#                   pbmc=c("data_pbmc_results", "descrip_pbmc"))


#####tailored demo
#load datasets
data_tcellcd8_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_full.rds")
data_pbmc_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_full.rds")

#dataset list for selection
data_full_list <- list(tcell_cd8="data_tcellcd8_full", pbmc="data_pbmc_full")



# user interface ----------------------------------------------------------

ui <- navbarPage(
  "Vitessce",
  
  ##### ui: basic demo ----------------------------------------------------
  tabPanel(
    "Pre-programmed demo",
    fluidPage(
      #select data
      h4("Dataset"),
      selectInput("dataset", label=NULL, choices=data_list),
      #selectInput("dataset_name", label=NULL, choices=data_names_list),
      
      #print dataset dimensions
      h4("Dataset dimensions"),
      htmlOutput("dataset_dimensions"),
      
      #create vitessce visualization
      h4("Vitessce visualization"),
      vitessce_output(output_id="vitessce_visualization", height="600px")
    )
  ),

  ##### ui: tailored demo -------------------------------------------------
  tabPanel("Tailored demo",
           fluidPage(
             #select data
             h4("Dataset"),
             selectInput("dataset_full", label=NULL, choices=data_full_list),
             
             #select filtering criteria
             h4("Data subsetting"),
             fluidRow(
               column(3, numericInput("user_min_cells", HTML("min.cells<br>(keep genes detected in at least <i>min.cells</i> cells)"), 100, min=0, max=NA)), #default value=100
               column(3, numericInput("user_min_features", HTML("min.features<br>(keep cells with at least <i>min.features</i> genes detected)"), 500, min=0, max=NA)), #default value=500
               column(3, numericInput("user_mt_gene_threshold", HTML("percent.mt<br>(keep cells with fewer than <i>percent.mt</i>% of genes mapping to mitochondrial genes)"), 5, min=0, max=100))
             ),
             
             #print dataset dimensions
             h4("Dataset dimensions"),
             #verbatimTextOutput("dataset_dimensions_tailored"),
             htmlOutput("dataset_dimensions_tailored"),
             
             #select panels to display in vitessce widget
             h4("Panels to include in Vitessce visualization"),
             fluidRow(
               column(3,
                      checkboxGroupInput("checkboxes_analyses", 
                                         label="Analyses",
                                         choices=list("PCA"="pca", "UMAP"="umap", "t-SNE"="tsne"),
                                         selected=c("pca", "umap", "tsne"))),
               column(3,
                      checkboxGroupInput("checkboxes_summaries", 
                                         label="Summaries",
                                         choices=list("Heatmap"="heatmap", "Cell set sizes"="cell_set_sizes"),
                                         selected=c("heatmap", "cell_set_sizes"))),
               column(3,
                      checkboxGroupInput("checkboxes_descrip", 
                                         label="Descriptions",
                                         choices=list("Dataset"="dataset_descrip", "Cell sets"="cell_sets", "Genes"="genes"),
                                         selected=c("dataset_descrip", "cell_sets", "genes"))),
               column(3,
                      checkboxGroupInput("checkboxes_view", 
                                         label="View options",
                                         choices=list("Link scatterplots"="link_scatterplots", "Light theme"="light_theme"),
                                         selected=c("link_scatterplots", "light_theme")))
               ), #end fluidRow
             
             # #test if data processing worked
             # h4("Test if data processing worked"),
             # htmlOutput("test_tailored"),
             
             #create vitessce visualization
             h4("Vitessce visualization"),
             vitessce_output(output_id="vitessce_visualization_tailored", height="600px")
             
             ) #end fluidPage
           ) #end tabPanel
)




# server ------------------------------------------------------------------

server <- function(input, output, session){
  
  ##### server: basic demo ------------------------------------------------
  data <- reactive({get(input$dataset)})
  # data_name <- reactive({get(input$dataset_name)})
  # data <- reactive({data_list[[data_name()]]})
  
  #print dimensions of dataset
  output$dataset_dimensions <- renderUI({
    dim(data())
    str_criteria <- "Filtering criteria <ul><li>min.cells = 100: keep genes detected in at least 100 cells</li><li>min.features = 500: keep cells with at least 500 genes detected</li></ul>"
    str_dim_data <- paste("Dataset dimensions:", dim(data())[1], "genes x", dim(data())[2], "cells")
    HTML(paste(str_criteria, str_dim_data, sep=""))
    })
  
  #vitessce visualization
  output$vitessce_visualization <- render_vitessce(expr={
    #create progress object
    progress <- shiny::Progress$new()
    progress$set(message = "", value = 0)
    on.exit(progress$close()) #close the progress bar when this reactive exits
    #function to update progress
    n <- 2
    updateProgress <- function(detail = NULL){
      progress$inc(amount = 1/n, detail = detail)
    }
    
    #set up widget
    updateProgress("Creating Vitessce visualization")
    vc <- VitessceConfig$new("My config")
    dataset <- vc$add_dataset("My dataset")
    dataset <- dataset$add_object(SeuratWrapper$new(data(), 
                                                    cell_set_meta_names=list("seurat_clusters"), 
                                                    num_genes=100))
    panel_scatterplot_pca <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="pca")
    panel_scatterplot_umap <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="umap")
    panel_scatterplot_tsne <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="tsne")
    panel_heatmap <- vc$add_view(dataset, Component$HEATMAP)
    #panel_status <- vc$add_view(dataset, Component$STATUS)
    panel_cellsets <- vc$add_view(dataset, Component$CELL_SETS)
    panel_cellset_sizes <- vc$add_view(dataset, Component$CELL_SET_SIZES)
    panel_genes <- vc$add_view(dataset, Component$GENES)
    panel_description <- vc$add_view(dataset, Component$DESCRIPTION)
    panel_description <- panel_description$set_props(description = "Test")
    vc$layout(hconcat(vconcat(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
                      vconcat(panel_heatmap, panel_cellset_sizes),
                      vconcat(panel_description, 
                              #panel_status, 
                              panel_cellsets, panel_genes)))

    vc$link_views(
      c(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
      c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
      c_values = c(1, 0, 0)
    )

    updateProgress("Complete!")
    vc$widget(theme="light")
    
  }) #end vitessce visualization output
  
  ##### server: tailored demo ---------------------------------------------
  #full dataset
  data_full <- reactive({get(input$dataset_full)})
  #subset dataset
  expr_matrix_subset <- reactive({GetAssayData(object=data_full(), slot="data")})
  data_subset <- reactive({
    data_subset_genes_and_cells <- CreateSeuratObject(counts=expr_matrix_subset(), project="subset", min.cells=input$user_min_cells, min.features=input$user_min_features) #subset data based on min.cells and min.features (user_min_cells and user_min_features)
    data_subset_genes_and_cells[["percent.mt"]] <- PercentageFeatureSet(data_subset_genes_and_cells, pattern="^MT-") #add mitochondrial genes column
    data_subset_mt_genes <- subset(data_subset_genes_and_cells, subset=percent.mt<input$user_mt_gene_threshold) #subset data based on mitochondrial genes (user_mt_gene_threshold)
    })
  
  #print dimensions of dataset
  output$dataset_dimensions_tailored <- renderUI({
    #print dimensions
    str_dim_data_full <- paste("Full dataset:", dim(data_full())[1], "genes x ", dim(data_full())[2], "cells")
    str_dim_data_subset <- paste("Subsetted dataset:", dim(data_subset())[1], "genes x ", dim(data_subset())[2], "cells")
    HTML(paste(str_dim_data_full, str_dim_data_subset, sep="<br/>"))
  })
  
  #analyze data, reactive
  data_tailored <- reactive({analyze_data(data_subset())})
  
  # #test if data analysis worked
  # output$test_tailored <- renderUI({
  #   #print dimensions
  #   str_test <- paste("Subsetted dataset:", dim(data_tailored())[1], "genes x ", dim(data_tailored())[2], "cells")
  #   HTML(paste(str_test, str_test, sep="<br/>"))
  # })
  
  #vitessce visualization
  output$vitessce_visualization_tailored <- render_vitessce(expr={
    #create progress object
    progress <- shiny::Progress$new()
    progress$set(message="", value = 0)
    on.exit(progress$close()) #close the progress bar when this reactive exits
    #function to update progress
    n <- 2
    updateProgress <- function(detail = NULL){
      progress$inc(amount = 1/n, detail = detail)
    }
    
    #set up widget
    updateProgress("Creating Vitessce visualization")
    vc <- VitessceConfig$new("My config")
    dataset <- vc$add_dataset("My dataset")
    
    #panels: analyses
    dataset <- dataset$add_object(SeuratWrapper$new(data_tailored(), 
                                                    cell_set_meta_names=list("seurat_clusters"), 
                                                    num_genes=100))
    #panels: analyses (column 1)
    reactive_column_analyses <- reactive({
      column_panels <- c()
      if("pca" %in% input$checkboxes_analyses){
        panel_scatterplot_pca <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="pca")
        column_panels <- append(column_panels, panel_scatterplot_pca)
      }
      if("umap" %in% input$checkboxes_analyses){
        panel_scatterplot_umap <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="umap")
        column_panels <- append(column_panels, panel_scatterplot_umap)
      }
      if("tsne" %in% input$checkboxes_analyses){
        panel_scatterplot_tsne <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="tsne")
        column_panels <- append(column_panels, panel_scatterplot_tsne)
      }
      column_panels
      })
    
    #panels: summaries (column 2)
    reactive_column_summaries <- reactive({
      column_panels <- c()
      if("heatmap" %in% input$checkboxes_summaries){
        panel_heatmap <- vc$add_view(dataset, Component$HEATMAP)
        column_panels <- append(column_panels, panel_heatmap)
      }
      if("cell_set_sizes" %in% input$checkboxes_summaries){
        panel_cellset_sizes <- vc$add_view(dataset, Component$CELL_SET_SIZES)
        column_panels <- append(column_panels, panel_cellset_sizes)
      }
      column_panels
    })
    
    #panels: description (column 3)
    reactive_column_descrip <- reactive({
      column_panels <- c()
      if("dataset_descrip" %in% input$checkboxes_descrip){
        panel_description <- vc$add_view(dataset, Component$DESCRIPTION)
        panel_description <- panel_description$set_props(description = "Test")
        column_panels <- append(column_panels, panel_description)
      }
      if("cell_sets" %in% input$checkboxes_descrip){
        panel_cellsets <- vc$add_view(dataset, Component$CELL_SETS)
        column_panels <- append(column_panels, panel_cellsets)
      }
      if("genes" %in% input$checkboxes_descrip){
        panel_genes <- vc$add_view(dataset, Component$GENES)
        column_panels <- append(column_panels, panel_genes)
      }
      column_panels
    })
    
    
    # #view options: link scatterplots
    # reactive_link_scatterplots <- reactive({
    #   if("XXX" %in% input$checkboxes_view){
    #     vc$layout(hconcat(do.call(vconcat, as.list(reactive_column_analyses())),
    #                       vconcat(reactive_column_summaries()[[1]], reactive_column_summaries()[[2]]),
    #                       vconcat(reactive_column_descrip()[[1]], reactive_column_descrip()[[2]], reactive_column_descrip()[[3]])
    #                       )
    #               )
    #   }
    #   if(!("XXX" %in% input$checkboxes_view)){
    #     vc$layout(hconcat(vconcat(reactive_column_analyses()[[1]], reactive_column_analyses()[[2]], reactive_column_analyses()[[3]]),
    #                       vconcat(reactive_column_summaries()[[1]], reactive_column_summaries()[[2]]),
    #                       vconcat(reactive_column_descrip()[[1]], reactive_column_descrip()[[2]], reactive_column_descrip()[[3]])
    #                       )
    #               )
    #     vc$link_views(c(reactive_column_analyses()),
    #                   c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
    #                   c_values = c(1, 0, 0)
    #     )}
    # 
    # })
    
    
    reactive_link_scatterplots <- reactive({
      if("link_scatterplots" %in% input$checkboxes_view){
        vc$link_views(reactive_column_analyses(),
                      c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
                      c_values=c(1, 0, 0)
        )}
      #if(!("link_scatterplots" %in% input$checkboxes_view)){}
      else{
        for(view in reactive_column_analyses()){
          vc$link_views(c(view),
                        c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
                        c(1, 0, 0)
                        )
          }
        }
    }) #end reactive: reactive_link_scatterplots
    

    #view options: theme
    reactive_light_theme <- reactive({
      if("light_theme" %in% input$checkboxes_view){vc$widget(theme="light")}
      else{vc$widget(theme="dark")}
    })
    
    #layout panels: run reactives to create/update columns
    vc$layout(hconcat(do.call(vconcat, as.list(reactive_column_analyses())),
                      do.call(vconcat, as.list(reactive_column_summaries())),
                      do.call(vconcat, as.list(reactive_column_descrip()))
                      )
              )
    
    #working version
    # vc$layout(hconcat(vconcat(reactive_column_analyses()[[1]], reactive_column_analyses()[[2]], reactive_column_analyses()[[3]]),
    #                   vconcat(reactive_column_summaries()[[1]], reactive_column_summaries()[[2]]),
    #                   vconcat(reactive_column_descrip()[[1]], reactive_column_descrip()[[2]], reactive_column_descrip()[[3]])
    #                   )
    # )
    
    # #reactive_link_scatterplots()
    # vc$link_views(
    #   reactive_column_analyses(),
    #   c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
    #   c_values = c(1, 0, 0)
    # )
    
    reactive_link_scatterplots()
    
    updateProgress("Complete!")
    reactive_light_theme()
    
  })
  
  
}



# app ---------------------------------------------------------------------

shinyApp(ui, server)
