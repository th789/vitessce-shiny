
library(shiny)
library(vitessce)
library(Seurat)

#####basic demo
#load datasets
data_tcellcd8_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_results.rds")
data_pbmc_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_results.rds")
#dataset list for selection
data_list <- list(tcell_cd8="data_tcellcd8_results", pbmc="data_pbmc_results")


#####tailored demo
#load datasets
data_tcellcd8_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_full.rds")
data_pbmc_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_full.rds")
#dataset list for selection
data_full_list <- list(tcell_cd8="data_tcellcd8_full", pbmc="data_pbmc_full")



# user interface ----------------------------------------------------------

ui <- navbarPage(
  "Vitessce",
  
  tabPanel(
    "Pre-programmed demo",
    fluidPage(
      #select data
      h4("Dataset"),
      selectInput("dataset", label=NULL, choices=data_list),
      
      #print dataset dimensions
      h4("Dataset dimensions"),
      p("Filtering criteria: genes detected in at least 100 cells x cells with at least 500 genes detected"),
      verbatimTextOutput("dataset_dimensions"),
      
      #create vitessce visualization
      h4("Vitessce visualization"),
      vitessce_output(output_id="vitessce_visualization", height="600px")
    )
  ),

  tabPanel("Tailored demo",
           fluidPage(
             #select data
             h4("Dataset"),
             selectInput("dataset_full", label=NULL, choices=data_full_list),
             
             #select filtering criteria
             numericInput("user_min_cells", "min.cells", 100, min=0, max=NA), #default value=100
             numericInput("user_min_features", "min.features", 500, min=0, max=NA), #default value=500
             
             #print dataset dimensions
             h4("Dataset dimensions"),
             #verbatimTextOutput("dataset_dimensions_tailored"),
             htmlOutput("dataset_dimensions_tailored")
             )
           )
)




# server ------------------------------------------------------------------

server <- function(input, output, session){
  #####basic demo
  data <- reactive({get(input$dataset)})
  
  #print dimensions of dataset
  output$dataset_dimensions <- renderPrint({
    dim(data())
    paste(dim(data())[1], "genes x", dim(data())[2], "cells")
    })
  
  #vitessce visualization
  output$vitessce_visualization <- render_vitessce(expr = {
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
    
  })
  
  #####tailored demo
  #full dataset
  data_full <- reactive({get(input$dataset_full)})
  #subset dataset
  expr_matrix_subset <- reactive({GetAssayData(object=data_full(), slot="data")})
  data_subset <- reactive({CreateSeuratObject(counts=expr_matrix_subset(), project="subset", min.cells=input$user_min_cells, min.features=input$user_min_features)})

  #print dimensions of dataset
  output$dataset_dimensions_tailored <- renderUI({
    #print dimensions
    str_dim_data_full <- paste("Full dataset dimensions:", dim(data_full())[1], "genes x ", dim(data_full())[2], "cells")
    str_dim_data_subset <- paste("Subsetted dataset dimensions:", dim(data_subset())[1], "genes x ", dim(data_subset())[2], "cells")
    HTML(paste(str_dim_data_full, str_dim_data_subset, sep="<br/>"))
  })
  
  
}

shinyApp(ui, server)
