
library(shiny)
library(vitessce)


#load datasets
data_tcellcd8_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_results.rds")
data_pbmc_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_results.rds")
#dataset list for selection
data_list <- list(tcell_cd8="data_tcellcd8_results", pbmc="data_pbmc_results")


#user interface
ui <- navbarPage(
  "Vitessce",
  
  tabPanel(
    "Pre-programmed demo",
    fluidPage(
      #select data
      selectInput("dataset", label="Dataset:", choices=data_list),

      "Dataset dimensions",
      verbatimTextOutput("dataset_dimensions"),

      "Vitessce visualization",
      vitessce_output(output_id="vitessce_visualization", height="600px")
    )
  ),

  tabPanel("Tailored demo")
)



#server
server <- function(input, output, session){
  data <- reactive({get(input$dataset)})
  
  #print dimensions of dataset
  output$dataset_dimensions <- renderPrint({
    dim(data())
    })
  
  #vitessce visualization
  output$vitessce_visualization <- render_vitessce(expr = {
    #create progress object
    progress <- shiny::Progress$new()
    progress$set(message = "", value = 0)
    on.exit(progress$close()) #close the progress bar when this reactive exits
    #function to update progress
    n <- 3
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

    updateProgress("Setting up views")
    panel_scatterplot_pca <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="pca")
    panel_scatterplot_umap <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="umap")
    panel_scatterplot_tsne <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="tsne")
    panel_heatmap <- vc$add_view(dataset, Component$HEATMAP)
    panel_status <- vc$add_view(dataset, Component$STATUS)
    panel_cellsets <- vc$add_view(dataset, Component$CELL_SETS)
    panel_cellset_sizes <- vc$add_view(dataset, Component$CELL_SET_SIZES)
    panel_genes <- vc$add_view(dataset, Component$GENES)
    vc$layout(hconcat(vconcat(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
                      vconcat(panel_heatmap, panel_cellset_sizes),
                      vconcat(panel_status, panel_cellsets, panel_genes)))

    vc$link_views(
      c(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
      c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
      c_values = c(1, 0, 0)
    )

    updateProgress("Complete!")
    vc$widget()
    
  })
  
}

shinyApp(ui, server)


