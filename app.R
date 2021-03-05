
library(shiny)
library(vitessce)


#load datasets
data_tcellcd8 <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/vitessce-shiny/data/data_tcellcd8.rds")
data_pbmc <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/vitessce-shiny/data/data_pbmc.rds")
#dataset list for selection
data_list <- list(tcell_cd8="data_tcellcd8", pbmc="data_pbmc")


#user interface
ui <- fluidPage(
  
  titlePanel("Vitessce demo"),
  
  "Select data:",
  selectInput("dataset", label="Dataset", choices=data_list),
  
  "Dataset dimensions",
  verbatimTextOutput("dataset_dimensions"),
  
  "Vitessce visualization",
  vitessce_output(output_id="vitessce_visualization", height="600px")
  
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
    #on.exit(progress$close()) #close the progress when this reactive exits (even if there's an error)
    #function to update progress
    n <- 6
    updateProgress <- function(detail = NULL) {
      progress$inc(amount = 1/n, detail = detail)
    }
    
    #data processing
    updateProgress("Step 1 of 3: Processing data")
    data <- data()
    data[["percent.mt"]] <- Seurat::PercentageFeatureSet(data, pattern = "^MT-")
    data <- subset(data, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)
    data <- Seurat::NormalizeData(data, normalization.method="LogNormalize", scale.factor=10000)
    data <- Seurat::FindVariableFeatures(data, selection.method="vst", nfeatures=2000)
    all.genes <- rownames(data)
    
    #data analysis
    updateProgress("Step 2 of 3: Analyzing data (PCA)")
    data <- Seurat::ScaleData(data, features=all.genes) #PCA
    data <- Seurat::RunPCA(data, features=Seurat::VariableFeatures(object=data)) 
    data <- Seurat::FindNeighbors(data, dims=1:10)
    data <- Seurat::FindClusters(data, resolution=0.5)
    updateProgress("Step 2 of 3: Analyzing data (UMAP)")
    data <- Seurat::RunUMAP(data, dims=1:20, min.dist=0.75) #UMAP
    updateProgress("Step 2 of 3: Analyzing data (t-SNE)")
    data <- Seurat::RunTSNE(data, dims=1:10, nthreads=4, max_iter=1000, k.seed=12345) #t-SNE
    
    #set up widget
    updateProgress("Step 3 of 3: Creating Vitessce visualization")
    vc <- VitessceConfig$new("My config")
    dataset <- vc$add_dataset("My dataset")
    dataset <- dataset$add_object(SeuratWrapper$new(data, cell_set_meta_names=list("seurat_clusters")))

    panel_scatterplot_pca <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="pca")
    panel_scatterplot_umap <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="umap")
    panel_scatterplot_tsne <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="tsne")
    panel_status <- vc$add_view(dataset, Component$STATUS)
    panel_cellsets <- vc$add_view(dataset, Component$CELL_SETS)
    vc$layout(hconcat(vconcat(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
                      vconcat(panel_cellsets, panel_status)))
    
    vc$link_views(
      c(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
      c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
      c_values = c(1, 0, 0)
    )
    
    updateProgress("Analysis complete!")
    vc$widget()
    #vc$widget(theme="light")
    
    
  })
  
}

shinyApp(ui, server)
