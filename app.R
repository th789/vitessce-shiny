#load libraries
library(shiny)
library(vitessce)



#load datasets
data_tcellcd8 <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/vitessce-shiny/data/data_tcellcd8.rds")
data_pbmc <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/vitessce-shiny/data/data_pbmc.rds")
#dataset list for selection
data_list <- list(tcell_cd8="data_tcellcd8", pbmc="data_pbmc")



#user interface
ui <- fluidPage(
  "Hello, world!",
  selectInput("dataset", label="Dataset", choices=data_list),

  "Dataset dimensions",
  verbatimTextOutput("dataset_dimensions"),
  #tableOutput("table"),

  "Vitessce visualization",
  vitessce_output("vitessce_visualization"),

)



#server
server <- function(input, output, session) {
  dataset <- reactive({get(input$dataset)})

  #print dimensions of dataset
  output$dataset_dimensions <- renderPrint({
    dim(dataset())
  })

  output$vitessce_visualization <- render_vitessce({
    #data processing
    data <- dataset()
    data[["percent.mt"]] <- Seurat::PercentageFeatureSet(data, pattern = "^MT-")
    data <- subset(data, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)
    data <- Seurat::NormalizeData(data, normalization.method="LogNormalize", scale.factor=10000)
    data <- Seurat::FindVariableFeatures(data, selection.method="vst", nfeatures=2000)
    all.genes <- rownames(data)
    data <- Seurat::ScaleData(data, features=all.genes)
    #PCA
    data <- Seurat::RunPCA(data, features=Seurat::VariableFeatures(object=data))
    data <- Seurat::FindNeighbors(data, dims=1:10)
    data <- Seurat::FindClusters(data, resolution=0.5)
    #UMAP
    data <- Seurat::RunUMAP(data, dims=1:20, min.dist=0.75)
    #t-SNE
    data <- Seurat::RunTSNE(data, dims=1:10, nthreads=4, max_iter=1000, k.seed=12345)

    #create Vitessce view config
    vc <- vitessce::VitessceConfig$new("My config")

    #add dataset
    dataset <- vc$add_dataset("My dataset")$add_object(vitessce::SeuratWrapper$new(data, cell_set_meta_names=list("seurat_clusters")))

    #create panels
    #PCA
    panel_scatterplot_pca <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="pca")
    #UMAP
    panel_scatterplot_umap <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="umap")
    #t-SNE
    panel_scatterplot_tsne <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="tsne")
    #status
    panel_status <- vc$add_view(dataset, Component$STATUS)
    #cell sets
    panel_cellsets <- vc$add_view(dataset, Component$CELL_SETS)

    #layout
    vc$layout(hconcat(vconcat(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
                      vconcat(panel_cellsets, panel_status)))

    #render the Vitessce widget
    vc$widget()

  })


}



shinyApp(ui, server)


