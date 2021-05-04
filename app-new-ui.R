# load data & packages ----------------------------------------------------

library(shiny)
library(vitessce)
library(Seurat)
source("tailored-demo-helpers.R")


#####basic demo
#load datasets
data_tcellcd4_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd4_results.rds")
data_tcellcd8_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_results.rds")
data_pbmc_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_results.rds")
data_lung_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_lung_results.rds")
data_nsclc_results <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_nsclc_results.rds")

#create pairwise lists
list_choices_names <- c("tcell_cd4"="tcell_cd4", "tcell_cd8"="tcell_cd8", "pbmc"="pbmc", "lung"="lung", "nsclc"="nsclc") #name-name
list_choices_names_dfs <- c("tcell_cd4"=data_tcellcd4_results, "tcell_cd8"=data_tcellcd8_results, "pbmc"=data_pbmc_results, "lung"=data_lung_results, "nsclc"=data_nsclc_results) #name-df
list_choices_names_descrip <- c("tcell_cd4"="CD4 T cells -- Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017).", 
                                "tcell_cd8"="CD8 T cells -- Zheng, G., Terry, J., Belgrader, P. et al. Massively parallel digital transcriptional profiling of single cells. Nat Commun 8, 14049 (2017).", 
                                "pbmc"="Peripheral blood mononuclear cells (PBMC) -- 10X Genomics \nhttps://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k", 
                                "lung"="Lung cells -- Travaglini, K.J., Nabhan, A.N., Penland, L. et al. A molecular cell atlas of the human lung from single-cell RNA sequencing. Nature 587, 619–625 (2020). ", 
                                "nsclc"="Non-small cell lung cancer -- 10X Genomics https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_nsclc_5gex") #name-description


#####tailored demo
#load datasets
data_tcellcd4_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd4_full.rds")
data_tcellcd8_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_full.rds")
data_pbmc_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_full.rds")
data_lung_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_lung_full.rds")
data_nsclc_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_nsclc_full.rds")

#dataset list for selection
data_full_list <- list(tcell_cd4="data_tcellcd4_full", tcell_cd8="data_tcellcd8_full", pbmc="data_pbmc_full", lung="data_lung_full", nsclc="data_nsclc_full")

list_choices_names_dfs_tailored <- c("tcell_cd4"=data_tcellcd4_full, "tcell_cd8"=data_tcellcd8_full, "pbmc"=data_pbmc_full, "lung"=data_lung_full, "nsclc"=data_nsclc_full) #name-df



# shiny app settings ------------------------------------------------------

options(shiny.maxRequestSize = 500*1024^2) #limit file size to 500MB (for file upload)


# dynamic ui tabs ---------------------------------------------------------


tabs_input_data <- tabsetPanel(
  id="tailored_demo_input_data",
  type="hidden",
  tabPanel("select_data",
           selectInput("dataset_full", label="Select example dataset", choices=list_choices_names)
  ),
  tabPanel("upload_data", 
           fileInput("user_dataset", "Upload dataset (SeuratObject in .rds file)", accept=".rds")
  )
)




# ui panels ---------------------------------------------------------------

## basic demo -------------------------------------------------------------

### sidebarpanel ----------------------------------------------------------

basic_demo_sidebarpanel <- sidebarPanel(
  
  #specify sidebarPanel features height, width, and scroll bar
  width=3,
  style = "position: fixed; height: 88.3vh; width: 42vh; overflow-y: auto;",
  
  #1. select dataset from list of examples
  h4("Dataset"),
  selectInput("dataset", label=NULL, choices=list_choices_names),
  
  #2. check dataset dimensions
  h4("Dataset dimensions"),
  htmlOutput("dataset_dimensions")
)


### mainpanel -------------------------------------------------------------

basic_demo_mainpanel <- mainPanel(
  #create vitessce visualization
  #h4("Vitessce visualization"),
  vitessce_output(output_id="vitessce_visualization", height="700px", width="1030px")
)

## tailored demo ----------------------------------------------------------

### sidebarpanel ----------------------------------------------------------

#sidebarpanel
tailored_demo_sidebarpanel <- sidebarPanel(
  
  #specify sidebarPanel features height, width, and scroll bar
  width=3,
  style = "position: fixed; height: 88.3vh; width: 42vh; overflow-y: auto;",
  
  ###1. specify dataset
  h4("1. Specify dataset"),
  #input data type: select dataset or upload dataset
  selectInput(inputId="tailored_demo_input", label="Input data", 
              choices = c("Select example dataset"="select_data", "Upload dataset"="upload_data")
              ),
  #based on input data type: drop-down list (select dataset) or data browser (upload dataset)
  tabs_input_data,
  
  ###2. perform quality control
  h4("2. Perform quality control (filter dataset)"),
  numericInput("user_min_cells", HTML("min.cells<br><span style='font-weight:normal'>keep genes detected in at least <i>min.cells</i> cells</span>"), 100, min=0, max=NA), #default value=100
  numericInput("user_min_features", HTML("min.features<br><span style='font-weight:normal'>keep cells with at least <i>min.features</i> genes detected</span>"), 500, min=0, max=NA), #default value=500
  numericInput("user_mt_gene_threshold", HTML("percent.mt<br><span style='font-weight:normal'>keep cells with less than <i>percent.mt</i>% of genes mapping to mitochondrial genes</span>"), 5, min=0, max=100), #default value=5
  
  ###3. check dataset dimensions
  h4("3. Check dataset dimensions"),
  htmlOutput("dataset_dimensions_tailored"),
  
  ###4. specify vitessce visualization parameters
  h4("4. Specify Vitessce visualization parameters"),
  #row1: analyses and summaries
  fluidRow(
    column(6, checkboxGroupInput("checkboxes_analyses", label="Analyses",
                                 choices=list("PCA"="pca", "UMAP"="umap", "t-SNE"="tsne"),
                                 selected=c("pca", "umap", "tsne"))
           ),
    column(6, checkboxGroupInput("checkboxes_summaries", label="Summaries",
                       choices=list("Heatmap"="heatmap", "Cell set sizes"="cell_set_sizes"),
                       selected=c("heatmap", "cell_set_sizes"))
           ),
    ), #end fluidRow (row1)
  #row2: description and view options
  fluidRow(
    column(6, checkboxGroupInput("checkboxes_descrip", label="Descriptions",
                       choices=list("Dataset"="dataset_descrip", "Cell sets"="cell_sets", "Genes"="genes"),
                       selected=c("dataset_descrip", "cell_sets", "genes"))
           ),
    column(6, checkboxGroupInput("checkboxes_view", label="View options",
                       choices=list("Link scatterplots"="link_scatterplots", "Light theme"="light_theme"),
                       selected=c("link_scatterplots", "light_theme"))
           ),
    ) #end fluidRow (row2)
  
  ) #end sidebarPanel (end tailored_demo_sidebarpanel)


### mainpanel -------------------------------------------------------------

#main panel
tailored_demo_mainpanel <- mainPanel(
  #h4("Vitessce visualization"),
  vitessce_output(output_id="vitessce_visualization_tailored", height="700px", width="1030px")
)

# ui ----------------------------------------------------------------------


# ui <- fluidPage(
#   sidebarLayout(
#     tailored_demo_sidebarpanel,
#     tailored_demo_mainpanel
#   )
# )



ui <- navbarPage(
  "Vitessce Shiny",
  
  ##### ui: basic demo ----------------------------------------------------
  tabPanel(
    "Demo",
    fluidPage(
      sidebarLayout(basic_demo_sidebarpanel, basic_demo_mainpanel)
    ) #end fluidPage
  ), #end tabPanel
  
  ##### ui: tailored demo -------------------------------------------------
  tabPanel("Run analysis",
           fluidPage(
             sidebarLayout(tailored_demo_sidebarpanel, tailored_demo_mainpanel)
             ) #end fluidPage
           ) #end tabPanel
) #end navbarPage (end ui)



# server ------------------------------------------------------------------



server <- function(input, output, session){
## basic demo -------------------------------------------------------------
  
  ###1. get data (based on selected dataset in input)
  data <- reactive({list_choices_names_dfs[[input$dataset]]}) #get dataset
  data_descrip <- reactive({list_choices_names_descrip[[input$dataset]]}) #get dataset description

  ###2. print dimensions of dataset
  output$dataset_dimensions <- renderUI({
    dim(data())
    str_criteria <- "Quality control (filtering criteria) <ul><li>min.cells = 100: keep genes detected in at least 100 cells</li><li>min.features = 500: keep cells with at least 500 genes detected</li><li>percent.mt = 5: keep cells with less than 5% of genes mapping to mitochondrial genes</li></ul>"
    str_dim_data <- paste("Dataset dimensions:", dim(data())[1], "genes x", dim(data())[2], "cells")
    HTML(paste(str_criteria, str_dim_data, sep=""))
  })
  
  ###3. vitessce visualization
  output$vitessce_visualization <- render_vitessce(expr={
    #create progress object
    progress <- shiny::Progress$new()
    progress$set(message = "", value = 0)
    on.exit(progress$close()) #close the progress bar when this reactive exits
    #function to update progress
    n <- 3
    updateProgress <- function(detail = NULL){
      progress$inc(amount = 1/n, detail = detail)
    }
    
    #vitessce --- set up widget
    updateProgress("Demo: creating individual visualizations")
    vc <- VitessceConfig$new("My config")
    dataset <- vc$add_dataset("My dataset")
    dataset <- dataset$add_object(SeuratWrapper$new(data(), 
                                                    cell_set_meta_names=list("seurat_clusters"), 
                                                    num_genes=100))
    
    
    #vitessce --- add views (panels)
    panel_scatterplot_pca <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="pca")
    panel_scatterplot_umap <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="umap")
    panel_scatterplot_tsne <- vc$add_view(dataset, Component$SCATTERPLOT, mapping="tsne")
    panel_heatmap <- vc$add_view(dataset, Component$HEATMAP)
    #panel_status <- vc$add_view(dataset, Component$STATUS)
    panel_cellsets <- vc$add_view(dataset, Component$CELL_SETS)
    panel_cellset_sizes <- vc$add_view(dataset, Component$CELL_SET_SIZES)
    panel_genes <- vc$add_view(dataset, Component$GENES)
    panel_description <- vc$add_view(dataset, Component$DESCRIPTION)
    panel_description <- panel_description$set_props(description=data_descrip())
    updateProgress("Demo: assembling Vitessce visualization")
    vc$layout(hconcat(vconcat(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
                      vconcat(panel_heatmap, panel_cellset_sizes),
                      vconcat(panel_description, 
                              #panel_status, 
                              panel_cellsets, panel_genes)))
    
    #vitessce --- link scatterplots
    vc$link_views(
      c(panel_scatterplot_pca, panel_scatterplot_umap, panel_scatterplot_tsne),
      c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
      c_values = c(1, 0, 0)
    )
    
    #update progress bar
    updateProgress("Demo: complete!") 
    
    #vitessce --- specify theme
    vc$widget(theme="light")
    
  }) #end vitessce visualization output
  
## tailored demo ----------------------------------------------------------

  ###1. obtain full dataset 
  #based on input$tailored_demo_input: selected dataset or uploaded dataset
  observeEvent(input$tailored_demo_input, {
    updateTabsetPanel(inputId="tailored_demo_input_data", selected=input$tailored_demo_input)
    }) 
  #create data_full() reactive by getting selected dataset or loading uploaded dataset
  data_full <- reactive({
    switch(input$tailored_demo_input,
           select_data = list_choices_names_dfs_tailored[[input$dataset_full]],
           upload_data=readRDS(input$user_dataset$datapath)
    )
  })
  
  #create data_descrip_tailored() reactive to get dataset descriptionn
  data_descrip_tailored <- reactive({
    switch(input$tailored_demo_input,
           select_data=list_choices_names_descrip[[input$dataset_full]],
           upload_data="My data"
    )
  })
  
  ###2. perform quality control: filter dataset
  expr_matrix_subset <- reactive({GetAssayData(object=data_full(), slot="data")})
  data_subset <- reactive({
    data_subset_genes_and_cells <- CreateSeuratObject(counts=expr_matrix_subset(), project="subset", min.cells=input$user_min_cells, min.features=input$user_min_features) #subset data based on min.cells and min.features (user_min_cells and user_min_features)
    data_subset_genes_and_cells[["percent.mt"]] <- PercentageFeatureSet(data_subset_genes_and_cells, pattern="^MT-") #add mitochondrial genes column
    data_subset_mt_genes <- subset(data_subset_genes_and_cells, subset=percent.mt<input$user_mt_gene_threshold) #subset data based on mitochondrial genes (user_mt_gene_threshold)
  })
  
  ###3. print dataset dimensions
  output$dataset_dimensions_tailored <- renderUI({
    str_dim_data_full <- paste("<b>Full dataset</b><br>", dim(data_full())[1], "genes x ", dim(data_full())[2], "cells<br> ")
    str_dim_data_subset <- paste("<b>Subsetted dataset</b><br>", dim(data_subset())[1], "genes x ", dim(data_subset())[2], "cells<br> <br>")
    #print dimensions
    HTML(paste(str_dim_data_full, str_dim_data_subset, sep="<br/>"))
  })

  ###4. create Vitessce visualization
  #analyze data
  data_tailored <- reactive({analyze_data(data_subset())})
  
  #vitessce visualization
  output$vitessce_visualization_tailored <- render_vitessce(expr={
    #create progress object
    progress <- shiny::Progress$new()
    progress$set(message="", value=0)
    on.exit(progress$close()) #close the progress bar when this reactive exits
    #function to update progress
    n <- 3
    updateProgress <- function(detail = NULL){
      progress$inc(amount = 1/n, detail = detail)
    }
    
    #vitessce --- set up widget
    updateProgress("Analysis: creating individual visualizations")
    vc <- VitessceConfig$new("My config")
    dataset <- vc$add_dataset("My dataset")
    dataset <- dataset$add_object(SeuratWrapper$new(data_tailored(), 
                                                    cell_set_meta_names=list("seurat_clusters"), 
                                                    num_genes=100))
    ###create reactives based on inputs
    #reactive: panels, analyses (column 1)
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
    
    #reactive: panels, summaries (column 2)
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
    
    #reactive: panels, description (column 3)
    reactive_column_descrip <- reactive({
      column_panels <- c()
      if("dataset_descrip" %in% input$checkboxes_descrip){
        panel_description <- vc$add_view(dataset, Component$DESCRIPTION)
        panel_description <- panel_description$set_props(description=data_descrip_tailored())
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
    
    updateProgress("Analysis: assembling Vitessce visualization")
    #reactive: view options, link or unlink scatterplots
    reactive_link_scatterplots <- reactive({
      #link scatterplots
      if("link_scatterplots" %in% input$checkboxes_view){
        vc$link_views(reactive_column_analyses(),
                      c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
                      c_values=c(1, 0, 0)
                      )
        } #end "if" statement
      #unlink scatterplots
      else{
        for(view in reactive_column_analyses()){
          vc$link_views(c(view),
                        c(CoordinationType$EMBEDDING_ZOOM, CoordinationType$EMBEDDING_TARGET_X, CoordinationType$EMBEDDING_TARGET_Y),
                        c(1, 0, 0)
                        )
          } #end for loop: for(view in reactive_column_analyses())
        } #end "else" statement
    }) #end reactive: reactive_link_scatterplots
    
    #reactive: view options, light or dark theme
    reactive_light_theme <- reactive({
      #light theme
      if("light_theme" %in% input$checkboxes_view){vc$widget(theme="light")}
      #dark theme
      else{vc$widget(theme="dark")}
    })
    
    
    #panel_columns: list of columns (in Vitessce visualization) that are not empty
    panel_columns <- list(analyses=reactive_column_analyses(), summaries=reactive_column_summaries(), descrip=reactive_column_descrip())
    #return error message if no panels (visualizations) are selected
    empty_visualization <- reactive({length(reactive_column_analyses())==0 & length(reactive_column_summaries())==0 & length(reactive_column_descrip())==0})
    validate(need(!empty_visualization(), "Please select at least one visualization"))
    #remove list of columns (reactive_column_analyses, reactive_column_summaries, and reactive_column_descrip) that are empty
    if(length(reactive_column_analyses())==0){panel_columns <- within(panel_columns, rm(analyses))}
    if(length(reactive_column_summaries())==0){panel_columns <- within(panel_columns, rm(summaries))}
    if(length(reactive_column_descrip())==0){panel_columns <- within(panel_columns, rm(descrip))}
    
    #apply do.call(vconcat) to every list (make as list) in panel_columns
    panel_columns_vconcat <- lapply(panel_columns, function(x){do.call(vconcat, as.list(x))})
    #apply do.call(hconcat) to the list 'panel_columns_vconcat'
    panel_columns_hconcat <- do.call(hconcat, panel_columns_vconcat)
    
    #vitessce --- add views: create/update columns and panels
    vc$layout(panel_columns_hconcat)
    
    # #vitessce --- add views: use reactives to create/update/layout columns and panels
    # vc$layout(hconcat(do.call(vconcat, as.list(reactive_column_analyses())),
    #                   do.call(vconcat, as.list(reactive_column_summaries())),
    #                   do.call(vconcat, as.list(reactive_column_descrip()))
    # )
    # )
    
    #vitessce --- link or unlink scatterplots
    reactive_link_scatterplots() 
    
    #update progress bar
    updateProgress("Analysis: complete!")
    
    #vitessce --- specify theme (light or dark)
    reactive_light_theme()
    
  }) #end vitessce visualization 
  

} #end server


# compile app -------------------------------------------------------------

shinyApp(ui=ui,server=server)


