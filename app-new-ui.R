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



#####tailored demo
#load datasets
data_tcellcd8_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_full.rds")
data_pbmc_full <- readRDS("~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_full.rds")

#dataset list for selection
data_full_list <- list(tcell_cd8="data_tcellcd8_full", pbmc="data_pbmc_full")



# shiny app settings ------------------------------------------------------

options(shiny.maxRequestSize = 500*1024^2) #limit file size to 500MB (for file upload)


# dynamic ui tabs ---------------------------------------------------------


tabs_input_data <- tabsetPanel(
  id="tailored_demo_input_data",
  type="hidden",
  tabPanel("select_data",
           selectInput("dataset_full", label="Select example dataset", choices=data_full_list)
  ),
  tabPanel("upload_data", 
           fileInput("user_dataset", "Upload dataset (SeuratObject in .rds file)", accept=".rds")
  )
)




# ui panels ---------------------------------------------------------------


#sidebarpanel
tailored_demo_sidebarpanel <- sidebarPanel(
  
  #specify sidebarPanel features height, width, and scroll bar
  style = "position: fixed; height: 87vh; width: 40vh; overflow-y: auto;",
  
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
  numericInput("user_min_cells", HTML("min.cells<br>(keep genes detected in at least <i>min.cells</i> cells)"), 100, min=0, max=NA), #default value=100
  numericInput("user_min_features", HTML("min.features<br>(keep cells with at least <i>min.features</i> genes detected)"), 500, min=0, max=NA), #default value=500
  numericInput("user_mt_gene_threshold", HTML("percent.mt<br>(keep cells with less than <i>percent.mt</i>% of genes mapping to mitochondrial genes)"), 5, min=0, max=100), #default value=5
  
  
  
  
  
  
  
  )


#main panel
tailored_demo_mainpanel <- mainPanel(
  htmlOutput("dataset_dimensions_tailored"),
)

# ui ----------------------------------------------------------------------


# ui <- fluidPage(
#   sidebarLayout(
#     tailored_demo_sidebarpanel,
#     tailored_demo_mainpanel
#   )
# )



ui <- navbarPage(
  "Vitessce",
  
  ##### ui: basic demo ----------------------------------------------------
  tabPanel(
    "Basic demo",
    fluidPage(

    )
  ),
  
  ##### ui: tailored demo -------------------------------------------------
  tabPanel("Tailored demo",
           fluidPage(
             sidebarLayout(tailored_demo_sidebarpanel, tailored_demo_mainpanel)
             ) #end fluidPage
           ) #end tabPanel
) #end navbarPage (end ui)



# server ------------------------------------------------------------------



server <- function(input, output, session) {
  observeEvent(input$tailored_demo_input, {
    updateTabsetPanel(inputId="tailored_demo_input_data", selected=input$tailored_demo_input)
  }) 
  
  
  
  data_full <- reactive({
    switch(input$tailored_demo_input,
           select_data=get(input$dataset_full),
           upload_data=readRDS(input$user_dataset$datapath)
    )
  })
  
  
  #print dimensions of dataset
  output$dataset_dimensions_tailored <- renderUI({
    #print dimensions
    str_dim_data_full <- paste("Full dataset:", dim(data_full())[1], "genes x ", dim(data_full())[2], "cells")
    HTML(paste(str_dim_data_full, str_dim_data_full, sep="<br/>"))
  })
}


shinyApp(ui=ui,server=server)


