
library(Seurat)

#setwd("~/Dropbox/ddesktop/lab-gehlenborg/vitessce-shiny")


# function to process 10x data --------------------------------------------

process_10x_data <- function(directory, data_name) {

  #read in 10x count data
  data_10x <- Read10X(data.dir = directory)
  #create SeuratObject
  data_seurat <- CreateSeuratObject(counts=data_10x, project=data_name, min.cells = 3, min.features = 200)
  #print data dimensions
  print(dim(data_seurat))
  #save data as RDS
  file_name <- paste0("data/", data_name, ".rds")
  saveRDS(data_seurat, file=file_name)
}


# pbmc --------------------------------------------------------------------

#read in 10x count data
counts_matrix_filename <- "~/Dropbox/ddesktop/lab-gehlenborg/seurat/filtered_gene_bc_matrices/hg19"
pbmc.data <- Read10X(data.dir = counts_matrix_filename)
#create SeuratObject
data_pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 100, min.features = 500)
#save data as RDS
saveRDS(data_pbmc, file = "data/data_pbmc.rds")
#dimensions
dim(data_pbmc) #4662 x 2482

# t-cells CD8 -------------------------------------------------------------

#read in 10x count data
counts_matrix_filename <- "~/Dropbox/ddesktop/lab-gehlenborg/datasets/tcellcd8/filtered_matrices_mex/hg19"
counts <- Read10X(data.dir = counts_matrix_filename)
#create SeuratObject
data_tcellcd8 <- CreateSeuratObject(counts = counts, min.cells = 100, min.features = 500, project = "tcellcd8")
#save data as RDS
saveRDS(data_tcellcd8, file = "data/data_tcellcd8.rds")
#dimensions
dim(data_tcellcd8) #6171 x 7856
