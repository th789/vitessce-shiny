
library(Seurat)
setwd("~/Dropbox/ddesktop/lab-gehlenborg/")


# function to process 10x data --------------------------------------------

read_10x_data <- function(directory, data_name, min_cells, min_features){
  #read in 10x count data
  data_10x <- Read10X(data.dir=directory)
  #create SeuratObject
  data_seurat <- CreateSeuratObject(counts=data_10x, project=data_name, min.cells=min_cells, min.features=min_features)
  #quality control: detect amount of mitochondrial DNA
  data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
  #print data dimensions
  print("dataset dimensions (genes x cells):")
  print(dim(data_seurat))
  
  return(data_seurat)
}



# function to analyze data ------------------------------------------------

#data in SeuratObject form
analyze_data <- function(data){
  #process data
  data <- NormalizeData(data, normalization.method="LogNormalize", scale.factor=10000)
  data <- FindVariableFeatures(data, selection.method="vst", nfeatures=2000)
  all.genes <- rownames(data)
  #data analysis
  #PCA
  data <- ScaleData(data, features=all.genes)
  data <- RunPCA(data, features=VariableFeatures(object=data)) 
  data <- FindNeighbors(data, dims=1:10)
  data <- FindClusters(data, resolution=0.5)
  #UMAP
  data <- RunUMAP(data, dims=1:20, min.dist=0.75)
  #t-SNE
  data <- RunTSNE(data, dims=1:10, nthreads=4, max_iter=1000, k.seed=12345)
  
  return(data)
}


# pbmc --------------------------------------------------------------------

#full dataset
data_pbmc_full <- read_10x_data(directory="~/Dropbox/ddesktop/lab-gehlenborg/seurat/filtered_gene_bc_matrices/hg19", 
                                data_name="pbmc", 
                                min_cells=0, 
                                min_features=0) #32738 x 2700
saveRDS(data_pbmc_full, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_full.rds")


#filtered dataset
data_pbmc_filtered <- read_10x_data(directory="~/Dropbox/ddesktop/lab-gehlenborg/seurat/filtered_gene_bc_matrices/hg19", 
                                    data_name="pbmc", 
                                    min_cells=100, 
                                    min_features=500) #4662 x 2482
data_pbmc_filtered <- subset(data_pbmc_filtered, subset=percent.mt<5) 
dim(data_pbmc_filtered) #4662 x 2439


#analyze data for vitessce visualization
data_pbmc_results <- analyze_data(data_pbmc_filtered)
saveRDS(data_pbmc_results, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_pbmc_results.rds")



# t-cells CD8 -------------------------------------------------------------

#full dataset
data_tcellcd8_full <- read_10x_data(directory="~/Dropbox/ddesktop/lab-gehlenborg/datasets/tcellcd8/filtered_matrices_mex/hg19", 
                                data_name="tcellcd8", 
                                min_cells=0, 
                                min_features=0) #32738 x 10209
saveRDS(data_tcellcd8_full, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_full.rds")


#filtered dataset
data_tcellcd8_filtered <- read_10x_data(directory="~/Dropbox/ddesktop/lab-gehlenborg/datasets/tcellcd8/filtered_matrices_mex/hg19", 
                                    data_name="tcellcd8", 
                                    min_cells=100, 
                                    min_features=500) #6171 x 7856
data_tcellcd8_filtered <- subset(data_tcellcd8_filtered, subset=percent.mt<5) 
dim(data_tcellcd8_filtered) #6171 x 7850


#analyze data for vitessce visualization
data_tcellcd8_results <- analyze_data(data_tcellcd8_filtered)
saveRDS(data_tcellcd8_results, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd8_results.rds")



# t-cells CD4 -------------------------------------------------------------

#full dataset
data_tcellcd4_full <- read_10x_data(directory="~/Dropbox/ddesktop/lab-gehlenborg/datasets/tcellcd4/filtered_matrices_mex/hg19", 
                                    data_name="tcellcd4", 
                                    min_cells=0, 
                                    min_features=0) #32738 x 11213
saveRDS(data_tcellcd4_full, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd4_full.rds")


#filtered dataset
data_tcellcd4_filtered <- read_10x_data(directory="~/Dropbox/ddesktop/lab-gehlenborg/datasets/tcellcd4/filtered_matrices_mex/hg19", 
                                        data_name="tcellcd4", 
                                        min_cells=100, 
                                        min_features=500) #5827 x 7279
data_tcellcd4_filtered <- subset(data_tcellcd4_filtered, subset=percent.mt<5) 
dim(data_tcellcd4_filtered) #5827 x 7276


#analyze data for vitessce visualization
data_tcellcd4_results <- analyze_data(data_tcellcd4_filtered)
saveRDS(data_tcellcd4_results, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_tcellcd4_results.rds")



# lung -------------------------------------------------------------

#full dataset
load("~/Dropbox/ddesktop/lab-gehlenborg/datasets/lung/droplet_normal_lung_blood_seurat_ntiss10x.P1.anno.20191002.RC4.Robj")
data_lung_full <- UpdateSeuratObject(ntiss10x.P1.anno) #update "old seurat object" #26485 x 9744
saveRDS(data_lung_full, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_lung_full.rds")


#filtered dataset
data_lung_filtered <- CreateSeuratObject(counts=data_lung_full@assays$RNA@counts, 
                                         project="lung", min.cells=100, min.features=500) #12514 x 9744
#quality control: detect amount of mitochondrial DNA
data_lung_filtered[["percent.mt"]] <- PercentageFeatureSet(data_lung_filtered, pattern = "^MT-")
data_lung_filtered <- subset(data_lung_filtered, subset=percent.mt<5) 
dim(data_lung_filtered) #12514 x 9744


#analyze data for vitessce visualization
data_lung_results <- analyze_data(data_lung_filtered)
saveRDS(data_lung_results, file="~/Dropbox/ddesktop/lab-gehlenborg/data/data_lung_results.rds")

