library(Seurat)

# function to analyze data ------------------------------------------------
#same function used in data processing

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