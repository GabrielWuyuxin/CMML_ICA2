#devtools::install_github("campbio/decontX")
#.libPaths(c("C:/rtools44", .libPaths()))
#library(Seurat)
library(decontX)
#library(dplyr)
#library(ggplot2)
library(celda)


# decontX

pbmc_1 <- readRDS("./pbmc_1.rds")
pbmc_10 <- readRDS("./pbmc_10.rds")
pbmc_100 <- readRDS("./pbmc_100.rds")

run_decontX <- function(seurat_obj){
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- decontX(sce)
  clean_counts <- decontXcounts(sce)
  seurat_obj$decontX_contamination <- colData(sce)$decontX_contamination
  
  seurat_obj[["decontX"]] <-
    CreateAssayObject(
      counts=clean_counts
    )
  
  
  #create new seurat object
  obj_decontX <- CreateSeuratObject(
    counts=clean_counts,
    meta.data=seurat_obj@meta.data
  )
  Idents(obj_decontX) <- Idents(seurat_obj)
  
  
  #run all downstream again
  obj_decontX <- NormalizeData(obj_decontX)
  
  obj_decontX <- FindVariableFeatures(obj_decontX)
  
  obj_decontX <- ScaleData(obj_decontX)
  
  obj_decontX <- RunPCA(obj_decontX)
  
  obj_decontX <- FindNeighbors(obj_decontX,dims=1:20
  )
  
  obj_decontX <- FindClusters(obj_decontX,resolution=.5
  )
  
  obj_decontX <- RunUMAP(obj_decontX,dims=1:20
  )
  
  return(list(seurat_obj = seurat_obj, 
              obj_decontX = obj_decontX))
  
}

#raw data
res_decontX <- run_decontX(pbmc_raw)
pbmc_raw <- res_decontX[[1]]
pbmc_raw_DX <- res_decontX[[2]]


#pbmc 1
res_1 <- run_decontX(pbmc_1)
pbmc_1      <- res_1[[1]]     
pbmc_1_DX   <- res_1[[2]]


#pbmc 10
# res_10 <- run_decontX(pbmc_10)
# pbmc_10      <- res_10[[1]]     
# pbmc_10_DX   <- res_10[[2]]  


#pbmc 100
res_100 <- run_decontX(pbmc_100)
pbmc_100      <- res_100[[1]]     
pbmc_100_DX   <- res_100[[2]]



  
DimPlot(pbmc_1_DX, 
        group.by  = "seurat_clusters")
DimPlot(pbmc_1_DX,
        group.by = "cell_type")
DimPlot(pbmc_1, group.by = "cell_type")




p1=DimPlot(pbmc_100, 
           reduction = "umap", 
           group.by  = "cell_type")+ggtitle("100 gene contamination")

p2=DimPlot(pbmc_100_DX, 
           reduction = "umap", 
           group.by  = "cell_type")+ggtitle("DecontX")

p1+p2










summary(pbmc_raw$decontX_contamination)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000415 0.0041834 0.0125866 0.0608095 0.0637045 0.9976453 
summary(pbmc_1$decontX_contamination)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000463 0.0047695 0.0148229 0.0647294 0.0686101 0.9971398
summary(pbmc_10$decontX_contamination)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000249 0.0040472 0.0114201 0.0626560 0.0591518 0.9978486 
summary(pbmc_100$decontX_contamination)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000364 0.0089631 0.0244452 0.0692973 0.0772402 0.9934024 


#查看decontX contamination
FeatureScatter(
  pbmc_100, 
  feature1="nCount_RNA",
  feature2="decontX_contamination",
  group.by="cell_type"
)

VlnPlot(
  pbmc_10,
  features="decontX_contamination",
  group.by="cell_type"
)
hist(
  pbmc_raw$decontX_contamination,
  breaks=100
)



