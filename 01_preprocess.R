.libPaths(c("C:/rtools44", .libPaths()))
library(Seurat)

preprocess_pbmc=function(raw_path){
  
  rawData <- Read10X_h5(raw_path)
  
  pbmc <- CreateSeuratObject(
    counts=rawData,
    min.cells=3,
    min.features=200
  )
  
  pbmc[["percent.mt"]] <-
    PercentageFeatureSet(
      pbmc,
      pattern="^MT"
    )
  
  pbmc <- subset(
    pbmc,
    subset=
      nFeature_RNA>200 &
      nFeature_RNA<6000 &
      percent.mt<10
  )
  
  pbmc <- NormalizeData(pbmc)
  
  pbmc <- FindVariableFeatures(pbmc)
  
  pbmc <- ScaleData(pbmc)
  
  pbmc <- RunPCA(pbmc)
  
  pbmc <- FindNeighbors(
    pbmc,
    dims=1:20
  )
  
  pbmc <- FindClusters(
    pbmc,
    resolution=.5
  )
  
  pbmc <- RunUMAP(
    pbmc,
    dims=1:20
  )
  
  return(pbmc)
  
}

pbmc_raw=preprocess_pbmc(
  "D:/wuyuxin/CMML/ICA2/10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5"
)




###细胞注释

pbmc_markers <- data.frame(
  cell_type = c(
    "CD4+ T cells",
    "CD14+ Monocytes",
    "CD8+ T cells",
    "NK cells",
    "B cells",
    "CD16+ Monocytes",
    "Dendritic cells",
    "Platelets"
  ),
  
  markers = c(
    
    # CD4 T
    "CD3D, CD3E, CD3G, CD4, CCR7",
    
    # classical mono
    "CD14, LYZ",
    
    # cytotoxic T
    "CD3D, CD3E, CD3G, CD8A, GZMB, PRF1",
    
    # NK
    "NCAM1, KLRD1",
    
    # B
    "CD19, MS4A1, CD79A",
    
    # CD16 Mono
    "FCGR3A, LYZ",
    
    # DC
    "FCER1A, CD1C, CLEC9A, IL3RA",
    
    # Platelet
    "ITGA2B, ITGB3, GP9, GP1BA"
    
  ),
  
  stringsAsFactors=FALSE
)


#看marker表达情况

FeaturePlot(
  pbmc_raw,
  features=c(
    "CD3D",
    "CD4",
    "CD8A",
    "LYZ",
    "FCGR3A",
    "NKG7",
    "MS4A1",
    "FCER1A",
    "ITGA2B"
  ),
  ncol=3
)


DotPlot(
  pbmc_raw,
  features=unique(
    unlist(
      strsplit(
        pbmc_markers$markers,
        ", "
      )
    )
  ),
  group.by = "seurat_clusters"
)+
  RotatedAxis()



#根据marker人工注释cluster

cluster2type <- c(
  
  "0"="CD14+ Monocytes",
  "1"="CD4+ T",
  "2"="CD8+ T",
  "3"="CD4+ T",
  "4"="CD14+ Monocytes",
  "5"="CD8+ T",
  "6"="B",
  "7"="B",
  "8"="CD16+ Monocytes",
  "9"="NK",
  "10"="CD14+ Monocytes",
  "11"="DC",
  "12"="pDC",
  "13"="Platelet"
  
)

pbmc_raw@meta.data$cell_type <- cluster2type[as.character(pbmc_raw@meta.data$seurat_clusters)]


Idents(pbmc_raw) <-
  pbmc_raw$cell_type


DimPlot(
  pbmc_raw,
  reduction="umap",
  group.by="cell_type",
  label=TRUE
)


saveRDS(
  pbmc_raw,
  "pbmc_raw.rds"
)
