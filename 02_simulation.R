library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
setwd("D:/wuyuxin/CMML/ICA2")
#read in data
pbmc_raw <- readRDS("pbmc_raw.rds")
#class(pbmc_raw)

simulate_contamination=function(
    seurat_obj,
    contam_genes,
    contam_fraction=.1
){
  
  counts=GetAssayData(
    seurat_obj,
    layer="counts"
  )
  
  counts <- as(
    counts,
    "dgCMatrix"
  )
  
  n_add=
    round(
      colSums(counts)*
        contam_fraction
    )
  
  for(g in contam_genes){
    
    if(g %in% rownames(counts)){
      
      counts[g,]=
        counts[g,]+
        rpois(
          ncol(counts),
          lambda=
            n_add/
            length(contam_genes)
        )
      
    }
    
  }
  
  new_obj=
    CreateSeuratObject(
      counts,
      meta.data = seurat_obj@meta.data
    )
  Idents(new_obj) <- Idents(seurat_obj)
  
  #go through pipeline again
  new_obj <- NormalizeData(new_obj)
  
  new_obj <- FindVariableFeatures(new_obj)
  
  new_obj <- ScaleData(new_obj)
  
  new_obj <- RunPCA(new_obj)
  
  new_obj <- FindNeighbors(new_obj,dims=1:20
  )
  
  new_obj <- FindClusters(new_obj,resolution=.5
  )
  
  new_obj <- RunUMAP(new_obj,dims=1:20
  )
  return(new_obj)
  
}



#run with noised data

pbmc_1=
  simulate_contamination(
    pbmc_raw,
    contam_genes="PPBP",
    contam_fraction=.1
  )
head(pbmc_1)

top10=c(
  "PPBP",
  "PF4",
  "GP9",
  "ITGA2B",
  "NRGN",
  "TUBB1",
  "SPARC",
  "RGS18",
  "CLU",
  "GNG11"
)

pbmc_10=
  simulate_contamination(
    pbmc_raw,
    top10,
    .1
  )

#100gene
platelet_cells=
  WhichCells(
    pbmc_raw,
    ident="Platelet"
  )

avg=
  rowMeans(
    GetAssayData(pbmc_raw)[
      ,
      platelet_cells
    ]
  )

top100=
  names(
    sort(
      avg,
      decreasing=T
    )
  )[1:100]
head(top100)
#[1] "TMSB4X" "FTL"    "B2M"    "ACTB"   "FTH1"   "MT-CO1"

pbmc100=
  simulate_contamination(
    pbmc_raw,
    top100,
    .1
  )



saveRDS(pbmc_1, "pbmc_1.rds")
saveRDS(pbmc_10, "pbmc_10.rds")
saveRDS(pbmc100, "pbmc_100.rds")

DimPlot(pbmc_raw, group.by = "cell_type")
DimPlot(pbmc_1, group.by = "cell_type")
DimPlot(pbmc_10, group.by = "cell_type")
DimPlot(pbmc100, group.by = "cell_type")
FeaturePlot(pbmc_raw, features = "PPBP")
FeaturePlot(pbmc_1, features = "PPBP")
FeaturePlot(pbmc_10, features = "PPBP")
FeaturePlot(pbmc100, features = "PPBP")
