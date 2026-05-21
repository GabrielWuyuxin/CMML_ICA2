#devtools::install_github("constantAmateur/SoupX", ref = "devel")
library(Seurat)
library(SoupX)
library(dplyr)
library(ggplot2)
library(celda)
library(tidyverse)

tod <- Read10X_h5("D:/wuyuxin/CMML/ICA2/10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5")
toc <- Read10X_h5("D:/wuyuxin/CMML/ICA2/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5")
sc <- SoupChannel(tod, toc)

##常规聚类和降维
toc <- CreateSeuratObject(toc)
toc <- NormalizeData(toc, normalization.method = "LogNormalize", scale.factor = 10000)
toc <- FindVariableFeatures(toc, selection.method = "vst", nfeatures = 3000)
toc.genes <- rownames(toc)
toc <- ScaleData(toc, features = toc.genes)
toc <- RunPCA(toc, features = VariableFeatures(toc), npcs = 40, verbose = F)
toc <- FindNeighbors(toc, dims = 1:30)
toc <- FindClusters(toc, resolution = 0.5)
toc <- RunUMAP(toc, dims = 1:30)


#cell annotation
FeaturePlot(
  toc,
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
    
    # CD14+ mono
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



DotPlot(
  toc,
  features=unique(
    unlist(
      strsplit(
        pbmc_markers$markers,
        ", "
      )
    )
  )
)+
  RotatedAxis()

cluster2type <- c(
  
  "0"="CD14+ Monocytes",
  "1"="CD4+ T",
  "2"="CD8+ T",
  "3"="CD4+ T",
  "4"="CD14+ Monocytes",
  "5"="CD8+ T",
  "6"="B",
  "7"="B",
  "8"="NK",
  "9"="CD14+ Monocytes",
  "10"="unknown",
  "11"="DC",
  "12"="pDC",
  "13"="unknown",
  "14"="Platelet"
  
)

toc@meta.data$cell_type <- cluster2type[as.character(toc@meta.data$seurat_clusters)]


Idents(toc) <-
  toc$cell_type

DimPlot(toc, group.by = "cell_type")

#提取聚类后的meta.data信息
matx <- toc@meta.data
sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
head(sc$metaData)

nonExpressedGeneList = list(HB = c("HBB","HBA2"), IG = c("IGKC"))

#自动生成目标污染率
sc = autoEstCont(sc)


#手动设置目标污染率
sc <- setContaminationFraction(sc, 0.1)


sc=adjustCounts(sc)


saveRDS(sc, "sc.rds")


sc <- readRDS("sc.rds")
class(sc)
# 假设您之前创建了 SoupChannel sc
# 去噪表达矩阵

#out=adjustCounts(sc)
pbmc_SX <- CreateSeuratObject(counts = sc)
# 执行降维（如 PCA、UMAP）
pbmc_SX <- NormalizeData(pbmc_SX)
pbmc_SX <- FindVariableFeatures(pbmc_SX)
pbmc_SX <- ScaleData(pbmc_SX)
pbmc_SX <- RunPCA(pbmc_SX)
pbmc_SX <- RunUMAP(pbmc_SX, dims = 1:20)

meta1 <- toc@meta.data  # 从 Seurat Object 1 提取 metadata
meta2 <- pbmc_SX@meta.data  # 从 Seurat Object 2 提取 metadata

# 使用行名作为映射键，将 Seurat Object 1 的 meta.data 信息添加到 Seurat Object 2
meta2$seurat_clusters <- meta1$seurat_clusters[match(rownames(meta2), rownames(meta1))]
meta2$cell_type <- meta1$cell_type[match(rownames(meta2), rownames(meta1))]

# 将更新后的 meta.data 写回到 Seurat Object 2
pbmc_SX@meta.data <- meta2

# DimPlot 可视化
DimPlot(pbmc_SX, reduction = "umap", group.by = "cell_type")
