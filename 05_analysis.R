library(dplyr)


#cell leakage 用背景表达水平（其他细胞类型的平均值）除以目标细胞类型中的表达值
calc_leakage=function(
    obj,
    gene,
    target_type
){
  
  exp=
    FetchData(
      obj,
      vars=c(
        gene,
        "cell_type"
      )
    )
  
  target=
    mean(
      exp[
        exp$cell_type==
          target_type,
        gene
      ]
    )
  
  off=
    mean(
      exp[
        exp$cell_type!=
          target_type,
        gene
      ]
    )
  
  score=
    off/target
  
  return(score)
  
}



## specificity用目标细胞类型中的基因表达值除以其他细胞类型中的表达值，得到目标基因在目标细胞类型上的特异性比例
specificity=function(
    obj,
    gene,
    target
){
  
  exp=
    FetchData(
      obj,
      vars=c(
        gene,
        "cell_type"
      )
    )
  
  tar=
    mean(
      exp[
        exp$cell_type==target,
        gene
      ])
  
  other=
    mean(
      exp[
        exp$cell_type!=target,
        gene
      ])
  
  tar/other
  
}






#run
calc_leakage(pbmc_raw, "CD3D", "CD4+ T")
calc_leakage(pbmc_raw_DX, "CD3D", "CD4+ T")
calc_leakage(pbmc_SX, "CD3D", "CD4+ T")

specificity(pbmc_raw, "CD3D", "CD4+ T")
specificity(pbmc_raw_DX, "CD3D", "CD4+ T")
specificity(pbmc_SX, "CD3D", "CD4+ T")

# Leakage for Classical Monocytes
calc_leakage(pbmc_raw, "CD14", "CD14+ Monocytes")
calc_leakage(pbmc_raw_DX, "CD14", "CD14+ Monocytes")
calc_leakage(pbmc_SX, "CD14", "CD14+ Monocytes")

# Specificity for Classical Monocytes
specificity(pbmc_raw, "CD14", "CD14+ Monocytes")
specificity(pbmc_raw_DX, "CD14", "CD14+ Monocytes")
specificity(pbmc_SX, "CD14", "CD14+ Monocytes")

# Leakage for Cytotoxic T cells
calc_leakage(pbmc_raw, "CD8A", "CD8+ T")
calc_leakage(pbmc_raw_DX, "CD8A", "CD8+ T")
calc_leakage(pbmc_SX, "CD8A", "CD8+ T")

# Specificity for Cytotoxic T cells
specificity(pbmc_raw, "CD8A", "CD8+ T")
specificity(pbmc_raw_DX, "CD8A", "CD8+ T")
specificity(pbmc_SX, "CD8A", "CD8+ T")

# Leakage for NK cells
calc_leakage(pbmc_raw, "NCAM1", "NK")
calc_leakage(pbmc_raw_DX, "NCAM1", "NK")
calc_leakage(pbmc_SX, "NCAM1", "NK")

# Specificity for NK cells
specificity(pbmc_raw, "NCAM1", "NK")
specificity(pbmc_raw_DX, "NCAM1", "NK")
specificity(pbmc_SX, "NCAM1", "NK")

# Leakage for DC
calc_leakage(pbmc_raw, "FCER1A", "DC")
calc_leakage(pbmc_raw_DX, "FCER1A", "DC")
calc_leakage(pbmc_SX, "FCER1A", "DC")

# Specificity for DC
specificity(pbmc_raw, "FCER1A", "DC")
specificity(pbmc_raw_DX, "FCER1A", "DC")
specificity(pbmc_SX, "FCER1A", "DC")

#load true object (filtered data)
filteredData <- Read10X_h5("./10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5")
true_obj <- CreateSeuratObject(filteredData, min.cells=3,
                               min.features=200)
true_obj<- NormalizeData(true_obj)
true_obj<- FindVariableFeatures(true_obj)
true_obj <- ScaleData(true_obj)
true_obj<- RunPCA(true_obj)
true_obj<- FindNeighbors(
  true_obj,
  dims=1:20
)

true_obj <- FindClusters(
  true_obj,
  resolution=.5
)

true_obj <- RunUMAP(
  true_obj,
  dims=1:20
)



p3 <- DotPlot(
  pbmc_raw,
  features=unique(
    unlist(
      strsplit(
        pbmc_markers$markers,
        ", "
      )
    )
  ),
  group.by = "cell_type"
)+
  RotatedAxis()
p4 <- DotPlot(
  pbmc_raw_DX,
  features=unique(
    unlist(
      strsplit(
        pbmc_markers$markers,
        ", "
      )
    )
  ),
  group.by = "cell_type"
)+
  RotatedAxis()
p5 <- DotPlot(
  pbmc_SX,
  features=unique(
    unlist(
      strsplit(
        pbmc_markers$markers,
        ", "
      )
    )
  ),
  group.by = "cell_type"
)+
  RotatedAxis()




calc_leakage(pbmc_raw, "PPBP", "Platelet")
calc_leakage(pbmc_1, "PPBP", "Platelet")
calc_leakage(pbmc_1_DX, "PPBP", "Platelet")
calc_leakage(pbmc_10, "PPBP", "Platelet")
calc_leakage(pbmc_10_DX, "PPBP", "Platelet")
calc_leakage(pbmc_100, "PPBP", "Platelet")
calc_leakage(pbmc_100_DX, "PPBP", "Platelet")

specificity(pbmc_raw, "PPBP", "Platelet")
specificity(pbmc_1, "PPBP", "Platelet")
specificity(pbmc_1_DX, "PPBP", "Platelet")
specificity(pbmc_10, "PPBP", "Platelet")
specificity(pbmc_10_DX, "PPBP", "Platelet")
specificity(pbmc_100, "PPBP", "Platelet")
specificity(pbmc_100_DX, "PPBP", "Platelet")