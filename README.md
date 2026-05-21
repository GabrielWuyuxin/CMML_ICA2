# CMML ICA2 miniproject 8: scRNA-seq Contamination and Noise Analysis

This repository contains R scripts and datasets for analyzing and assessing contamination and noise in single-cell RNA-sequencing (scRNA-seq) data. The main workflows cover data preprocessing, contamination simulation, and post-analysis using tools such as DecontX and SoupX.

---

## **Content Overview**

### **Scripts**
1. **`01_preprocess.R`**  
   - Preprocesses raw scRNA-seq data using the Seurat workflow.  
   - Includes steps such as normalization, feature selection, scaling, and clustering.  
   
2. **`02_simulation.R`**  
   - Introduces synthetic noise into the dataset by injecting random contamination based on specified parameters (e.g., `contam_fraction` and `contam_genes`).  

3. **`03_runDecontX.R`**  
   - Runs contamination removal analysis using the DecontX tool with default settings.

4. **`04_runSoupX.R`**  
   - Applies the SoupX tool to estimate and remove ambient RNA contamination.  
   - Uses SoupX’s automatic contamination rate estimation pipeline.

5. **`05_analysis.R`**  
   - Conducts downstream analysis, including visualization with UMAP and computing evaluation metrics such as **leakage** and **specificity** for marker genes.

### **Datasets**
- **`10k_PBMC_3p_nextgem_Chromium_v3.rds`**  
  - Processed scRNA-seq data from 10k PBMC (Peripheral Blood Mononuclear Cells) generated using 10x Genomics Chromium v3 chemistry.  

### **Additional Files**
- **`.gitattributes`**  
  - Git metadata file for ensuring consistent line endings and repository behavior.

---

## **Usage Instructions**

### **1. Preprocessing Data**
Run the following script to preprocess the raw data:
```bash
Rscript 01_preprocess.R
```

### **2. Simulate Contamination**
Add synthetic noise to the output of preprocessing using the simulation script:
```bash
Rscript 02_simulation.R
```

### **3. Remove Contamination**
- DecontX:
  ```bash
  Rscript 03_runDecontX.R
  ```
- SoupX:
  ```bash
  Rscript 04_runSoupX.R
  ```

### **4. Analyze Results**
Generate UMAP visualization and compute evaluation metrics by running:
```bash
Rscript 05_analysis.R
```

---

## **Key Metrics**
- **Leakage**  
  Ratio of the background expression level of a marker gene (mean expression in non-representative cell types) to its expression in the target cell type.

- **Specificity**  
  Ratio of a marker gene’s expression in the target cell type to its expression across all other cell types.

---

## **Dependencies**
The analysis relies on the following R packages:
- **[Seurat (v5)](https://satijalab.org/seurat/)** for preprocessing and clustering.
- **[celda (DecontX)](https://bioconductor.org/packages/release/bioc/html/celda.html)** for contamination removal based on clustering results.
- **[SoupX](https://github.com/constantAmateur/SoupX)** for ambient RNA decontamination.

Install required packages using:
```R
install.packages(c("Seurat", "SoupX"))
BiocManager::install("celda")
```

---

## **Code and Data Availability**
All datasets and source code in this repository are publicly available. For questions or reporting issues, please create an issue in this GitHub repository.  

---

## **License**
This project is licensed under the MIT License.

---
