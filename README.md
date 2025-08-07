# Spatial Transcriptomics Analysis - Mouse Brain (Sagittal-Anterior)

This project replicates the **Spatial Transcriptomics analysis** from [Satija Lab's Seurat vignette](https://github.com/satijalab/seurat/blob/master/vignettes/spatial_vignette.Rmd) using the **Mouse Brain Serial Section 1 (Sagittal-Anterior)** dataset from **10x Genomics**.

## Dataset
Dataset used:
Mouse Brain Serial Section 1 (Sagittal-Anterior)

Download via SeuratData:
library(SeuratData)
InstallData("stxBrain")


## Required R packages
Install all dependencies by running:
```r
install.packages("devtools")
devtools::install_github("satijalab/seurat")
devtools::install_github("satijalab/seurat-data")
install.packages(c("ggplot2", "patchwork", "dplyr", "htmltools", "Rfast2"))
install.packages("BiocManager")
BiocManager::install("glmGamPoi")
```


## Summary of Analysis
1. Quality Control Plots
Violin plot and spatial map of UMI (Unique Molecular Identifier) counts per spot, show spatial heterogeneity indicating different cell densities.

2. Normalization Comparision 
Log normalization vs SCTransform. The normalization comparison plot shows how normalization reduces variability due to sequencing depth or UMI counts, making gene expression across spatial spots more comparable and reliable for downstream analysis.

3. Spatially Variable gene expression
The plot displays spatially variable gene expression of Hpca(Hippocalcin) and Ttr(Transthyretin) in the anterior cortex of the brain.

4. Clustering and Dimesionality Reduction 
PCA + UMAP is used to identify transcriptionally distinct clusters.Spatially mapped clusters reveal correspondence to anatomical brain structures.

5. Cortex subset
Visualization of subset of clusters enriched in cortex region.


## Conclusion 
Spatial transcriptomics is a powerful technique that combines gene expression profiling with spatial information, allowing researchers to map where genes are active within a tissue.
This project demonstrated the application of spatial transcriptomics to study gene expression in the mouse brain, specifically focusing on the sagittal-anterior section. By integrating transcriptomic data with spatial information using the Seurat package, we visualized spatially variable gene expression patterns and identified region-specific markers(Hpca and Ttr).


## Acknowledgements
Satija Lab
10x Genomics dataset
Seurat R package

