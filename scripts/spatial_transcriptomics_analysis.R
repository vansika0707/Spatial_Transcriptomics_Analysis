install.packages("devtools")
devtools::install_github("satijalab/seurat-data")
install.packages('BiocManager')
BiocManager::install('glmGamPoi')
install.packages('Rfast2')
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratData)
InstallData("pbmc3k")
library("htmltools")

InstallData("stxBrain")
brain <- LoadData('stxBrain', type = 'anterior1')
## Data preprocessing
#The initial preprocessing steps that we perform on the spot by gene expression data 
#are similar to a typical scRNA-seq experiment. 
#We first need to normalize the data in order to account for variance 
#in sequencing depth across data points. 

plot1 <- VlnPlot(brain, features = "nCount_Spatial", layer = "counts") + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = 'nCount_Spatial') + theme(legend.position = "right")
wrap_plots(plot1, plot2)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# normalization
# rerun normalization to store sctransform residuals for all genes
brain <- SCTransform(brain, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# also run standard log normalization for comparison
brain <- NormalizeData(brain, verbose = FALSE, assay = "Spatial")
## norm test 2
# Computes the correlation of the log normalized data and sctransform residuals with the number of UMIs
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
# norm test 3
p1 <- GroupCorrelationPlot(brain, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") + theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(brain, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") + theme(plot.title = element_text(hjust = 0.5))
p1 + p2
## Gene expression visualization
# Visualize Gene Expression
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
#Dimensionality reduction + Clustering
brain <- RunPCA(brain, assay = "SCT")
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
DimPlot(brain, reduction = "umap")
SpatialDimPlot(brain, label = TRUE)
# Identify Spatially Variable Features
head(GetTissueCoordinates(brain))
coords <- GetTissueCoordinates(brain)
brain <- LoadData("stxBrain", type = "anterior1")
brain <- FindSpatiallyVariableFeatures(
  brain,
  assay = "SCT",
  selection.method = "moransi",
  features = VariableFeatures(brain)[1:1000]
)
#spatially variable gene expression across the brain tissue section
top.features <- head(SpatiallyVariableFeatures(brain), 6)
SpatialFeaturePlot(brain, features = top.features)
#cortex-enriched clusters
cortex <- subset(brain, idents = c(1,2,3,4,6,7))
SpatialDimPlot(cortex, crop = TRUE)

