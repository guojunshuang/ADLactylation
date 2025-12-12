library(Seurat)
library(jsonlite)
library(png)
library(tidyverse)
library(ggpubr)
library(patchwork)
setwd("~/singlecell/Lactylation/AD/GSE220442/1-1")
spe1 = Read10X("./filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("~/singlecell/Lactylation/AD/GSE220442/1-1", 
                                              "spatial"), filter.matrix = TRUE)
spe1 <- CreateSeuratObject(counts = spe1, assay = "Spatial")
image2 <- image2[Cells(x = spe1)]
DefaultAssay(spe1 = image2) <- "Spatial"
spe1[["slice1"]] <- image2
for (i in colnames((spe1@images$slice1@coordinates))) {
  spe1@images$slice1@coordinates[[i]] <- as.integer(spe1@images$slice1@coordinates[[i]])
}
SpatialFeaturePlot(spe1, features = "nFeature_Spatial")
spe1@meta.data<-cbind(spe1@meta.data,image2@coordinates)
spe1 <- SCTransform(spe1, assay = "Spatial", verbose = FALSE)
spe1 <- RunPCA(spe1, assay = "SCT", verbose = FALSE)
spe1 <- FindNeighbors(spe1, reduction = "pca", dims = 1:30)
spe1 <- FindClusters(spe1, verbose = FALSE)
spe1 <- RunUMAP(spe1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(spe1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(spe1, label = TRUE, label.size = 3)
p1 + p2
SpatialDimPlot(spe1, 
               cells.highlight = CellsByIdentities(object = spe1, 
                                                   idents = c(0,1,2,3,4,5,6,7)), 
               facet.highlight = TRUE, ncol = 4)
de_markers <- FindMarkers(spe1, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = spe1, features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
SpatialFeaturePlot(spe1, features = c("ARF1"))
SpatialFeaturePlot(spe1, features = c("ACSF2"))
SpatialFeaturePlot(spe1, features = c("LDHA"))
SpatialFeaturePlot(spe1, features = c("MECP2"))
setwd("~/singlecell/Lactylation/AD/GSE220442/1-1")
save(spe1,file = '1_1.rdata')
load('~/singlecell/Lactylation/AD/GSE220442/1-1/1_1.rdata')
ldha_expression <- GetAssayData(spe1, assay = "SCT", slot = "data")["LDHA", ]
ldha_df <- data.frame(Cell = names(ldha_expression), LDHA_Expression = as.numeric(ldha_expression))
write.table(ldha_df, file = "1_1LDHA_expression.txt", sep = "\t", row.names = FALSE, quote = FALSE)
