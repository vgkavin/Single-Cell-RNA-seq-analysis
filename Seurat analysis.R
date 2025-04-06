#Downloaded the seurat object from public source.
#updated the seurat object to latest version
#subsetted the raw counts and features matrix data from the integrated seurat object
#performed SC transform and integration

install.packages("Seurat")
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggplot2)





hnscc.int <- readRDS("hnscc.gene.expression.integrated.rds")

hnscc.int <- UpdateSeuratObject(hnscc.int)

View(hnscc.int@meta.data)


cts <- hnscc.int@assays #subset assays from seurat object
cts <- cts[["RNA"]] #subset the SCT transform slot
hnscc.int.seurat.obj <- CreateSeuratObject(cts) #make a seurat object with RNA (Raw data) slot

View(hnscc.int.seurat.obj@meta.data)

#edit metadata
hnscc.int.seurat.obj$sample <- rownames(hnscc.int.seurat.obj@meta.data)
hnscc.int.seurat.obj@meta.data <- separate(hnscc.int.seurat.obj@meta.data, col = 'sample', into = c('details', 'barcode'), 
                                           sep = '_')
hnscc.int.seurat.obj@meta.data <- separate(hnscc.int.seurat.obj@meta.data, col = 'details', into = c('patient', 'treat'), 
                                           sep = '\\.')
View(hnscc.int.seurat.obj@meta.data)

#Data integration steps
#split the by patient numbers
obj.list <- SplitObject(hnscc.int.seurat.obj, split.by = 'patient')


#perform normalization and find DEGs for the split up dataset
for(i in 1:length(obj.list)){
  obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = TRUE)
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], nfeatures = 3000, verbose = TRUE )
}

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

#Find anchors
anchors <- FindIntegrationAnchors(object = obj.list,
                                  anchor.features = features)

hnscc.integrated <- IntegrateData(anchorset = anchors)


#Downstream analysis with integrated data

# Scaling- To reduce/normalize unwanted noise/variation- 
#Helps reduce errors while clustering

all_genes <- rownames(hnscc.integrated)
hnscc.integrated <- ScaleData(hnscc.integrated, features = all_genes)

# Linear Dimensionality reduction- Princple component analysis

DEGs <- VariableFeatures(hnscc.integrated)
hnscc.integrated <- RunPCA(hnscc.integrated, features = DEGs)

#to view the top 5 positive and negative sources of heterogeniety in all 5 principle components
print(hnscc.integrated[["pca"]], dims = 1:20, nfeatures = 5) 

#to visualize the sources of heterogeinety in each cell.
DimHeatmap(hnscc.integrated, dims = 1, cells = 500, balanced = TRUE) 
VizDimLoadings(hnscc.integrated, dims = 1:2, reduction = "pca")

#to determine the dimensionality and filter PCs based on variance.
ElbowPlot(hnscc.integrated) 

#Clustering

#First, identify cells that have similar gene(feature) expression pattern
hnscc.integrated <- FindNeighbors(hnscc.integrated, dims = 1:50)

#Cluster those cells with similar expression

hnscc.integrated <- FindClusters(hnscc.integrated, resolution = c(seq(0, 1, by = 0.01))) #play around with resolution to find the best fit
View(hnscc.integrated@meta.data)
PCA_HNSCC <- DimPlot(hnscc.integrated, reduction = "pca")

Idents(hnscc.integrated) <- "integrated_snn_res.0.29"
Idents(hnscc.integrated)

#non linear dimensionality reduction using umap

hnscc.integrated <- RunUMAP(hnscc.integrated, dims = 1:50)
UMAP_HNSCC <- DimPlot(hnscc.integrated, reduction = "umap")

grid.arrange(PCA_HNSCC, UMAP_HNSCC, ncol=2)

UMAP
saveRDS(hnscc.integrated, "hnscc_KVG.rds")

