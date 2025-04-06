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
library(SingleR)
library(celldex)
library(Seurat)
library(SeuratObject)
library(pheatmap)



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

#Annotation of celltypes using SingleR (Source- Blurprint Encode project)
Blueprint<- celldex::BlueprintEncodeData()

#Extract counts data and make SCE dataset
hnscc.KVG.counts <- GetAssayData(hnscc.integrated, assay = "integrated", layer = "counts")
hnscc.KVG.counts <- as.SingleCellExperiment(hnscc.integrated)

#Annotate
pred <- SingleR(test = hnscc.KVG.counts,
                ref = Blueprint,
                labels = Blueprint$label.main)

#update seurat metadata with annotation data
hnscc.integrated$blueprint.main <- pred$labels[match(rownames(hnscc.integrated@meta.data), rownames(pred))]
UMAP_SingleR <- DimPlot(hnscc.integrated, reduction = "umap", group.by = "blueprint.main") + ggtitle("SingleR")

#change idents for further analysis
Idents(hnscc.integrated) <- "blueprint.main"

#verify the correctness of Blueprint encode based annotation
#identify celltypes assigned and markers used for annotation

celltypes_blueprint <- unique(hnscc.integrated@meta.data$blueprint.main)
markers_blueprint <- FindAllMarkers(hnscc.integrated)

top_markers_bp <- markers_blueprint %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%   # Select top 20 based on avg_log2FC
  arrange(cluster, desc(avg_log2FC)) 

#identify top markers for each cells and store it for ease of access
Cells <- c("Fibroblast", "Keratinocytes", "Endothelial cell")
Hits <- c(5, 12, 9)
Genes <- c("THY1, PCOLCE, COL5A3, GREM1, MMP11",
           "KRT14, KRT17, KRT8, PKP3, S100A2, EPPK1, TP63, CALML3, LAMC2, LAMA3, LAMB3, COL17A1",
           "CDH5, KDR, CLDN5, TIE1, ROBO4, CLEC14A, DLL4, PLVAP, and ANGPT2",
)

markers_bp <- data.frame(Cells = Cells, Hits = Hits, Genes = Genes)

#visualize marker gene expression in UMAP to confirm presence in corresponsing cells
FeaturePlot(hnscc.integrated, features = c("THY1","PCOLCE", "COL5A3", "GREM1", "MMP11"), min.cutoff = "q10", max.cutoff = "q90")
