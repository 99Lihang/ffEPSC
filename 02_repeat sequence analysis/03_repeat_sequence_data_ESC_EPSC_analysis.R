library(Seurat)
library(ggplot2)
library(data.table)

# =======================================
# Subset data for ffEPS cells and normalize

EPSC_cell_rpmk <- subset(seurat_obj,subset = cell_type%in%c("ffEPS"))

EPSC_cell_rpmk <- NormalizeData(EPSC_cell_rpmk, normalization.method = "LogNormalize", scale.factor = 10000)
EPSC_cell_rpmk <- FindVariableFeatures(EPSC_cell_rpmk, selection.method = "vst", nfeatures = 1500)
# Scale the data using all genes

# Scale the data using all genes
all.genes <- rownames(EPSC_cell_rpmk)
EPSC_cell_rpmk <- ScaleData(EPSC_cell_rpmk, features = all.genes)

# Perform PCA (determine optimal number of PCs using ElbowPlot)
EPSC_cell_rpmk <- RunPCA(EPSC_cell_rpmk, features = VariableFeatures(object = EPSC_cell_rpmk), npcs = 24)

# Find neighbors and clusters
EPSC_cell_rpmk <- FindNeighbors(EPSC_cell_rpmk, dims = 1:24)
EPSC_cell_rpmk <- FindClusters(EPSC_cell_rpmk, resolution = 1)

# Perform dimensionality reduction using t-SNE or UMAP
EPSC_cell_rpmk <- RunTSNE(EPSC_cell_rpmk,perplexity =8)
EPSC_cell_rpmk <- RunUMAP(EPSC_cell_rpmk,dims = 1:20,n.neighbors = 24)

EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_1"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_11"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_13"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_9"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_2"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_6"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_25"),"cell_subtype"]<-"EPS.1"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_3"),"cell_subtype"]<-"EPS.1"

EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_21"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_19"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_24"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_17"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_14"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_5"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_16"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_18"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_23"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_15"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_20"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_10"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_4"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_12"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_8"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_22"),"cell_subtype"]<-"EPS.2"
EPSC_cell_rpmk@meta.data[which(EPSC_cell_rpmk@meta.data$cell_name=="ffEPS_7"),"cell_subtype"]<-"EPS.2"

Idents(EPSC_cell_rpmk) <- EPSC_cell_rpmk@meta.data$cell_subtype
EPS.markers <- FindAllMarkers(EPSC_cell_rpmk,logfc.threshold = 0.1,only.pos = T)



intersected_genes <- c("LTR78","LTR5","LTR5-Hs","LTR5A","LTR62","LTR17","LTR06","LTR14B","LTR12C")


for (gene in intersected_genes) {
  # Extract expression data for the specific gene
  gene_data <- FetchData(EPSC_cell_rpmk, vars = c(gene, "cell_subtype"))
  
  # Print the plot to the PDF
  boxplot_plot <- ggplot(gene_data, aes(x = cell_subtype, y = get(gene), fill = cell_subtype)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = gene, y = "Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Print the plot to the PDF
  print(boxplot_plot)
}






# =======================================
# Subset data for ESC cells and normalize
ESC_cell_rpmk <- subset(seurat_obj, subset = cell_type %in% c("ESC"))

ESC_cell_rpmk <- NormalizeData(ESC_cell_rpmk, normalization.method = "LogNormalize", scale.factor = 10000)
ESC_cell_rpmk <- FindVariableFeatures(ESC_cell_rpmk, selection.method = "vst", nfeatures = 1500)
# Scale the data using all genes
all.genes <- rownames(ESC_cell_rpmk)
ESC_cell_rpmk <- ScaleData(ESC_cell_rpmk, features = all.genes)

# Perform PCA (determine optimal number of PCs using ElbowPlot)
ESC_cell_rpmk <- RunPCA(ESC_cell_rpmk, features = VariableFeatures(object = ESC_cell_rpmk), npcs = 23)


# Find neighbors and clusters
ESC_cell_rpmk <- FindNeighbors(ESC_cell_rpmk, dims = 1:20)
ESC_cell_rpmk <- FindClusters(ESC_cell_rpmk, resolution = 1.03)


# Perform dimensionality reduction using t-SNE or UMAP
ESC_cell_rpmk <- RunTSNE(ESC_cell_rpmk, dims = 1:15, perplexity = 5)
ESC_cell_rpmk <- RunUMAP(ESC_cell_rpmk, dims = 1:15, n.neighbors = 16)


