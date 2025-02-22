library(Seurat)
library(ggplot2)
library(data.table)

# =======================================
# Subset data for ffEPS cells and normalize
EPSC_cell <- subset(seurat_obj, subset = cell_type %in% c("ffEPS"))
EPSC_cell <- NormalizeData(EPSC_cell, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features (feature selection)
EPSC_cell <- FindVariableFeatures(EPSC_cell, selection.method = "vst", nfeatures = 5000)
# Scale the data using all genes
all.genes <- rownames(EPSC_cell)
EPSC_cell <- ScaleData(EPSC_cell, features = all.genes)

# Perform PCA (determine optimal number of PCs using ElbowPlot)
EPSC_cell <- RunPCA(EPSC_cell, features = VariableFeatures(object = EPSC_cell), npcs = 23)

# Find neighbors and clusters
EPSC_cell <- FindNeighbors(EPSC_cell, dims = 1:20)
EPSC_cell <- FindClusters(EPSC_cell, resolution = 1.03)


# Perform dimensionality reduction using t-SNE or UMAP
EPSC_cell <- RunTSNE(EPSC_cell, dims = 1:15, perplexity = 5)
EPSC_cell <- RunUMAP(EPSC_cell, dims = 1:15, n.neighbors = 10)

# =======================================
# Subset data for ESC cells and normalize
ESC_cell <- subset(seurat_obj, subset = cell_type %in% c("ESC"))
ESC_cell <- NormalizeData(ESC_cell, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features for ESC cells
ESC_cell <- FindVariableFeatures(ESC_cell, selection.method = "vst", nfeatures = 5000)
top10 <- head(VariableFeatures(ESC_cell), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(ESC_cell)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the data using all genes
all.genes <- rownames(ESC_cell)
ESC_cell <- ScaleData(ESC_cell, features = all.genes)

# Perform PCA
ESC_cell <- RunPCA(ESC_cell, features = VariableFeatures(object = ESC_cell), npcs = 15)

# Find neighbors and clusters
ESC_cell <- FindNeighbors(ESC_cell, dims = 1:15)
ESC_cell <- FindClusters(ESC_cell, resolution = 1.03)

# Display clustering results
table(ESC_cell$seurat_clusters)

# Perform dimensionality reduction using t-SNE or UMAP
ESC_cell <- RunTSNE(ESC_cell, dims = 1:15, perplexity = 5)
ESC_cell <- RunUMAP(ESC_cell, dims = 1:15, n.neighbors = 10)


ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_1"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_2"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_3"),"cell_subtype"]<-"ESC.2"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_4"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_5"),"cell_subtype"]<-"ESC.2"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_6"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_7"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_8"),"cell_subtype"]<-"ESC.2"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_9"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_10"),"cell_subtype"]<-"ESC.2"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_11"),"cell_subtype"]<-"ESC.2"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_12"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_13"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_14"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_15"),"cell_subtype"]<-"ESC.1"
ESC_cell@meta.data[which(ESC_cell@meta.data$cell_name=="ESC_16"),"cell_subtype"]<-"ESC.1"


# Visualize clustering with UMAP
DimPlot(ESC_cell, reduction = "umap", group.by = "cell_subtype")
Idents(ESC_cell) <- ESC_cell@meta.data$cell_subtype
ESC.markers <- FindAllMarkers(ESC_cell,logfc.threshold = 0.1)

# =======================================
# Gene set-based visualization
folder_path <- "/geneset/"
plot_path <- "/plot/"

# List all .txt files containing gene sets
txt_files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)

# Extract expression data and metadata
expr_data <- GetAssayData(object = ESC_cell, assay = "RNA", slot = "data")
expr_data <- as.matrix(expr_data)
meta <- data.table(cell = rownames(ESC_cell@meta.data), ESC_cell@meta.data)
setnames(meta, "cell", "cell_name")
new_celltype <- sort(ESC_cell$cell_subtype)

# Iterate through gene set files
for (txt_file in txt_files) {
  # Read gene markers from file
  cellmarkers <- read.table(txt_file, stringsAsFactors = FALSE)
  
  # Ensure genes exist in dataset
  interesting_genes <- rownames(ESC_cell) %in% cellmarkers$V1
  selected_data <- expr_data[interesting_genes, names(new_celltype), drop = FALSE]
  
  # Ensure gene order is maintained
  gene_order <- cellmarkers$V1
  existing_genes <- gene_order[gene_order %in% rownames(selected_data)]
  selected_data <- selected_data[existing_genes, , drop = FALSE]
  
  # Create annotation for cell types
  cell_annotation <- data.frame(cluster = new_celltype)
  
  # Generate output filename
  base_name <- sub("\\.txt$", "", basename(txt_file))
  pdf_filename <- paste0(plot_path, "ESC_", base_name, "_gene.pdf")
  
  # Generate and save DotPlot
  p <- DotPlot(ESC_cell, features = rownames(selected_data), group.by = "cell_subtype") +
    scale_color_gradient(low = "lightgrey", high = "#006995") + coord_flip()
  
  ggsave(pdf_filename, plot = p, width = 8, height = 5)
}
