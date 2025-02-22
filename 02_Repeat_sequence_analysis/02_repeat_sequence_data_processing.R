# Load required libraries
library(Seurat)
library(dplyr)
library(plotrix)
library(parallel)
library(gplots)

# Read counts matrix
counts <- read.table("/home/T2T/counts/counts.txt")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts,
                                 min.features = 200,
                                 min.cells = 3, 
                                 project = "owndata")
### ================================

# Normalize data using log-normalization
EPS_ESC_cell_rpmk <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
EPS_ESC_cell_rpmk <- FindVariableFeatures(EPS_ESC_cell_rpmk, selection.method = "vst", nfeatures = 1500)

# Scale data for PCA
all.genes <- rownames(EPS_ESC_cell_rpmk)
EPS_ESC_cell_rpmk <- ScaleData(EPS_ESC_cell_rpmk, features = all.genes)

# Perform PCA for dimensionality reduction (default npcs = 50, here set to 40)
EPS_ESC_cell_rpmk <- RunPCA(EPS_ESC_cell_rpmk, features = VariableFeatures(object = EPS_ESC_cell_rpmk), npcs = 40)

# Construct nearest neighbor graph and cluster cells
EPS_ESC_cell_rpmk <- FindNeighbors(EPS_ESC_cell_rpmk, dims = 1:20)
EPS_ESC_cell_rpmk <- FindClusters(EPS_ESC_cell_rpmk, resolution = 1)

# Perform t-SNE and UMAP for visualization
EPS_ESC_cell_rpmk <- RunTSNE(EPS_ESC_cell_rpmk,perplexity =13,dims = 1:5)
EPS_ESC_cell_rpmk <- RunUMAP(EPS_ESC_cell_rpmk,dims = 1:8)

# Find marker genes for each cluster
Idents(EPS_ESC_cell_rpmk) <- EPS_ESC_cell_rpmk@meta.data$cell_type
markers <- FindAllMarkers(EPS_ESC_cell_rpmk, logfc.threshold = 0.1, only.pos = TRUE, min.pct = 0.1)

# Select top 500 marker genes with log2FC > 0.25
top500 <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 500) %>%
  ungroup()

# Heatmap for top marker genes
DoHeatmap(EPS_ESC_cell_rpmk, features = top500$gene, group.by = "cell_type")


### Heterogeneity Analysis using Silhouette Score

# Extract UMAP and t-SNE embeddings
umap1 = EPS_ESC_cell_rpmk@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(tx = EPS_ESC_cell_rpmk@meta.data$cell_type)
tsne1 = EPS_ESC_cell_rpmk@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(tx = EPS_ESC_cell_rpmk@meta.data$cell_type)

# Ensure cell type labels are factors
umap1$tx <- as.factor(umap1$tx)
tsne1$tx <- as.factor(tsne1$tx)

# Compute silhouette scores for UMAP and t-SNE
sil_umap <- silhouette(as.numeric(umap1$tx), dist(umap1[, 1:2]))
sil_tsne <- silhouette(as.numeric(tsne1$tx), dist(tsne1[, 1:2]))

# Convert silhouette scores into data frames
sil_df_umap <- as.data.frame(sil_umap)
sil_df_tsne <- as.data.frame(sil_tsne)

# Calculate average silhouette score per cluster
silhouette_means_umap <- sil_df_umap %>%
  mutate(cluster = as.factor(cluster)) %>%
  group_by(cluster) %>%
  summarize(mean_silhouette = mean(sil_width))

silhouette_means_tsne <- sil_df_tsne %>%
  mutate(cluster = as.factor(cluster)) %>%
  group_by(cluster) %>%
  summarize(mean_silhouette = mean(sil_width))

# Print silhouette scores
print("Silhouette Score for UMAP:")
print(silhouette_means_umap)

print("Silhouette Score for t-SNE:")
print(silhouette_means_tsne)




intersected_genes <- c("LTR78","LTR5","LTR5-Hs","LTR62","LTR14B","LTR12C")

  
  # Print the plot to the PDF
  boxplot_plot <- ggplot(gene_data, aes(x = cell_type, y = get(gene), fill = cell_type)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = gene, y = "Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Print the plot to the PDF
  print(boxplot_plot)
