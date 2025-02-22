# Load required libraries
library(Seurat)
library(dplyr)
library(plotrix)
library(parallel)
library(gplots)

# Read counts matrix
counts <- read.table("/home/Gh38/counts/counts.txt")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts,
                                 min.features = 200,
                                 min.cells = 3, 
                                 project = "owndata")

### ================================

# Normalize data using log-normalization
EPS_ESC_cell <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
EPS_ESC_cell <- FindVariableFeatures(EPS_ESC_cell, selection.method = "vst", nfeatures = 4500)

# Scale data for PCA
all.genes <- rownames(EPS_ESC_cell)
EPS_ESC_cell <- ScaleData(EPS_ESC_cell, features = all.genes)

# Perform PCA for dimensionality reduction (default npcs = 50, here set to 40)
EPS_ESC_cell <- RunPCA(EPS_ESC_cell, features = VariableFeatures(object = EPS_ESC_cell), npcs = 40)

# Construct nearest neighbor graph and cluster cells
EPS_ESC_cell <- FindNeighbors(EPS_ESC_cell, dims = 1:20)
EPS_ESC_cell <- FindClusters(EPS_ESC_cell, resolution = 1.3)


# Perform t-SNE and UMAP for visualization
EPS_ESC_cell <- RunTSNE(EPS_ESC_cell, perplexity = 13)
EPS_ESC_cell <- RunUMAP(EPS_ESC_cell, dims = 1:15)

# Find marker genes for each cluster
Idents(EPS_ESC_cell) <- EPS_ESC_cell@meta.data$cell_type
markers <- FindAllMarkers(EPS_ESC_cell, logfc.threshold = 0.1, only.pos = TRUE, min.pct = 0.1)

# Select top 500 marker genes with log2FC > 0.25
top500 <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 500) %>%
  ungroup()

# Heatmap for top marker genes
DoHeatmap(EPS_ESC_cell, features = top500$gene, group.by = "cell_type")

### Heterogeneity Analysis using Silhouette Score

# Extract UMAP and t-SNE embeddings
umap1 = EPS_ESC_cell@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(tx = EPS_ESC_cell@meta.data$cell_type)
tsne1 = EPS_ESC_cell@reductions$tsne@cell.embeddings %>% as.data.frame() %>% cbind(tx = EPS_ESC_cell@meta.data$cell_type)

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

### Correlation Analysis of Gene Expression

# Extract expression matrix and metadata
a <- as.matrix(EPS_ESC_cell@assays$RNA@data)
a.meta <- data.table(EPS_ESC_cell@meta.data, keep.rownames = "cell")
a.meta <- a.meta[, c("cell", "cell_name")]

# Filter low-expression genes
a <- a[rowMeans(a) > 0.1, ]

# Compute mean expression for each cell type
meantable <- mclapply(unique(a.meta$cell_name), function(x) {
  meta1 <- a.meta[cell_name == x, cell]
  rowMeans(subset(a, select = meta1))
}, mc.cores = 1)

mergemean <- do.call(cbind, meantable)
colnames(mergemean) <- unique(a.meta$cell_name)

# Filter genes with mean expression > 0.1
mergemean <- mergemean[rowMeans(mergemean) > 0.1, ]

# Compute Spearman correlation matrix
spearcol1 <- cor(mergemean, method = "spearman")


heatmap.2(spearcol1,
          dendrogram = "both",  
          trace = "none",
          density.info = "none",
          Colv = TRUE,  
          Rowv = TRUE, 
          #cellnote = round(spearcol1, 4),
          #notecex = 0.5,
          #notecol = "black",
          col = colorRampPalette(c("#D83830","#FAF6DA", "#4471A7"))(100),
          keysize = 0.5)  


# Read cell marker data
cellmarkers <- read.table("ZGA_totiptent_naive_primed_markers.txt")

# Extract gene names from the dataset
gene_names <- rownames(EPS_ESC_cell)

# Select genes from cellmarkers that are present in the EPS_ESC_cell dataset
cellmarkers <- cellmarkers[cellmarkers$V1 %in% gene_names, , drop = FALSE]

# Number of genes per plot
genes_per_plot <- 8

# Loop through and generate violin plots
for (i in 1:num_plots) {
  # Define the range of genes to be plotted in each iteration
  start_gene <- (i - 1) * genes_per_plot + 1
  end_gene <- min(i * genes_per_plot, nrow(cellmarkers))
  
  # Extract the specific genes to be plotted
  genes_to_plot <- cellmarkers$V1[start_gene:end_gene]
  
  # Generate violin plots grouped by cell type
  vln_plot <- VlnPlot(EPS_ESC_cell, features = genes_to_plot, group.by = "cell_type", ncol = 4,
                      cols = custom_colors)
  
  # Output the plot to a PDF document
  print(vln_plot)
}
# Close the PDF device (if applicable)




