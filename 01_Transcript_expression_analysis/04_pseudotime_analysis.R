############### Load Required Libraries ###############
library(monocle)

# Retrieve the raw counts matrix using the layer argument
raw_counts_matrix <- GetAssayData(object = seurat_obj, layer = "counts", assay = "RNA")

# Convert the data to a sparse matrix
sparse_raw_counts_matrix <- as(as.matrix(raw_counts_matrix), 'sparseMatrix')

# Extract metadata for cells
pdata <- seurat_obj@meta.data

# Add cluster annotation information
pdata$group <- seurat_obj@meta.data$finalcelltype

# Extract gene name information
fdata <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))

# Construct AnnotatedDataFrame objects
pd <- new('AnnotatedDataFrame', data = pdata)
fd <- new('AnnotatedDataFrame', data = fdata)

# Create CellDataSet object
cds <- newCellDataSet(sparse_raw_counts_matrix, 
                      phenoData = pd, 
                      featureData = fd, 
                      expressionFamily = VGAM::negbinomial.size())

############### Normalization and Dispersion Estimation ###############
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Define ordering genes (Ensure express_genes is defined beforehand)
cds <- setOrderingFilter(cds, express_genes$gene)
plot_ordering_genes(cds)

# Generate dispersion table to select ordering genes for trajectory analysis
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.5 & 
                           dispersion_empirical >= 1.5 * dispersion_fit)$gene_id

# Apply ordering genes filter
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

############### Dimensionality Reduction and Trajectory Inference ###############
cds <- reduceDimension(cds, max_components = 3, method = 'DDRTree')

# Plot cell trajectory colored by different features
plot_cell_trajectory(cds, color_by = "subcelltype")
plot_cell_trajectory(cds, color_by = "Pseudotime")

############### Differential Gene Expression Analysis Along Pseudotime ###############
diffgene <- differentialGeneTest(cds[ordering_genes, ], 
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)", 
                                 cores = 4)  # Time-consuming step

# Save differential gene results
write.csv(diffgene, file = "pseudotime_degs.csv")

############### Pseudotime Expression Heatmap and Individual Gene Plots ###############
# Set directory paths for gene sets and output plots
geneset_dir <- "/geneset/"
output_dir <- "/plot/"

# List all gene set files (.txt) in the directory
gene_files <- list.files(geneset_dir, pattern = "\\.txt$", full.names = TRUE)

# Define custom colors for different groups
group_colors <- c("ffEPSC_C1" = "#953751", "ffEPSC_C2" = "#F1A69E", 
                  "ESC_C1" = "#006995", "ESC_C2" = "#B7CBC2")

# Loop through each gene file and generate plots
for (gene_file in gene_files) {
  # Extract the file name without extension
  gene_file_name <- tools::file_path_sans_ext(basename(gene_file))
  
  # Read gene list from file
  gene_data <- read.table(gene_file)
  gene_data$V1 <- gsub("−", "-", gene_data$V1)  # Replace special characters
  
  # Select genes that are present in the dataset
  genes <- row.names(subset(fData(cds), gene_short_name %in% gene_data$V1))
  
  # Generate and save pseudotime heatmap
  heatmap_file <- paste0(output_dir, "Pseudotemporal_Expression_", gene_file_name, ".pdf")
  pdf(heatmap_file)
  plot_pseudotime_heatmap(cds[genes, ], cores = 1, cluster_rows = TRUE, show_rownames = TRUE)
  dev.off()
  
  # Generate individual gene expression plots over pseudotime
  pseudotime_plot_file <- paste0(output_dir, "traject_", gene_file_name, "_Pseudotime.pdf")
  pdf(pseudotime_plot_file, height = 1.5, width = 4)
  for (gene in genes) {
    cds_subset <- cds[gene, ]  # Subset the gene
    p <- plot_genes_in_pseudotime(cds_subset, color_by = "finalcelltype") +
      scale_color_manual(values = group_colors)  # Apply custom colors
    print(p)
  }
  dev.off()
}