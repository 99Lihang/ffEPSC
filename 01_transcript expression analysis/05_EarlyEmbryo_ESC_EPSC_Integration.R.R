# Load necessary libraries
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)

# Load datasets
Oocyte_Morulae <- readRDS("Oocyte_Morulae.raw.rds")
# Keep only specific cell types
Oocyte_Morulae <- subset(Oocyte_Morulae, celltype %in% c("Zygote", "2_cell", "4_cell", "8_cell", "hESC", "Morulae", "Late blastocyst"))
# Rename cells to ensure unique cell names when merging
Oocyte_Morulae <- RenameCells(Oocyte_Morulae, add.cell.id = "Oocyte_Morulae")

E3_7 <- readRDS("E3_7.raw.rds")
E3_7 <- RenameCells(E3_7, add.cell.id = "E3_7")

D6_12 <- readRDS("D6_12.raw.rds")
D6_12 <- RenameCells(D6_12, add.cell.id = "D6_12")

TEMP.OBJ <- readRDS("EPS+ESC.rds")

# Merge all datasets into a single Seurat object
obj <- merge(Oocyte_Morulae, y = list(E3_7, D6_12, TEMP.OBJ))

# Correct metadata: update sample identity for EPS and ESC
obj@meta.data[obj@meta.data$celltype == "ffEPS", "orig.ident"] <- "own.EPS"
obj@meta.data[obj@meta.data$celltype == "ESC", "orig.ident"] <- "own.ESC"

# Calculate mitochondrial gene percentage
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Plot QC metrics
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")

# Filter out low-quality cells
obj <- subset(obj, subset = nFeature_RNA > 200)

# Normalize and scale data
merged_data <- obj
merged_data <- NormalizeData(merged_data)
merged_data <- FindVariableFeatures(merged_data)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)

# Perform batch correction using RPCA integration
merged_data <- FindIntegrationAnchors(object.list = list(merged_data), reduction = "rpca", dims = 1:30)
merged_data <- IntegrateData(anchorset = merged_data, dims = 1:30)

# Compute neighbors and clusters
merged_data <- FindNeighbors(merged_data, dims = 1:10)
merged_data <- FindClusters(merged_data, resolution = 2)

# Run UMAP for visualization
merged_data <- RunUMAP(merged_data, dims = 1:10, reduction.name = "umap.rpca")

# Standardize final cell type annotations
merged_data@meta.data$finalcelltype <- NA
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "8_cell"] <- "E3"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "Morulae"] <- "E4"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "D6"] <- "E5"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "D7"] <- "E6"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "Late blastocyst"] <- "E5"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "D8"] <- "E8"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "D9"] <- "E9"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "D10"] <- "E10"
merged_data@meta.data$finalcelltype[merged_data@meta.data$celltype1 == "D12"] <- "E12"

# Define factor levels for ordered cell type visualization
merged_data@meta.data$finalcelltype <- factor(merged_data@meta.data$finalcelltype,
                                              levels = c("Zygote", "2_cell", "4_cell", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E12",
                                                         "EPS.1", "EPS.2", "ESC.1", "ESC.2", "hESC"))

# Define colors for UMAP visualization
gradient_cell_types <- c("Zygote", "2_cell", "4_cell", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E12")
start_color <- "sienna1" # Start color (orange)
mid_color <- "darkolivegreen1" # Middle color (green)
mid_color2 <- "lightseagreen" # Second middle color
end_color <- "purple4" # End color (purple)

# Generate gradient colors
gradient_colors <- colorRampPalette(c(start_color, mid_color, mid_color2, end_color))(length(gradient_cell_types))

# Define specific colors for EPS and ESC populations
specific_colors <- c(
  "EPS.1" = "palevioletred4",
  "EPS.2" = "palevioletred1",
  "ESC.1" = "#FF9933",
  "ESC.2" = "#FFFF00",
  "hESC" = "lightslateblue"
)

# Merge color mappings
cell_type_colors <- setNames(c(gradient_colors, specific_colors), c(gradient_cell_types, names(specific_colors)))

# Check if all cell types have assigned colors
missing_types <- setdiff(unique(merged_data@meta.data$finalcelltype), names(cell_type_colors))
if (length(missing_types) > 0) {
  stop(paste("Missing colors for:", paste(missing_types, collapse = ", ")))
}

# Add color annotations to metadata
merged_data@meta.data$CellTypeColor <- cell_type_colors[merged_data@meta.data$finalcelltype]

# UMAP visualization
DimPlot(merged_data, reduction = "umap.rpca", group.by = "finalcelltype", cols = cell_type_colors, combine = FALSE, label = TRUE)
