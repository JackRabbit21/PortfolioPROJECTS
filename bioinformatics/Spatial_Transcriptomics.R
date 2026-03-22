library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
set.seed(42)

# Define input directories for the datasets
data_dir_sample1 <- "/home/kech00002/projects/scbi_p3/project_3_dataset/Section_1"
data_dir_sample2 <- "/home/kech00002/projects/scbi_p3/project_3_dataset/Section_2"

#output directory 
output_dir <- "/home/kech00002/projects/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load spatial data
sample1 <- Load10X_Spatial(
  data.dir = data_dir_sample1,
  filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5"
)

sample2 <- Load10X_Spatial(
  data.dir = data_dir_sample2,
  filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5"
)

#WEEK1
# Q1.1: Properties of the Slides
num_spots_sample1 <- length(colnames(sample1))
num_spots_sample2 <- length(colnames(sample2))
spot_spacing <- sample1@images$slice1@scale.factors$spot
cat("Number of spots in Sample 1:", num_spots_sample1, "\n")
cat("Number of spots in Sample 2:", num_spots_sample2, "\n")
cat("Spot spacing:", spot_spacing, "\n")

# Q1.2: Resolution of the spatial transcriptomics technology
avg_cell_size <- 10 
resolution_comparison <- spot_spacing / avg_cell_size
cat("Resolution comparison (spot spacing to cell size):", resolution_comparison, "\n")

# Q1.3: Output of Space Ranger
counts_matrix_sample1 <- GetAssayData(sample1, assay = "Spatial", slot = "counts")
counts_matrix_sample2 <- GetAssayData(sample2, assay = "Spatial", slot = "counts")
cat("Dimensions of the gene-expression matrix (Sample 1):", dim(counts_matrix_sample1), "\n")
cat("Dimensions of the gene-expression matrix (Sample 2):", dim(counts_matrix_sample2), "\n")

# Save tissue Sample 1
jpeg(file.path(output_dir, "tissue_image_sample1.jpg"))
SpatialFeaturePlot(sample1, features = "nCount_Spatial")
dev.off()

# Save tissue Sample 2
jpeg(file.path(output_dir, "tissue_image_sample2.jpg"))
SpatialFeaturePlot(sample2, features = "nCount_Spatial")
dev.off()

sample1 <- SCTransform(sample1, assay = "Spatial", verbose = FALSE)
sample2 <- SCTransform(sample2, assay = "Spatial", verbose = FALSE)

# 3.1
SpatialFeaturePlot(sample1, features = c("Hpca", "Ttr"))
SpatialFeaturePlot(sample2, features = c("Hpca", "Ttr"))

VlnPlot(sample1, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2)
sample1 <- subset(sample1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 2500 & nCount_Spatial < 10000)

VlnPlot(sample2, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2)
sample2 <- subset(sample2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 2500 & nCount_Spatial < 10000)
#3.2
sample1 <- SCTransform(sample1, assay = "Spatial", verbose = FALSE)
sample2 <- SCTransform(sample2, assay = "Spatial", verbose = FALSE)

# 4.1& 4.2 (SAMPLE1)
sample1 <- RunPCA(sample1, assay = "SCT", verbose = FALSE)
sample1 <- FindNeighbors(sample1, reduction = "pca", dims = 1:30)
sample1 <- FindClusters(sample1, verbose = FALSE)
sample1 <- RunUMAP(sample1, reduction = "pca", dims = 1:30)
ElbowPlot(sample1, ndims = 50)
# 4.1& 4.2 Sample 2
sample2 <- RunPCA(sample2, assay = "SCT", verbose = FALSE)
sample2 <- FindNeighbors(sample2, reduction = "pca", dims = 1:30)
sample2 <- FindClusters(sample2, verbose = FALSE)
sample2 <- RunUMAP(sample2, reduction = "pca", dims = 1:30)
ElbowPlot(sample2, ndims = 50)

#UMAP plots
jpeg(file.path(output_dir, "umap_sample1.jpg"))
DimPlot(sample1, reduction = "umap", label = TRUE)
dev.off()

jpeg(file.path(output_dir, "umap_sample2.jpg"))
DimPlot(sample2, reduction = "umap", label = TRUE)
dev.off()

# Save Seurat objects for further analysis
saveRDS(sample1, file = file.path(output_dir, "sample1_processed.rds"))
saveRDS(sample2, file = file.path(output_dir, "sample2_processed.rds"))
cat("Seurat objects saved in:", output_dir, "\n")

#WEEK2
# Load the Seurat object
rds_file_path <- "/home/kech00002/projects/output/sample2_processed.rds"
loaded_sample2 <- readRDS(rds_file_path)

# Confirm the object is loaded
print(loaded_sample2)

DefaultAssay(loaded_sample2) <- "Spatial"

# Normalize and scale the data 
loaded_sample2 <- NormalizeData(loaded_sample2, verbose = FALSE)
loaded_sample2 <- ScaleData(loaded_sample2, verbose = FALSE)

# 5.1 DEG Analysis Based on Clustering
# Perform DEG analysis to compare gene expression between clusters
deg_results <- FindAllMarkers(loaded_sample2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(deg_results, "deg_clustering_results.csv")  
head(deg_results)

# Identify significant DEGs with avg_log2FC > 1
significant_deg <- deg_results %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1)

print("Significant DEGs (avg_log2FC > 1):")
print(significant_deg)

#5.2
DefaultAssay(loaded_sample2) <- "Spatial"

# Normalize the data 
loaded_sample2 <- NormalizeData(loaded_sample2, verbose = FALSE)

# Find variable features
loaded_sample2 <- FindVariableFeatures(loaded_sample2, 
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)

# Scale the data
loaded_sample2 <- ScaleData(loaded_sample2, 
                           features = VariableFeatures(loaded_sample2),
                           verbose = FALSE)

# Find spatially variable features 
spatial_features <- FindVariableFeatures(loaded_sample2,
                                       selection.method = "vst",
                                       nfeatures = 2000,
                                       verbose = FALSE)

# Get top 3 variable features
top3_spatial_features <- head(VariableFeatures(loaded_sample2), 3)

print("Top 3 variable features:")
print(top3_spatial_features)

#6 Merging Data
# Define paths
output_dir <- "/home/kech00002/projects/output"

# Load spatial datasets
sample1 <- Load10X_Spatial(
    data.dir = "/home/kech00002/projects/scbi_p3/project_3_dataset/Section_1",
    filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5"
)

sample2 <- Load10X_Spatial(
    data.dir = "/home/kech00002/projects/scbi_p3/project_3_dataset/Section_2",
    filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5"
)

# Normalize and scale each dataset 
sample1 <- NormalizeData(sample1, verbose = FALSE)
sample1 <- FindVariableFeatures(sample1, selection.method = "vst", nfeatures = 2000)
sample1 <- ScaleData(sample1, verbose = FALSE)

sample2 <- NormalizeData(sample2, verbose = FALSE)
sample2 <- FindVariableFeatures(sample2, selection.method = "vst", nfeatures = 2000)
sample2 <- ScaleData(sample2, verbose = FALSE)

# Perform PCA 
sample1 <- RunPCA(sample1, verbose = FALSE)
sample2 <- RunPCA(sample2, verbose = FALSE)

# Merging Without Batch Correction
merged_no_correction <- merge(sample1, y = sample2, project = "Merged_No_Correction")
merged_no_correction <- ScaleData(merged_no_correction, verbose = FALSE)
merged_no_correction <- RunPCA(merged_no_correction, verbose = FALSE)
merged_no_correction <- FindNeighbors(merged_no_correction, dims = 1:30)
merged_no_correction <- FindClusters(merged_no_correction, resolution = 0.5)
merged_no_correction <- RunUMAP(merged_no_correction, dims = 1:30)

# Plot UMAP and spatial clustering without batch correction
DimPlot(merged_no_correction, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
    ggtitle("UMAP without Batch Correction")

# Spatial plot showing clustering
SpatialDimPlot(merged_no_correction, images = NULL, group.by = "seurat_clusters")

# Optionally, saving plots to output directory
ggsave(file.path(output_dir, "UMAP_no_batch_correction.png"))
ggsave(file.path(output_dir, "Spatial_clustering_no_batch_correction.png"))

table(merged_no_correction$seurat_clusters, merged_no_correction$orig.ident)

#batch correction
# Normalize and scale datasets 
sample1 <- NormalizeData(sample1, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

sample2 <- NormalizeData(sample2, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# Perform integration (batch correction)
anchors <- FindIntegrationAnchors(object.list = list(sample1, sample2), dims = 1:30)
merged_with_correction <- IntegrateData(anchorset = anchors, dims = 1:30)

# Preprocess the integrated data
merged_with_correction <- ScaleData(merged_with_correction, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:30)

# Plot UMAP with batch correction
DimPlot(merged_with_correction, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP with Batch Correction")
SpatialDimPlot(merged_with_correction, images = NULL, group.by = "seurat_clusters")

#WEEK3

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)

set.seed(42)
# Set options to manage memory usage
options(future.globals.maxSize = 3 * 1024^3)


reference_data_path <- "/home/kech00002/projects/scbi_p3/project_3_dataset/allen_cortex.rds"
sample_file_path <- "/home/kech00002/projects/output/merged_with_correction_seurat_object.rds"

# Load datasets
reference_data <- readRDS(reference_data_path)
spatial_data <- readRDS(sample_file_path)

print("Initial assays:")
print(Assays(spatial_data))


DefaultAssay(spatial_data) <- "Spatial"

# Process reference data
reference_data <- SCTransform(reference_data, verbose = FALSE)
reference_data <- RunPCA(reference_data, npcs = 30, verbose = FALSE)
reference_data <- FindNeighbors(reference_data, dims = 1:30)
reference_data <- FindClusters(reference_data, resolution = 0.5)

reference_data <- RunUMAP(reference_data, dims = 1:30, return.model = TRUE)

# Processing spatial data with standard normalization
spatial_data <- NormalizeData(spatial_data)
spatial_data <- FindVariableFeatures(spatial_data)
spatial_data <- ScaleData(spatial_data)
spatial_data <- RunPCA(spatial_data, npcs = 30)

# Find transfer anchors
anchors <- FindTransferAnchors(
    reference = reference_data,
    query = spatial_data,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:30
)

# Transfer labels
spatial_data <- MapQuery(
    anchorset = anchors,
    query = spatial_data,
    reference = reference_data,
    refdata = list(cell_type = "cluster"),
    reference.reduction = "pca",
    reduction.model = "umap"
)

# Plot UMAP with annotations
DimPlot(spatial_data, reduction = "ref.umap", group.by = "predicted.cell_type") +
    ggtitle("Automatic Annotation in UMAP Space")


#8 Deconvolution

# Define file paths
reference_data_path <- "/home/kech00002/projects/scbi_p3/project_3_dataset/allen_cortex.rds"
sample_file_path <- "/home/kech00002/projects/output/merged_with_correction_seurat_object.rds"

# Load datasets
reference_data <- readRDS(reference_data_path)
spatial_data <- readRDS(sample_file_path)

#8.2
# Create a data frame with cell names
cell_data <- data.frame(
  cell_name = colnames(reference_data),
  cluster = reference_data$cluster
)

# Function to sample cells from each cluster
sample_cluster <- function(cells, n_target = 250) {
  n_available <- length(cells)
  n_sample <- min(n_target, n_available)
  sample(cells, size = n_sample, replace = FALSE)
}

# Perform downsampling
downsampled_cells <- cell_data %>%
  group_by(cluster) %>%
  summarise(
    sampled_cells = list(sample_cluster(cell_name)),
    .groups = "drop"
  ) %>%
  pull(sampled_cells) %>%
  unlist()

# Subset the Seurat object using the correct cell names
reference_data_downsampled <- subset(reference_data, cells = downsampled_cells)

#PART 8.3: DEG Analysis

#list of genes present in both datasets
common_genes <- intersect(rownames(reference_data_downsampled), rownames(spatial_data))

# Subset the reference data to include only common genes
reference_data_downsampled <- subset(reference_data_downsampled, features = common_genes)


# Perform DEG analysis
Idents(reference_data_downsampled) <- "cluster"
deg_results <- FindAllMarkers(
  reference_data_downsampled,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Get the top 20 DEGs with the lowest p-values for each cluster
top_degs <- deg_results %>%
  group_by(cluster) %>%          # Group by cluster
  slice_min(order_by = p_val, n = 20) %>%  # Select top 20 with smallest p_val
  ungroup()                      

# View the top DEGs
print(top_degs)


# Save results
write.csv(deg_results, "all_deg_results.csv", row.names = FALSE)
write.csv(top_degs, "top_20_degs_per_cluster.csv", row.names = FALSE)

# Create a summary using dplyr
deg_summary <- top_degs %>%
  group_by(cluster) %>%
  summarise(n_genes = n()) %>%
  arrange(cluster)


cat("\nNumber of top DEGs selected per cluster:\n")
print(deg_summary)

selected_genes <- unique(top_degs$gene)
saveRDS(selected_genes, "selected_degs_for_deconvolution.rds")


cat("\nTotal number of unique genes selected:", length(selected_genes), "\n")
cat("Number of clusters with DEGs:", nrow(deg_summary), "\n")