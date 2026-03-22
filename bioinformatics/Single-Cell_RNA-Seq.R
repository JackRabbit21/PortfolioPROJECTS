library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(fields)         
library(ggplot2)

set.seed(42)
# Set output directory
output_dir <- "results/"
dir.create(output_dir)

BMMC_D1T1_data <- readRDS("D:\\unecessary stuff\\Desktop\\BIOINFO WS24\\SINGLE CELL\\PROJECT\\scbi_ds1\\GSM4138872_scRNA_BMMC_D1T1.rds")
CD34_D2T1_data <- readRDS("D:\\unecessary stuff\\Desktop\\BIOINFO WS24\\SINGLE CELL\\PROJECT\\scbi_ds1\\GSM4138874_scRNA_CD34_D2T1.rds")
BMMC_D1T2_data <-readRDS("D:\\unecessary stuff\\Desktop\\BIOINFO WS24\\SINGLE CELL\\PROJECT\\scbi_ds1\\GSM4138873_scRNA_BMMC_D1T2.rds")
CD34_D3T1_data <- readRDS("D:\\unecessary stuff\\Desktop\\BIOINFO WS24\\SINGLE CELL\\PROJECT\\scbi_ds1\\GSM4138875_scRNA_CD34_D3T1.rds")

# Create Seurat objects
BMMC_1_seurat <- CreateSeuratObject(counts = BMMC_D1T1_data)
CD34_1_seurat <- CreateSeuratObject(counts = CD34_D2T1_data)
BMMC_2_seurat <- CreateSeuratObject(counts = BMMC_D1T2_data)
CD34_2_seurat <- CreateSeuratObject(counts = CD34_D3T1_data)


#Adding metadata for each sample
BMMC_1_seurat <- AddMetaData(object = BMMC_1_seurat, metadata = "D1",  col.name = 'Donor')
BMMC_1_seurat <- AddMetaData(object = BMMC_1_seurat, metadata = "T1", col.name = 'Replicate')
BMMC_1_seurat <- AddMetaData(object = BMMC_1_seurat, metadata =  "F", col.name = 'Sex')

BMMC_2_seurat <- AddMetaData(object = BMMC_2_seurat, metadata = "D1", col.name = 'Donor')
BMMC_2_seurat <- AddMetaData(object = BMMC_2_seurat, metadata = "T2", col.name = 'Replicate')
BMMC_2_seurat <- AddMetaData(object = BMMC_2_seurat, metadata = "F", col.name = 'Sex')

CD34_1_seurat <- AddMetaData(object = CD34_1_seurat, metadata = "D2", col.name = 'Donor')
CD34_1_seurat <- AddMetaData(object = CD34_1_seurat, metadata = "T1", col.name = 'Replicate')
CD34_1_seurat <- AddMetaData(object = CD34_1_seurat, metadata = "M", col.name = 'Sex')

CD34_2_seurat <- AddMetaData(object = CD34_2_seurat, metadata = "D3", col.name = 'Donor')
CD34_2_seurat <- AddMetaData(object = CD34_2_seurat, metadata = "T1", col.name = 'Replicate')
CD34_2_seurat <- AddMetaData(object = CD34_2_seurat, metadata = "F",col.name = 'Sex')

# How many cells are in each sample?
ncol(BMMC_1_seurat)
ncol(BMMC_2_seurat)
ncol(CD34_1_seurat)
ncol(CD34_2_seurat)

# How many genes are part of the expression-matrix?
nrow(BMMC_1_seurat)
nrow(BMMC_2_seurat)
nrow(CD34_1_seurat)
nrow(CD34_2_seurat)

names(BMMC_1_seurat@meta.data)

# Loading required libraries
library(Seurat)
library(ggplot2)

# Step 1: Add percent.mt metadata
BMMC_1_seurat[["percent.mt"]] <- PercentageFeatureSet(BMMC_1_seurat, pattern = "^MT-")
CD34_1_seurat[["percent.mt"]] <- PercentageFeatureSet(CD34_1_seurat, pattern = "^MT-")
BMMC_2_seurat[["percent.mt"]] <- PercentageFeatureSet(BMMC_2_seurat, pattern = "^MT-")
CD34_2_seurat[["percent.mt"]] <- PercentageFeatureSet(CD34_2_seurat, pattern = "^MT-")

# Step 2: Visualization and Filtering for Each Dataset

# For BMMC_1
vln_BMMC_1 <- VlnPlot(BMMC_1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("VlnPlot_BMMC_1.png", plot = vln_BMMC_1, width = 10, height = 5, units = "in", dpi = 300)

scatter1_BMMC_1 <- FeatureScatter(BMMC_1_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2_BMMC_1 <- FeatureScatter(BMMC_1_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_BMMC_1 <- CombinePlots(plots = list(scatter1_BMMC_1, scatter2_BMMC_1))
ggsave("FeatureScatter_BMMC_1.png", plot = combined_BMMC_1, width = 10, height = 5, units = "in", dpi = 300)

BMMC_1_seurat <- subset(BMMC_1_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# For CD34_1
vln_CD34_1 <- VlnPlot(CD34_1_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("VlnPlot_CD34_1.png", plot = vln_CD34_1, width = 10, height = 5, units = "in", dpi = 300)

scatter1_CD34_1 <- FeatureScatter(CD34_1_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2_CD34_1 <- FeatureScatter(CD34_1_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_CD34_1 <- CombinePlots(plots = list(scatter1_CD34_1, scatter2_CD34_1))
ggsave("FeatureScatter_CD34_1.png", plot = combined_CD34_1, width = 10, height = 5, units = "in", dpi = 300)

CD34_1_seurat <- subset(CD34_1_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# For BMMC_2
vln_BMMC_2 <- VlnPlot(BMMC_2_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("VlnPlot_BMMC_2.png", plot = vln_BMMC_2, width = 10, height = 5, units = "in", dpi = 300)

scatter1_BMMC_2 <- FeatureScatter(BMMC_2_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2_BMMC_2 <- FeatureScatter(BMMC_2_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_BMMC_2 <- CombinePlots(plots = list(scatter1_BMMC_2, scatter2_BMMC_2))
ggsave("FeatureScatter_BMMC_2.png", plot = combined_BMMC_2, width = 10, height = 5, units = "in", dpi = 300)

BMMC_2_seurat <- subset(BMMC_2_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# For CD34_2
vln_CD34_2 <- VlnPlot(CD34_2_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("VlnPlot_CD34_2.png", plot = vln_CD34_2, width = 10, height = 5, units = "in", dpi = 300)

scatter1_CD34_2 <- FeatureScatter(CD34_2_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2_CD34_2 <- FeatureScatter(CD34_2_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_CD34_2 <- CombinePlots(plots = list(scatter1_CD34_2, scatter2_CD34_2))
ggsave("FeatureScatter_CD34_2.png", plot = combined_CD34_2, width = 10, height = 5, units = "in", dpi = 300)

CD34_2_seurat <- subset(CD34_2_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#######################################################################################################################################


#DOUBLETFILTER PROCESS
# Data Normalization for each dataset (Use LogNormalize with a reduced scale factor to save memory)
BMMC_1_seurat <- NormalizeData(BMMC_1_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
CD34_1_seurat <- NormalizeData(CD34_1_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
BMMC_2_seurat <- NormalizeData(BMMC_2_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
CD34_2_seurat <- NormalizeData(CD34_2_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Preprocess 
BMMC_1_seurat <- FindVariableFeatures(BMMC_1_seurat, selection.method = "vst", nfeatures = 1000)
BMMC_1_seurat <- ScaleData(BMMC_1_seurat)
BMMC_1_seurat <- RunPCA(BMMC_1_seurat)

CD34_1_seurat <- FindVariableFeatures(CD34_1_seurat, selection.method = "vst", nfeatures = 1000)
CD34_1_seurat <- ScaleData(CD34_1_seurat)
CD34_1_seurat <- RunPCA(CD34_1_seurat)

BMMC_2_seurat <- FindVariableFeatures(BMMC_2_seurat, selection.method = "vst", nfeatures = 1000)
BMMC_2_seurat <- ScaleData(BMMC_2_seurat)
BMMC_2_seurat <- RunPCA(BMMC_2_seurat)

CD34_2_seurat <- FindVariableFeatures(CD34_2_seurat, selection.method = "vst", nfeatures = 1000)
CD34_2_seurat <- ScaleData(CD34_2_seurat)
CD34_2_seurat <- RunPCA(CD34_2_seurat)


# Create the ElbowPlot
elbow_plot <- ElbowPlot(BMMC_1_seurat, ndims = 50)
ggsave("ElbowPlot_BMMC_1.png", plot = elbow_plot, width = 10, height = 8, units = "in", dpi = 300)

elbow_plot <- ElbowPlot(CD34_1_seurat, ndims = 50)
ggsave("ElbowPlot_CD34_1_seurat.png", plot = elbow_plot, width = 10, height = 8, units = "in", dpi = 300)

elbow_plot <- ElbowPlot(BMMC_2_seurat, ndims = 50)
ggsave("ElbowPlot_BMMC_2_seurat.png", plot = elbow_plot, width = 10, height = 8, units = "in", dpi = 300)

elbow_plot <- ElbowPlot(CD34_2_seurat, ndims = 50)
ggsave("ElbowPlot_CD34_2_seurat.png", plot = elbow_plot, width = 10, height = 8, units = "in", dpi = 300)

#UMAP
BMMC_1_seurat <- RunUMAP(BMMC_1_seurat, dims = 1:10)
DimPlot(BMMC_1_seurat, reduction = "umap")
CD34_1_seurat <- RunUMAP(CD34_1_seurat, dims = 1:10)
DimPlot(CD34_1_seurat, reduction = "umap")
BMMC_2_seurat <- RunUMAP(BMMC_2_seurat, dims = 1:10)
DimPlot(BMMC_2_seurat, reduction = "umap")
CD34_2_seurat <- RunUMAP(CD34_2_seurat, dims = 1:10)
DimPlot(CD34_2_seurat, reduction = "umap")
#pK Identification (ground-truth)

sweep.res.BMMC_1_seurat <- paramSweep(BMMC_1_seurat, PCs = 1:10, sct = FALSE)
gt.calls <- BMMC_1_seurat@meta.data[rownames(sweep.res.BMMC_1_seurat[[1]]), "GT"]  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_BMMC_1_seurat <- summarizeSweep(sweep.res.BMMC_1_seurat, GT = FALSE)
bcmvn_BMMC_1_seurat <- find.pK(sweep.stats_BMMC_1_seurat)

sweep.res.CD34_1_seurat <- paramSweep(CD34_1_seurat, PCs = 1:10, sct = FALSE)
gt.calls <- CD34_1_seurat@meta.data[rownames(sweep.res.CD34_1_seurat[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_CD34_1_seurat <- summarizeSweep(sweep.res.CD34_1_seurat, GT = FALSE)
bcmvn_CD34_1_seurat <- find.pK(sweep.stats_CD34_1_seurat)

sweep.res.BMMC_2_seurat <- paramSweep(BMMC_2_seurat, PCs = 1:10, sct = FALSE)
gt.calls <- BMMC_2_seurat@meta.data[rownames(sweep.res.BMMC_2_seurat[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_BMMC_2_seurat<- summarizeSweep(sweep.res.BMMC_2_seurat, GT = FALSE)
bcmvn_BMMC_2_seurat <- find.pK(sweep.stats_BMMC_2_seurat)

sweep.res.CD34_2_seurat <- paramSweep(CD34_2_seurat, PCs = 1:10, sct = FALSE)
gt.calls <- CD34_2_seurat@meta.data[rownames(sweep.res.CD34_2_seurat[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_CD34_2_seurat<- summarizeSweep(sweep.res.CD34_2_seurat, GT = FALSE)
bcmvn_CD34_2_seurat <- find.pK(sweep.stats_CD34_2_seurat)

doublet_rate <- 0.05
nExp_poi <- round(ncol(BMMC_1_seurat) * doublet_rate)
nExp_poi <- round(ncol(CD34_1_seurat) * doublet_rate)
nExp_poi <- round(ncol(BMMC_2_seurat) * doublet_rate)
nExp_poi <- round(ncol(CD34_2_seurat) * doublet_rate)

## Run DoubletFinder with varying classification stringencies
BMMC_1_seurat <- doubletFinder(BMMC_1_seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
CD34_1_seurat <- doubletFinder(CD34_1_seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
BMMC_2_seurat <- doubletFinder(BMMC_2_seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
CD34_2_seurat <- doubletFinder(CD34_2_seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)



#Merging Dataset (without batch correction)
data_all <- merge(
  BMMC_1_seurat,
  y = c(BMMC_2_seurat, CD34_1_seurat, CD34_2_seurat),
  add.cell.ids = list("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"),
  project = "Project1"
)

#Normalization
data_all <- NormalizeData(data_all, normalization.method = "LogNormalize", scale.factor = 10000)
data_all <- NormalizeData(data_all)

#Feature Selection
data_all <- FindVariableFeatures(data_all, selection.method = "vst", nfeatures = 1000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data_all), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

data_all <- RunPCA(data_all, npcs = 30)
data_all <- RunUMAP(data_all, dims = 1:10)
DimPlot(data_all, reduction = "umap", group.by = "orig.ident")

#################################################################################################################
#Batch Correction (changed to avoid memory overload)

library(Seurat)
library(dplyr)

# Step 1: Assign batch labels to each dataset
BMMC_1_seurat$batch <- "BMMC_1"
CD34_1_seurat$batch <- "CD34_1"
BMMC_2_seurat$batch <- "BMMC_2"
CD34_2_seurat$batch <- "CD34_2"

# Combine all datasets into a list
seurat_list <- list(BMMC_1_seurat, CD34_1_seurat, BMMC_2_seurat, CD34_2_seurat)

# Step 2: Downsample cells to reduce memory usage
seurat_list <- lapply(seurat_list, function(x) {
  x <- subset(x, cells = sample(Cells(x), min(1000, ncol(x)))) 
  return(x)
})

# Step 3: Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 1000) 

# Step 4: Process in chunks to avoid overloading memory
# Integrate the first two datasets
anchors_1 <- FindIntegrationAnchors(object.list = seurat_list[1:2], anchor.features = features, dims = 1:10) 
combined_seurat_1 <- IntegrateData(anchorset = anchors_1, dims = 1:10)

# Integrate the next dataset with the integrated object
anchors_2 <- FindIntegrationAnchors(object.list = list(combined_seurat_1, seurat_list[[3]]), anchor.features = features, dims = 1:10)
combined_seurat_2 <- IntegrateData(anchorset = anchors_2, dims = 1:10)

# Integrate the final dataset
anchors_3 <- FindIntegrationAnchors(object.list = list(combined_seurat_2, seurat_list[[4]]), anchor.features = features, dims = 1:10)
combined_seurat <- IntegrateData(anchorset = anchors_3, dims = 1:10)

# Switch the default assay to the integrated assay
DefaultAssay(combined_seurat) <- "integrated"

# Normalize, scale, and run PCA on the integrated dataset
combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(combined_seurat, npcs = 15, verbose = FALSE)  

# Find neighbors and clusters
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:15)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Run UMAP
combined_seurat <- RunUMAP(combined_seurat, dims = 1:15)

# Plot UMAP by batch
umap_by_batch <- DimPlot(combined_seurat, reduction = "umap", group.by = "batch", label = TRUE) +
  ggtitle("UMAP by Batch")
ggsave("UMAP_by_Batch.png", plot = umap_by_batch, width = 8, height = 6, dpi = 300)

# Plot UMAP by clusters
umap_by_clusters <- DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP by Clusters")
ggsave("UMAP_by_Clusters.png", plot = umap_by_clusters, width = 8, height = 6, dpi = 300)


##########################################################################################################################33
head(combined_seurat)
###WEEK 3
#Cell type annotation

##AUTOMATIC ANNOTATION
# Load necessary libraries
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

annotate <- function(data, name) {
  # Load reference data
  ref_data <- HumanPrimaryCellAtlasData()
  
  # Run SingleR
  sce <- as.SingleCellExperiment(data)
  predictions <- SingleR(test = sce,
                         ref = ref_data,
                         labels = ref_data$label.main)
  
  # Add predictions to Seurat object
  data$cell_type <- predictions$labels
  
  # Create and save UMAP plot
  p <- DimPlot(data, 
               reduction = "umap", 
               group.by = "cell_type", 
               label = TRUE) +
    ggtitle(name) +
    theme(legend.position = "right")
  
  # Save individual plot
  ggsave(filename = file.path(output_dir, paste0(name, "_singleR_annotation.png")),
         plot = p,
         width = 12,
         height = 8)
  
  return(list(object = data, plot = p))
}
results_bmmc_d1t1 <- annotate(BMMC_1_seurat, "BMMC_1")
results_bmmc_d1t2 <- annotate(BMMC_2_seurat, "BMMC_2")
results_cd34_d2t1 <- annotate(CD34_1_seurat, "CD34_1")
results_cd34_d3t1 <- annotate(CD34_2_seurat, "CD34_2")


# Create and save combined plot
combined_plot <- (results_bmmc_d1t1$plot + results_bmmc_d1t2$plot) /
  (results_cd34_d2t1$plot + results_cd34_d3t1$plot)
ggsave(filename = file.path(output_dir, "combined_singleR_annotation.png"),
       plot = combined_plot,
       width = 20,
       height = 16)


###################################################################################################
#Manual Annotation

perform_manual_annotation <- function(seurat_obj, name) {
  # Find clusters
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  # Find markers for each cluster
  all_markers <- FindAllMarkers(seurat_obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
  
  # Save markers
  write.csv(all_markers, 
            file = file.path(output_dir, paste0(name, "_Differential_expression.csv")))
  
  # Manual annotation based on Table 2 markers
  
  cell_type_markers <- list(
    "HSC" = c("CD34", "CD38", "Scal", "Kit"),
    "LMPP" = c("CD38", "CD52", "CSF3R", "c-fms", "Kit", "CD34", "CD38"),
    "CLP" = c("FLT3"),
    "GMP_Neutrophiles" = c("ELANE"),
    "Common_Myeloid_Progenitor" = c("IL3", "GM-CSF", "M-CSF"),
    "B_cells" = c("CD19", "CD20", "CD38"),
    "Pre_B_cells" = c("CD19", "CD34"),
    "Plasma" = c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN"),
    "T_cells_CD8" = c("CD3D", "CD3E", "CD8A", "CD8B"),
    "T_cells_CD4" = c("CD3D", "CD3E", "CD4"),
    "NK_cells" = c("FCGR3A", "NCAM1", "NKG7", "KLRB1"),
    "Erythrocytes" = c("GATA1", "HBB", "HBA1", "HBA2"),
    "PDC" = c("IRF8", "IRF4", "IRF7"),
    "EDC" = c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA"),
    "Monocytes_CD14" = c("CD14", "CCL2", "CCL4", "IL1B"),
    "Monocytes_CD16" = c("FCGR3A", "CD68", "S100A12"),
    "Basophils" = c("GATA2")
  )
  
  # Function to assign cell types based on marker expression
  assign_cell_types <- function(cluster_markers, cell_type_markers) {
    cluster_assignments <- character(length(unique(cluster_markers$cluster)))
    names(cluster_assignments) <- unique(cluster_markers$cluster)  # Assign cluster names to assignments
    
    for(cluster in unique(cluster_markers$cluster)) {
      cluster_genes <- cluster_markers$gene[cluster_markers$cluster == cluster]
      if (length(cluster_genes) == 0) next  # Skip if no genes in the cluster
      
      best_match <- "Unknown"
      max_overlap <- 0
      for(cell_type in names(cell_type_markers)) {
        overlap <- length(intersect(cluster_genes[1:50], cell_type_markers[[cell_type]]))
        if(overlap > max_overlap) {
          max_overlap <- overlap
          best_match <- cell_type
        }
      }
      cluster_assignments[as.character(cluster)] <- best_match  # Ensure assignment by cluster name
    }
    
    # Check if any assignments are still NA and set them to "Unknown"
    cluster_assignments[is.na(cluster_assignments)] <- "Unknown"
    
    return(cluster_assignments)
  }
  
  # Assign cell types
  cluster_assignments <- assign_cell_types(all_markers, cell_type_markers)
  names(cluster_assignments) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, cluster_assignments)
  
  # Create UMAP plots
  manual_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle(paste0(name, " - Manual Annotation"))
  
  # Plot marker genes (select 3 key markers)
  key_markers <- c("CD34", "CD19", "CD3E")  # Example markers
  marker_plots <- lapply(key_markers, function(marker) {
    FeaturePlot(seurat_obj, features = marker) +
      ggtitle(paste0(name, " - ", marker))
  })
  
  # Violin plots for the same markers
  violin_plots <- VlnPlot(seurat_obj, features = key_markers, ncol = 3)
  
  # Save plots
  ggsave(filename = file.path(output_dir, paste0(name, "_manual_annotation_umap.png")),
         plot = manual_umap,
         width = 12,
         height = 8)
  
  lapply(1:length(marker_plots), function(i) {
    ggsave(filename = file.path(output_dir, paste0(name, "_", key_markers[i], "_feature_plot.png")),
           plot = marker_plots[[i]],
           width = 12,
           height = 8)
  })
  
  ggsave(filename = file.path(output_dir, paste0(name, "_violin_plots.png")),
         plot = violin_plots,
         width = 12,
         height = 8)
  
  return(seurat_obj)
}

# Apply manual annotation to all Seurat objects
seurat_bmmc_d1t1 <- perform_manual_annotation(results_bmmc_d1t1$object, "BMMC_D1T1")
seurat_bmmc_d1t2 <- perform_manual_annotation(results_bmmc_d1t2$object, "BMMC_D1T2")
seurat_cd34_d2t1 <- perform_manual_annotation(results_cd34_d2t1$object, "CD34_D2T1")
seurat_cd34_d3t1 <- perform_manual_annotation(results_cd34_d3t1$object, "CD34_D3T1")