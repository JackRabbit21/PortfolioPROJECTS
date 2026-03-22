# Load the ArchR package
library(ArchR)
library(ggplot2)
library(ComplexHeatmap)
# Set a seed for reproducibility
set.seed(42)


# Specify the genome
addArchRGenome("hg38")

addArchRThreads(threads = 16)
###############################################################################################################################
#WEEK1

# Define the input fragment files and sample names
inputFiles <- c(
 "/home/kevin-koshy/Desktop/Academics/Single-Cell/scbi_p2/scbi_p2/hft_ctx_w21_dc1r3_r1_atac_fragments.tsv.gz",
  "/home/kevin-koshy/Desktop/Academics/Single-Cell/scbi_p2/scbi_p2/hft_ctx_w21_dc2r2_r1_atac_fragments.tsv.gz",
  "/home/kevin-koshy/Desktop/Academics/Single-Cell/scbi_p2/scbi_p2/hft_ctx_w21_dc2r2_r2_atac_fragments.tsv.gz"
)

sampleNames <- c("sample1", "sample2", "sample3")



#Arrow Creation
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,              
  sampleNames = sampleNames,            
  addTileMat = TRUE,                    
  addGeneScoreMat = TRUE,               
  minFrags = 500,                    
  minTSS = 4,                        
  
)
# Print the created Arrow files
ArrowFiles

#Create arch project
proj <- ArchRProject(
 ArrowFiles = ArrowFiles,                    # Arrow files created earlier
  outputDirectory = "/home/kevin-koshy/Desktop/Academics/Single-Cell",  # Fixed output directory
  copyArrows = FALSE                           
)
print(proj)


cat("Dimension of peak set:\n")
print(dim(getPeakSet(proj)))

cat("Number of cells per sample:\n")
cellCounts <- table(proj$Sample)
print(cellCounts)

# Plot the fragment length distribution of all samples in a single plot
fragmentSizePlot <- plotFragmentSizes(proj, groupBy = "Sample")
ggsave("fragment_size_distribution.pdf", plot = fragmentSizePlot, width = 8, height = 6)

# Plot the TSS enrichment scores for each sample
tssEnrichmentPlot <- plotTSSEnrichment(proj, groupBy = "Sample")
ggsave("tss_enrichment.pdf", plot = tssEnrichmentPlot, width = 8, height = 6)

# Plot the number of fragments vs TSS enrichment for each sample
fragmentsVsTSSPlot <- ggplot(proj@cellColData, aes(x = log10(nFrags), y = TSSEnrichment, color = Sample)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Number of Fragments vs TSS Enrichment",
       x = "Log10(Number of Fragments)",
       y = "TSS Enrichment") +
  scale_color_manual(values = c("blue", "red", "green"))
ggsave("fragments_vs_tss.pdf", plot = fragmentsVsTSSPlot, width = 8, height = 6)


#Doublet Screening 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)
proj <- filterDoublets(ArchRProj = proj)
###############################################################################################################################
#WEEK2
#2 Dimensionality Reduction
#2.1 Iterative LSI
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI"
)

proj <- addUMAP(
 ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

#2.2 UMAP with sample annotation and QC metrics
p1 <- plotEmbedding(
 ArchRProj = proj,
 colorBy = "cellColData",
 name = "Sample",
 embedding = "UMAP"
)

p3 <- plotEmbedding(
 ArchRProj = proj,
 colorBy = "cellColData",
 name = "nFrags",
 embedding = "UMAP"
)

plotPDF(p1, p3, name = "UMAP_QC_Metrics.pdf", ArchRProj = proj, width = 5, height = 5)


#2.3 Dealing with batch effects
proj <- addHarmony(
 ArchRProj = proj,
 reducedDims = "IterativeLSI",
 name = "Harmony",
 groupBy = "Sample"
)

proj <- addUMAP(
 ArchRProj = proj,
 reducedDims = "Harmony",
 name = "UMAP_Harmony",
 nNeighbors = 30,
 minDist = 0.5
)

umapHarmony <- plotEmbedding(
 ArchRProj = proj,
 colorBy = "cellColData",
 name = "Sample",
 embedding = "UMAP_Harmony"
)

ggsave("UMAP_Harmony_BatchEffects.pdf", plot = umapHarmony, width = 8, height = 6)

#3 Clustering
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

umapClusters <- plotEmbedding(
ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

ggsave("UMAP_Clusters.pdf", plot = umapClusters, width = 8, height = 6)
clusterCounts <- table(proj$Clusters)
print(clusterCounts)

sampleProportions <- table(proj$Sample, proj$Clusters)
print(sampleProportions)

#PEAKS
#4.1 Peak calling

pathToMacs2 <- findMacs2()

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)

proj1 <- addPeakMatrix(proj, force = TRUE)

#cell types that we are working with 
table(proj1$Clusters)

#4.2 Cluster marker peaks

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj1, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)
###################################################################################################################
# WEEK3 
proj1 <- addImputeWeights(proj1)
markerGenes  <- c(
    "TOP2A", 
    "MKI67", 
    "AURKA", "SATB2", "SLC12A7" 
  )

p <- plotEmbedding(
    ArchRProj = proj1, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj1)
)
p_no_smooth <- plotEmbedding(
    ArchRProj = proj1, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = NULL  
)
# Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = "none", fill = "none") +  # Use "none" instead of FALSE
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()
    )
})
p_no_smooth2 <- lapply(p_no_smooth, function(x){
    x + guides(color = "none", fill = "none") +  
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()
    )
})


# Combine the two plots (with and without smoothing)
plot_grid_combined <- do.call(cowplot::plot_grid, c(list(ncol = 2), p_no_smooth2, p2))

# Save the combined plot to a file (e.g., PNG, PDF, etc.)
ggsave("markergenes_with_and_without_smoothing.png", plot = plot_grid_combined, width = 12, height = 8, dpi = 300)



#6 Transcription Factor motif activity
proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")

proj1 <- addBgdPeaks(proj1)  
proj1 <- addDeviationsMatrix(
  ArchRProj = proj1, 
  peakAnnotation = "Motif",
  force = TRUE
)

#6.2 Plot UMAP embeddings for marker TFs
motifVariability <- getVarDeviations(proj1, name = "MotifMatrix", plot = TRUE)
plotPDF(motifVariability, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)

#6.3 Motif activity
motifs <- c("GATA1", "CEBPA")
markerMotifs <- getFeatures(proj1, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

# Plot UMAP embeddings for marker motifs
p <- plotEmbedding(
    ArchRProj = proj1, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs),  
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj1)  
)

# Modify individual plots
p2 <- lapply(p, function(x) {
    x + guides(color = "none", fill = "none") +  
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()
    )
})

combined_plot <- do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
print(combined_plot)

# Retrieve the peak set
getPeakSet(proj)

ggsave(
    filename = "MarkerMotifs_UMAP.png",  # File name with extension
    plot = combined_plot,               # The plot object to save
    width = 12,                         # Width of the image in inches
    height = 4,                         # Height of the image in inches
    dpi = 300                           # Resolution in dots per inch

#7 Integration with gene expression
#7.1 Data integration

scRNA <- readRDS("/home/kevin-koshy/Desktop/Academics/Single-Cell/scbi_p2/scbi_p2/new_pbmc.rds")

colnames(scRNA@meta.data)

proj1 <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = scRNA,
    addToArrow = TRUE,
    groupRNA = "Cluster.Name",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
)

getAvailableMatrices(proj)

markerGenes <- c("ID4", "NFIA", "ASCL1") 

for (gene in markerGenes) {
    # Generate the UMAP plot
    p <- plotEmbedding(
        ArchRProj = proj1, 
        colorBy = "GeneIntegrationMatrix", 
        name = gene, 
        embedding = "UMAP"
    )

ggsave(filename = paste0("UMAP_", gene, ".png"), plot = p, width = 6, height = 5, dpi = 300)
print(paste("Saved UMAP plot for", gene, "to UMAP_", gene, ".png"))
}

#7.3 Cluster labels from gene expression
proj1$predictedGroup_Un <- factor(proj1$predictedGroup_Un)
cM <- as.matrix(confusionMatrix(proj1$Clusters, proj1$predictedGroup_Un))
# Assign labels based on the order
preClust <- colnames(cM) # Use the labels in order of appearance
assignedClusters <- rownames(cM) # Corresponding ATAC-seq clusters

# Combine the assignments into a table
assignments <- cbind(preClust, assignedClusters)

# Print the assignments
print(assignments)

#8 Peak-gene linkage
proj1 <- addPeak2GeneLinks(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj1,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

print(p2g)
p <- plotPeak2GeneHeatmap(ArchRProj = proj1, groupBy = "Clusters")

pdf("Peak2GeneHeatmap.pdf", width = 12, height = 8)
draw(p)
dev.off()
############################################################################################################################
#WEEK4
#9 Differential accessibility
#9.1 Differential peak accessibility
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C17",
  bgdGroups = "C9"
)


pma <- plotMarkers(
  seMarker = markerTest, 
  name = "C17", 
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", 
  plotAs = "MA"  # Generate MA plot
)
pma

pv <- plotMarkers(
  seMarker = markerTest, 
  name = "C17", 
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", 
  plotAs = "Volcano"  # Generate Volcano plot
)
pv

plotPDF(
  pma, pv, 
  name = "C17-MA-Volcano-Plots", 
  width = 5, 
  height = 5, 
  ArchRProj = proj1, 
  addDOC = FALSE
)
############################################################################################################################
#9.2 TF motif enrichment

# Add Motif Annotations
proj1 <- addMotifAnnotations(ArchRProj = proj1, motifSet = "cisbp", name = "Motif")

# Enrichment for GluN5
motifGluN5 <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj1,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

# Extract and rank GluN5 motifs
df <- data.frame(TF = rownames(motifGluN5), mlog10Padj = assay(motifGluN5)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
print("Top TF motifs for GluN5-specific peaks:")
head(df)

# Enrichment for Cyc. Prog.
motifsCycprog <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj1,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.2 & Log2FC <= -0.25"
)

# Extract and rank CycProg motifs
df2 <- data.frame(TF = rownames(motifsCycprog), mlog10Padj = assay(motifsCycprog)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))
print("Top TF motifs for Cyc. Prog.-specific peaks:")
head(df2)

# Plot heatmaps
heatmapGluN5 <- plotEnrichHeatmap(
    motifGluN5, 
    n = 10,  # Adjust to show top N motifs
    transpose = TRUE
)

heatmapCycprog <- plotEnrichHeatmap(
    motifsCycprog, 
    n = 10, 
    transpose = TRUE
)

# Draw and save heatmaps
ComplexHeatmap::draw(heatmapGluN5, heatmap_legend_side = "bot", annotation_legend_side = "bot")
ComplexHeatmap::draw(heatmapCycprog, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(
    heatmapGluN5, 
    name = "GluN5-Motif-Heatmap", 
    width = 8, 
    height = 6, 
    ArchRProj = proj1, 
    addDOC = FALSE
)

plotPDF(
    heatmapCycprog, 
    name = "CycProg-Motif-Heatmap", 
    width = 8, 
    height = 6, 
    ArchRProj = proj1, 
    addDOC = FALSE
)
#########################################################################################################################