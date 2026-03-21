## ============================================================================
## Visium Spatial Transcriptomics pipeline - Mouse Brain Coronal Section
## End-to-end analysis: QC → SCTransform → Clustering → Spatial Variable Genes
## Requires: Seurat v5, SeuratData, ggplot2, patchwork
## ============================================================================

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

dir.create("figures", showWarnings = FALSE)

## --- Load data -----------------------------------------------------------

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

## --- QC ---------------------------------------------------------------------

brain[["percent_mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")

p_qc_counts <- SpatialFeaturePlot(brain, features = "nCount_Spatial") +
  theme(legend.position = "right")
p_qc_genes  <- SpatialFeaturePlot(brain, features = "nFeature_Spatial") +
  theme(legend.position = "right")
p_qc_mt     <- SpatialFeaturePlot(brain, features = "percent_mt") +
  theme(legend.position = "right")

ggsave("figures/qc_spatial.png", p_qc_counts | p_qc_genes | p_qc_mt,
       width = 18, height = 6, dpi = 300)

brain <- subset(brain,
                nFeature_Spatial > 200 &
                nCount_Spatial > 500 &
                percent_mt < 25)

## --- Normalisation (SCTransform) --------------------------------------------

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

## --- Dimensionality reduction & clustering ----------------------------------

brain <- RunPCA(brain, verbose = FALSE)
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 0.8, verbose = FALSE)
brain <- RunUMAP(brain, dims = 1:30, verbose = FALSE)

p_umap <- DimPlot(brain, reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend()
ggsave("figures/umap_clusters.png", p_umap, width = 7, height = 6, dpi = 300)

p_spatial_clusters <- SpatialDimPlot(brain, label = TRUE, label.size = 3) +
  theme(legend.position = "right")
ggsave("figures/spatial_clusters.png", p_spatial_clusters,
       width = 8, height = 7, dpi = 300)

## --- Spatially variable gene detection --------------------------------------

brain <- FindSpatiallyVariableFeatures(brain,
                                       assay = "SCT",
                                       selection.method = "moransi",
                                       features = VariableFeatures(brain),
                                       verbose = FALSE)

sv_genes <- SpatiallyVariableFeatures(brain, selection.method = "moransi")
top_sv   <- head(sv_genes, 6)

p_sv <- SpatialFeaturePlot(brain, features = top_sv, ncol = 3) &
  theme(legend.position = "right")
ggsave("figures/top_spatially_variable_genes.png", p_sv,
       width = 18, height = 12, dpi = 300)

## --- Known marker genes on tissue -------------------------------------------

markers <- c("Hpca", "Ttr", "Mbp")
markers_present <- markers[markers %in% rownames(brain)]

if (length(markers_present) > 0) {
  p_markers <- SpatialFeaturePlot(brain, features = markers_present, ncol = 3) &
    theme(legend.position = "right")
  ggsave("figures/known_markers_spatial.png", p_markers,
         width = 6 * length(markers_present), height = 6, dpi = 300)
}

## --- Session info -----------------------------------------------------------

writeLines(capture.output(sessionInfo()), "figures/session_info.txt")

cat("Done. Figures written to figures/\n")
