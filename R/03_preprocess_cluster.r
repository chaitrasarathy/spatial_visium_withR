# Preprocessing, dimensionality reduction, and clustering

source("R/00_setup.R")

visium <- readRDS(file.path(output_dir, "02_visium_hd_qc.rds"))
cluster_dir <- file.path(output_dir, "clustering")

if (config$analysis$normalize_method == "SCT") {
  visium <- SCTransform(
    visium,
    assay = "Spatial",
    clip.range = c(-10, 10),
    verbose = FALSE
  )
} else {
  visium <- NormalizeData(visium, verbose = FALSE)
  visium <- FindVariableFeatures(
    visium,
    nfeatures = config$analysis$n_variable_features,
    verbose = FALSE
  )
  visium <- ScaleData(visium, verbose = FALSE)
}

visium <- RunPCA(visium, npcs = config$analysis$pca_dims, verbose = FALSE)
visium <- RunUMAP(
  visium,
  dims = 1:config$analysis$pca_dims,
  verbose = FALSE
)
visium <- FindNeighbors(visium, dims = 1:config$analysis$pca_dims, verbose = FALSE)
visium <- FindClusters(
  visium,
  resolution = config$analysis$cluster_resolution,
  verbose = FALSE
)

p_umap <- DimPlot(visium, reduction = "umap", group.by = "seurat_clusters")
save_plot(p_umap, file.path(cluster_dir, "umap_clusters.png"))

p_spatial <- SpatialDimPlot(
  visium,
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3
)
save_plot(p_spatial, file.path(cluster_dir, "spatial_clusters.png"), width = 10, height = 8)

save_rds(
  visium,
  file.path(output_dir, "03_visium_hd_clustered.rds")
)

message("Clustering complete.")
