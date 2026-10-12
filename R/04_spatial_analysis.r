# Marker-gene annotation and spatial analysis

source("R/00_setup.R")

visium <- readRDS(file.path(output_dir, "03_visium_hd_clustered.rds"))
spatial_dir <- file.path(output_dir, "spatial")

# Spatially variable features
visium <- SpatialFeaturePlot(
  visium,
  features = c("nCount_RNA"),
  alpha = 0.8
)

# Marker-gene spatial expression
marker_list <- config$marker_genes

for (domain in names(marker_list)) {
  genes <- marker_list[[domain]]
  genes <- intersect(genes, rownames(visium))

  if (length(genes) == 0) {
    warning("No marker genes found for domain: ", domain)
    next
  }

  p <- SpatialFeaturePlot(
    visium,
    features = genes[seq_len(min(4, length(genes)))],
    alpha = 0.8,
    ncol = 2
  )

  save_plot(
    p,
    file.path(spatial_dir, paste0("markers_", domain, ".png")),
    width = 10,
    height = 8
  )
}

# Neighborhood enrichment among clusters
visium <- FindNeighbors(
  visium,
  reduction = "pca",
  dims = 1:config$analysis$pca_dims,
  verbose = FALSE
)

# A simple cluster-co-occurrence summary can be extended with Squidpy
# or Giotto if Python/R interoperability is added later.
cluster_table <- table(
  visium$seurat_clusters
)
write.csv(
  as.data.frame(cluster_table),
  file.path(spatial_dir, "cluster_bin_counts.csv"),
  row.names = FALSE
)

save_rds(
  visium,
  file.path(output_dir, "04_visium_hd_spatial.rds")
)

message("Spatial marker and neighborhood analysis complete.")
