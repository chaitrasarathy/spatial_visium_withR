# Cell-type deconvolution with RCTD

source("R/00_setup.R")

visium <- readRDS(file.path(output_dir, "04_visium_hd_spatial.rds"))
decon_dir <- file.path(output_dir, "deconvolution")

reference_path <- config$paths$reference_rds
cell_type_col <- config$paths$reference_celltype_column

if (!file.exists(reference_path)) {
  warning(
    "Reference RDS not found: ", reference_path,
    ". Skipping deconvolution. Add an annotated scRNA-seq reference and update config.yaml."
  )
  quit(save = "no", status = 0)
}

reference <- readRDS(reference_path)
Idents(reference) <- cell_type_col

reference <- subset(
  reference,
  idents = names(which(table(Idents(reference)) >= config$deconvolution$min_reference_cells_per_type))
)

message("Running RCTD deconvolution...")

rctd_result <- run.RCTD(
  spatialRNA = visium,
  reference = reference,
  max_cores = config$deconvolution$max_cores,
  UMI_min = config$analysis$min_counts,
  verbose = FALSE
)

# Extract cell-type proportions
proportions <- as.data.frame(rctd_result@results$weights)

write.csv(
  proportions,
  file.path(decon_dir, "rctd_cell_type_proportions.csv")
)

# Add dominant cell type to the Seurat object
if (ncol(proportions) > 0) {
  visium$RCTD_dominant_celltype <- colnames(proportions)[
    max.col(proportions, ties.method = "first")
  ]

  p <- SpatialDimPlot(
    visium,
    group.by = "RCTD_dominant_celltype",
    label = TRUE,
    label.size = 3
  )

  save_plot(
    p,
    file.path(decon_dir, "rctd_dominant_celltype_map.png"),
    width = 10,
    height = 8
  )
}

save_rds(
  visium,
  file.path(output_dir, "05_visium_hd_deconvoluted.rds")
)

message("RCTD deconvolution complete.")
