# Utility functions for the Visium HD pipeline

suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(ggplot2)
})

#' Read pipeline configuration
read_config <- function(path = "config/config.yaml") {
  stopifnot(file.exists(path))
  yaml::read_yaml(path)
}

#' Ensure required directories exist
ensure_dirs <- function(output_dir) {
  dirs <- file.path(
    output_dir,
    c("qc", "clustering", "spatial", "deconvolution", "logs")
  )
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dirs)
}

#' Save a ggplot figure
save_plot <- function(plot, path, width = 9, height = 7, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = width, height = height, dpi = dpi)
}

#' Write a simple QC summary table
write_qc_summary <- function(object, path) {
  summary_df <- data.frame(
    metric = c(
      "n_cells",
      "n_genes",
      "median_counts_per_cell",
      "median_genes_per_cell",
      "median_mito_percent"
    ),
    value = c(
      ncol(object),
      nrow(object),
      median(object$nCount_RNA),
      median(object$nFeature_RNA),
      median(object$percent.mt)
    )
  )
  write.csv(summary_df, path, row.names = FALSE)
  invisible(summary_df)
}

#' Save Seurat object safely
save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
}
