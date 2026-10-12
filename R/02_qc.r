# Quality control for Visium HD data

source("R/00_setup.R")

visium <- readRDS(file.path(output_dir, "01_visium_hd_loaded.rds"))
qc_dir <- file.path(output_dir, "qc")

qc_before <- VlnPlot(
  visium,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
) +
  theme(axis.text.x = element_text(angle = 0))

save_plot(
  qc_before,
  file.path(qc_dir, "qc_metrics_before_filtering.png")
)

visium <- subset(
  visium,
  subset = nCount_RNA > config$analysis$min_counts &
    nFeature_RNA > config$analysis$min_genes &
    percent.mt < config$analysis$max_mito_percent
)

qc_after <- VlnPlot(
  visium,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
) +
  theme(axis.text.x = element_text(angle = 0))

save_plot(
  qc_after,
  file.path(qc_dir, "qc_metrics_after_filtering.png")
)

write_qc_summary(
  visium,
  file.path(qc_dir, "qc_summary.csv")
)

save_rds(
  visium,
  file.path(output_dir, "02_visium_hd_qc.rds")
)

message("Retained ", ncol(visium), " bins after QC.")
