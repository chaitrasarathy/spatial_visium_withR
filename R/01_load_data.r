# Load Visium HD Space Ranger outputs

source("R/00_setup.R")

space_ranger_dir <- config$paths$space_ranger_dir
bin_size <- config$analysis$bin_size_um

stopifnot(dir.exists(space_ranger_dir))

message("Loading Visium HD data from: ", space_ranger_dir)

# For Visium HD, 10x commonly provides binned outputs such as:
# binned_outputs/square_016um/
# Adjust this path if your Space Ranger output uses another structure.
binned_dir <- file.path(
  space_ranger_dir,
  "binned_outputs",
  sprintf("square_%03dum", bin_size)
)

if (dir.exists(binned_dir)) {
  data_dir <- binned_dir
} else {
  data_dir <- space_ranger_dir
}

visium <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = config$sample$id,
  filter.matrix = TRUE,
  to.upper = FALSE
)

DefaultAssay(visium) <- "Spatial"

visium$percent.mt <- PercentageFeatureSet(
  visium,
  pattern = "^MT-"
)

save_rds(
  visium,
  file.path(output_dir, "01_visium_hd_loaded.rds")
)

message("Loaded object: ")
print(visium)
