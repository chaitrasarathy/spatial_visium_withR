# Load packages and configure the Visium HD pipeline

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(yaml)
  library(future)
})

source("R/utils.R")

config <- read_config("config/config.yaml")
output_dir <- config$paths$output_dir
ensure_dirs(output_dir)

plan(multisession, workers = config$deconvolution$max_cores)

message("Pipeline configuration loaded for sample: ", config$sample$id)
