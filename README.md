# Spatial Transcriptomics Analysis — 10x Visium Mouse Brain
End-to-end spatial transcriptomics workflow on the 10x Visium mouse brain anterior coronal section. Analysis covers quality control, SCTransform normalisation, unsupervised clustering, spatially variable gene detection (Moran's I), and projection of clusters and marker genes onto tissue. Built with Seurat v5.

# Data
The 10x Visium mouse brain anterior coronal section is loaded via the SeuratData package (stxBrain dataset). Raw data are also available directly from 10x Genomics. No data files are committed to this repository.

# Reproducibility
Requirements: R ≥ 4.3, Seurat ≥ 5.0, SeuratData, ggplot2, patchwork
```bash
# Install the dataset (once)
Rscript -e 'SeuratData::InstallData("stxBrain")'

# Run the full analysis
Rscript analysis/visium_mouse_brain_analysis.R
```

# Run the full analysis
Rscript analysis/visium_mouse_brain_analysis.R
All figures are written to figures/. A full sessionInfo() is saved to figures/session_info.txt.

# Repository structure
```
├── README.md
├── analysis/
│   └── visium_mouse_brain_analysis.R
├── figures/
│   ├── qc_spatial.png
│   ├── umap_clusters.png
│   ├── spatial_clusters.png
│   ├── top_spatially_variable_genes.png
│   ├── known_markers_spatial.png
│   └── session_info.txt
├── .gitignore
└── LICENSE
```
