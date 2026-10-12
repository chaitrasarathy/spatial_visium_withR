# Methods

## Data processing

Visium HD Space Ranger outputs are loaded into R using `Seurat::Load10X_Spatial()`.
For high-resolution Visium HD data, binned outputs are used where available.

## Quality control

Bins are filtered using configurable thresholds for:

- Minimum UMI counts
- Minimum detected genes
- Maximum mitochondrial transcript percentage

QC distributions are plotted before and after filtering.

## Normalization and clustering

The default workflow uses SCTransform normalization, PCA, UMAP, nearest-neighbor graph construction, and Leiden/Louvain-style clustering through Seurat.

## Spatial analysis

Marker-gene expression is projected onto tissue space to support annotation of tumor epithelium, stroma, immune, endothelial, goblet-cell, and smooth-muscle domains.

The pipeline is structured to support later extension to spatially variable gene testing, neighborhood enrichment, and tissue-domain segmentation.

## Cell-type deconvolution

RCTD is used to estimate cell-type composition at each spatial bin using an annotated single-cell RNA-seq reference.

The dominant inferred cell type is visualized spatially, while full proportion estimates are exported as a CSV file.

## Reporting

The final report is generated with Quarto and includes methods, QC, figures, interpretation, limitations, and session information.
