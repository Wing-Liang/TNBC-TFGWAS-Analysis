# Single-cell processing and annotation

Quality control, peak calling, cell-type annotation and RNA-ATAC integration of joint scRNA-seq and scATAC-seq from normal breast tissue of four donors (GEO GSE162798).

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_import_geo_data.R` | Import the GSE162798 scRNA-seq and scATAC-seq matrices, fragments and cell metadata. | input to all single-cell steps |
| `02_call_peaks_macs2.R` | Call scATAC-seq peaks with MACS2 (Signac CallPeaks), discovery cohort. | peak set (Fig 2E, Table S1) |
| `03_qc_doublets.R` | ATAC and RNA quality control and doublet removal (scDblFinder). | final cell set |
| `04_annotate_cell_types.R` | Cluster the scRNA-seq data and annotate cell types from canonical markers. | Fig 2A |
| `05_atac_umap_figures.R` | scATAC-seq dimensionality reduction, chromVAR motif activity and Fig 2 panels. | Fig 2C-E |
| `06_wnn_integration.R` | Weighted-nearest-neighbour integration of RNA and ATAC. | Fig 2 |
| `07_markers_and_da_peaks.R` | Differentially accessible peaks, peak annotation and top-5 marker genes per cell type. | Fig 2F, Table S1; DA peaks used for Fig 3F |
| `08_marker_plots.R` | Canonical marker feature plots for the scRNA-seq cell types. | Fig 2B |
| `09_cell_counts.R` | Number of cells per cell type. | Fig 2D |
| `10_motif_annotation_chromvar.R` | JASPAR 2022 motif annotation, chromVAR activity, differential motif activity and motif enrichment in cell-type-unique DA peaks. | input to Figs 3, S2-S8; motif set of Table S4 |
| `11_export_umap_embedding.R` | Export RNA UMAP coordinates and cell types of the post-QC cells. | input to Figs 3D, S9 |
