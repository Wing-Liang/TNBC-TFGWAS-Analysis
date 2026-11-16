# TNBC TF-GWAS

Analysis code for:

> Liang Y, Wang X, Lee M, Yang Y, Cai Q, Long J, Chandarlapaty S, Guo X, Buas MF, Shu X.
> **Integration of single-cell multi-omics data with GWAS reveals key transcriptional regulators
> underlying susceptibility to triple-negative breast cancer.** Manuscript in preparation (2026).
> DOI: _to be added on publication_

## Overview

The TF-GWAS framework tests, for each of 727 JASPAR 2022 transcription-factor (TF) motifs,
whether GWAS variants in accessible chromatin carrying the motif show stronger association with
breast-cancer risk:

```
{subtype}_chisquare ~ motif_in_peak + (1 | LD block)        (lme4 / lmerTest)
```

fitted genome-wide for overall breast cancer and five intrinsic subtypes, using peaks from joint
scRNA-/scATAC-seq of normal breast tissue, and repeated in an independent single-cell cohort.

| Folder (`scripts/`) | Analysis | Figures / tables |
|---|---|---|
| `01_singlecell_annotation` | scRNA-/scATAC-seq QC, peak calling, cell-type annotation, motif annotation and chromVAR (GSE162798) | Fig 2, Table S1 |
| `02_gwas_ld_blocks` | BCAC subtype GWAS filtering; variants lifted GRCh37 → GRCh38 and assigned to 1,336 European LD blocks | Table S2 |
| `03_tf_gwas` | Variant × peak × motif matrix; TF-GWAS model for overall BC and five subtypes | Fig 3A/C, Tables S3–S4 |
| `04_celltype_enrichment` | Subtype-specific motifs; cell-type enrichment (AddModuleScore, chromVAR) | Fig 3B/D–F, Figs S2–S8, Table S5 |
| `05_replication` | Independent replication (GSE244585, 15 donors, 10x Multiome); TF-family grouping | Tables S6–S8 |
| `06_dissociation_ap1` | Dissociation-stress and AP-1 scores; cohort composition | Fig S9 |
| `07_scdrs` | Single-cell disease relevance score | Fig S10 |
| `08_tcga_viper_depmap` | TCGA-BRCA DESeq2; msVIPER TF activity; DepMap 25Q3 gene effects | Fig 4, Tables S9–S11, S14 |
| `09_motifbreakR_gtex_coloc` | motifbreakR; colocalisation with GTEx v8 breast cis-eQTL; GTEx v11 genotype–expression | Fig 5, Fig S11, Tables S12–S13 |

Each folder's `README.md` lists its scripts in run order with what each one produces;
[`docs/FIGURES_AND_TABLES.md`](docs/FIGURES_AND_TABLES.md) maps every figure panel and table to
its script. Every script starts with a header giving its purpose, inputs, outputs and run command.

## Running the analysis

Run every script from the repository root. Inputs are read from `data/`, outputs are written to
`results/`; the replication cohort uses its own working folder `replication/` (with its raw
inputs in `replication/data/`). None of these folders is tracked by git.

```
TNBC-TFGWAS-Analysis/
├── data/                 public inputs (see "Data" below)
├── replication/          GSE244585 working folder: data/, qc/, peaks/, tfgwas/, ...
├── results/              all discovery-cohort outputs
│   ├── singlecell/   ld_blocks/   gwas/   tfgwas/   motif_sets/   motifbreakR/
│   └── coloc/   target_genes/   tcga/   tf_activity/   scdrs/
├── reference/            small lookup tables shipped with the code
└── scripts/              analysis code, folders 01–09
```

Folders can be relocated with environment variables:

| Variable | Contents | Default |
|---|---|---|
| `TNBC_ROOT` | project root (contains `data/`, `results/`, `replication/`) | `.` |
| `BCAC_DIR` | BCAC/CIMBA subtype and overall breast-cancer GWAS summary statistics | `data/BCAC` |
| `LDBLOCK_DIR` | `pyrho_EUR_LD_blocks.bed`, `hg19ToHg38.over.chain` | `data/LDblocks` |
| `LIFTOVER_DIR` | `hg19ToHg38.over.chain` | `data/liftover` |
| `KG_DIR` | 1000 Genomes EUR reference: PLINK panel `1kg.v3/EUR` (motifbreakR) and `1000G.EUR.QC.plink.simple.vcf` (rsIDs for MAGMA) | `data/1KG` |
| `GTEX_DIR` | GTEx / eQTL Catalogue files | `data/GTEx` |
| `TCGA_DIR` | TCGA-BRCA expression object and GENCODE v36 annotation | `data/TCGA` |
| `DEPMAP_DIR` | DepMap `CRISPRGeneEffect_25Q3.csv`, `Model.csv` | `data/DepMap` |
| `MAGMA_DIR` | MAGMA inputs and gene-level results (scDRS) | `data/MAGMA` |
| `LDLINK_TOKEN` | personal LDlink API token (locus plots only) | — |

Run-time options used by individual scripts: `VIPER_MINSIZE` (default 5, as reported),
`MOTIF_SET`, `OVERALLBC_COHORT`, `GROUP_REF` (default `reference/`) and `OBS_RESULTS`; see each
script's header.

## Requirements

- R 4.4.3 with Signac 1.16.0, Seurat 5.5.1, chromVAR 1.30.1, motifmatchr 1.30.0, JASPAR2022,
  TFBSTools 1.44.0, lme4, lmerTest, coloc, motifbreakR and data.table; the TCGA, VIPER and DepMap
  scripts (`08_tcga_viper_depmap`) were run with R 4.5.1, DESeq2 1.48.2 and viper 1.44.0.
  Full list: `environment/R_packages.txt`.
- MACS2 2.2.9.1, PLINK 1.9 (LD clumping), MAGMA and scDRS 1.0.3
- Python ≥ 3.9 with h5py, numpy, pandas, matplotlib, openpyxl (`environment/requirements.txt`)
- SLURM for the `*.sbatch` wrappers (optional; set the partition line for your cluster)

## Data

No data are included in this repository. Download each dataset from its source and save it in
the local folder given in the last column (the defaults can be changed with the environment
variables above):

| Data | Source | Local folder |
|---|---|---|
| Discovery scRNA-/scATAC-seq, 4 donors | GEO [GSE162798](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162798) | `data/GSE162798/` |
| Replication 10x Multiome, pools 9/11/12, 15 donors | GEO [GSE244585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244585); CELLxGENE collection `77446b76-1c2d-4a71-8e59-0efd4374d98e` | `replication/data/` (downloaded by `05_replication/01_download_geo_data.sh`) |
| Breast-cancer subtype GWAS (`icogs_onco_meta_intrinsic_subtypes_summary_stats_zhang_2020.zip`) and overall-BC GWAS (Zhang et al. 2020) | Breast Cancer Association Consortium (BCAC), [GWAS summary associations](https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations) | `$BCAC_DIR` |
| LD blocks (pyrho EUR, GRCh38) | [jmacdon/LDblocks_GRCh38](https://github.com/jmacdon/LDblocks_GRCh38) | `$LDBLOCK_DIR` |
| GRCh37 → GRCh38 chain | UCSC `hg19ToHg38.over.chain.gz` | `$LDBLOCK_DIR`, `$LIFTOVER_DIR` |
| 1000 Genomes Phase 3 EUR panel | 1000 Genomes Project | `$KG_DIR/1kg.v3/` |
| 1000 Genomes EUR SNP reference VCF (rsIDs) | [scPagwas](https://github.com/sulab-wmu/scPagwas) reference data, `1000G.EUR.QC.plink.simple.vcf` | `$KG_DIR` |
| GTEx v8 breast cis-eQTL summary statistics | eQTL Catalogue `QTD000211` | `$GTEX_DIR/eQTL_v8/` |
| GTEx v11 breast expression and v9 WGS genotypes (Fig 5D–H) | dbGaP `phs000424` (controlled access) | `$GTEX_DIR/v11/` |
| TCGA-BRCA RNA-seq and PAM50 calls | Genomic Data Commons (Data Release 44.0), downloaded by `08_tcga_viper_depmap/01_download_tcga_brca.R` | `results/tcga/`, `$TCGA_DIR` |
| CRISPR gene effects | DepMap Public 25Q3 | `$DEPMAP_DIR` |
| TF motifs | JASPAR 2022 (`JASPAR2022` Bioconductor package) | — |
| LD for locus plots | LDlink (`LDlinkR`) | — |

`reference/` holds small lookup tables used by the scripts (the 727 tested motif IDs, motif →
TF family and JASPAR 2022 cluster assignments, promoter-SNP target genes), derived from
JASPAR 2022 (CC BY 4.0).

## Citation

If you use this code, please cite the article above (see also `CITATION.cff`).

## License

MIT — see `LICENSE`.
