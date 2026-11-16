# Software environment

Most analyses were run in one conda environment (`signac_env`) with R 4.4.3, Signac 1.16.0,
Seurat 5.5.1, chromVAR 1.30.1, motifmatchr 1.30.0, JASPAR2022 0.99.8, TFBSTools 1.44.0,
lme4 2.0-6, lmerTest 3.2-1, coloc 5.2.3, BSgenome.Hsapiens.UCSC.hg38 1.4.5 / hg19 1.4.3,
EnsDb v86 / v75, and MACS2 2.2.9.1.

| File | Contents |
|---|---|
| `signac_env.yml` | exported `signac_env` conda environment (exact versions) |
| `sessionInfo.txt` | R `sessionInfo()` of `signac_env` with the main analysis packages loaded |
| `R_packages.txt` | every R package the scripts load |
| `requirements.txt` | Python packages |

Some steps use packages that are not part of `signac_env` and were run in separate
environments: motifbreakR, TCGA download and DESeq2 (TCGAbiolinks,
DESeq2), VIPER (viper, aracne.networks; R 4.5.1), and scDRS (Python). Install these from
Bioconductor, CRAN or PyPI as listed in `R_packages.txt` and `requirements.txt`.

To recreate the main environment:

    conda env create -f environment/signac_env.yml
