# Software environment

Analyses were run in a conda environment with R 4.4.3, Signac 1.16.0, Seurat 5.5.1,
chromVAR 1.30.1, motifmatchr 1.30.0, JASPAR2022, TFBSTools 1.44.0,
BSgenome.Hsapiens.UCSC.hg38 1.4.5 / hg19 1.4.3, EnsDb v86 / v75, and MACS2 2.2.9.1.

| File | Contents |
|---|---|
| `R_packages.txt` | every R package the scripts load |
| `requirements.txt` | Python packages |
| `signac_env.yml` | exported conda environment (exact versions) |
| `sessionInfo.txt` | R `sessionInfo()` of the analysis environment |

To recreate the environment:

    conda env create -f environment/signac_env.yml
