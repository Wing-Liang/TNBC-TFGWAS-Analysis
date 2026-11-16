#!/bin/bash
# ==============================================================================
# 01_download_geo_data.sh
# Download GSE244585 pools 9, 11 and 12 (European ancestry) from GEO.
# Produces: replication input
# Run     : bash scripts/05_replication/01_download_geo_data.sh   (from the repository root)
# ==============================================================================
# 01 - Download GEO data for the three European (non-Ashkenazi) pools of GSE244585
#      Pools 9, 11, 12  (GEO "White American" = paper "European-non-Ashkenazi")
# Genome build of source data: GRCh38/hg38
set -u
cd "${TNBC_ROOT:?set TNBC_ROOT}/replication"
SAMPLE_BASE='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7821nnn'
SERIES_BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244585/suppl'

# ATAC fragments (+ tabix index)  GSM: pool9=7821199 pool11=7821200 pool12=7821201
for spec in "GSM7821199 KTBpool9" "GSM7821200 KTBpool11" "GSM7821201 KTBpool12"; do
  set -- $spec; gsm=$1; nm=$2
  wget -c "$SAMPLE_BASE/$gsm/suppl/${gsm}_${nm}_atac_fragments.tsv.gz"
  wget -c "$SAMPLE_BASE/$gsm/suppl/${gsm}_${nm}_atac_fragments.tsv.gz.tbi.gz"
  gunzip -kf "${gsm}_${nm}_atac_fragments.tsv.gz.tbi.gz"     # Signac needs .tbi (not .tbi.gz)
  gzip -t "${gsm}_${nm}_atac_fragments.tsv.gz" && echo "$nm fragments OK" || echo "$nm fragments TRUNCATED"
done

# RNA multiome h5   GSM: pool9=7821215 pool11=7821216 pool12=7821217
for spec in "GSM7821215 KTBpool9" "GSM7821216 KTBpool11" "GSM7821217 KTBpool12"; do
  set -- $spec; gsm=$1; nm=$2
  wget -c "$SAMPLE_BASE/$gsm/suppl/${gsm}_${nm}_filtered_feature_bc_matrix.h5"
done
# NOTE: the pool11 h5 (GSM7821216) is truncated at the GEO source; RNA counts for all three pools
# are therefore taken from the CELLxGENE object below (03a_extract_cellxgene_data.py).

# per-barcode metrics (series-level) -> source of is_cell cell calls
mkdir -p data/geo_cells data/fragments_hg38 data/cellxgene
for p in 9 11 12; do
  wget -c "$SERIES_BASE/GSE244585_KTBpool${p}_per_barcode_metrics.csv.gz" -O "data/geo_cells/GSE244585_KTBpool${p}_per_barcode_metrics.csv.gz"
  gzip -t "data/geo_cells/GSE244585_KTBpool${p}_per_barcode_metrics.csv.gz" && echo "pool${p} metrics OK"
done

# CELLxGENE object of the source atlas (Bhat-Nakshatri et al., Nat Med 2024): author cell types,
# raw RNA counts and embeddings (dataset 61af564d-e5ea-4d34-a0f3-2668a00db376)
wget -c "https://datasets.cellxgene.cziscience.com/31798b24-255e-4726-8805-66507589ccb8.h5ad" -O data/cellxgene/breast_atlas.h5ad
