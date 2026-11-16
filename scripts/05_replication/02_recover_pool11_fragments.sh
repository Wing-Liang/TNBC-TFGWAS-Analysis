#!/bin/bash
# ==============================================================================
# 02_recover_pool11_fragments.sh
# Recover the truncated pool-11 fragment file deposited on GEO.
# Produces: replication input
# Run     : bash scripts/05_replication/02_recover_pool11_fragments.sh   (from the repository root)
# ==============================================================================
cd "${TNBC_ROOT:?set TNBC_ROOT}/replication"
LOG=data/fragments_hg38/pool11_prep.log
echo "[$(date)] pool11: tolerant decompress -> bgzip (cellranger-sorted, no sort needed)" > $LOG
zcat GSM7821200_KTBpool11_atac_fragments.tsv.gz 2>>$LOG | bgzip > data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz
echo "[$(date)] pool11: tabix" >> $LOG
tabix -p bed data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz
n=$(zcat data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz | wc -l)
echo "[$(date)] pool11: DONE fragments=$n" >> $LOG
tabix -l data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz | head -3 | tr '\n' ' ' >> $LOG
echo "" >> $LOG
