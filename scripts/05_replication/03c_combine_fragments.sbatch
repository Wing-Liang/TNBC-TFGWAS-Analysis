#!/bin/bash
# 03c_combine_fragments.sbatch
# Combine the three pools' ATAC fragments into one file: barcodes prefixed with the pool,
# CELLxGENE cells only, autosomes chr1-22 only, fully re-sorted and tabix-indexed.
# Inputs : GEO fragment files (01), data/fragments_hg38/pool11 (02), data/cellxgene/pool*_cells.txt (03a)
# Output : combined/combined_fragments.tsv.gz (+ .tbi)
# Run    : sbatch scripts/05_replication/03c_combine_fragments.sbatch   (or bash; from the repository root)
#SBATCH --job-name=combine_frag
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=replication/logs/slurm_combine_%j.out
#SBATCH --error=replication/logs/slurm_combine_%j.err
set -u
cd "${TNBC_ROOT:?set TNBC_ROOT}/replication"
mkdir -p combined logs
export LC_ALL=C                    # deterministic byte-order collation for sort/tabix
LOG=logs/03c_combine_fragments.log
{
echo "[$(date)] combine 3 pools -> ONE fragment file: pool-prefixed barcodes, VALID CELLS ONLY,"
echo "          AUTOSOMES chr1-22 only (drops chrX/Y + alt contigs), FULL re-sort (not merge) for tabix"
TMP=${TMPDIR:-/tmp}/combsort_${SLURM_JOB_ID:-$$}; mkdir -p "$TMP"
# each stream: keep chr1-22 + valid CELLxGENE barcodes, prefix barcode with pool id.
# NOTE: full `sort` (no -m) is required -- inputs are in cellranger order (chr1,chr10,...,chrY,alts)
# which is NOT monotonic under the sort key, so a merge would interleave chromosomes.
AUTO='^chr([1-9]|1[0-9]|2[0-2])$'
sort -k1,1 -k2,2n --parallel=8 -T "$TMP" -S 40G \
  <(zcat GSM7821199_KTBpool9_atac_fragments.tsv.gz    | awk -v A="$AUTO" 'BEGIN{OFS="\t"; while((getline b < "data/cellxgene/pool9_cells.txt")>0)  k[b]=1} /^#/{next} ($1 ~ A && $4 in k){$4="pool9_"$4;  print}') \
  <(zcat data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz | awk -v A="$AUTO" 'BEGIN{OFS="\t"; while((getline b < "data/cellxgene/pool11_cells.txt")>0) k[b]=1} /^#/{next} ($1 ~ A && $4 in k){$4="pool11_"$4; print}') \
  <(zcat GSM7821201_KTBpool12_atac_fragments.tsv.gz   | awk -v A="$AUTO" 'BEGIN{OFS="\t"; while((getline b < "data/cellxgene/pool12_cells.txt")>0) k[b]=1} /^#/{next} ($1 ~ A && $4 in k){$4="pool12_"$4; print}') \
  | bgzip > combined/combined_fragments.tsv.gz
echo "[$(date)] tabix index"
tabix -p bed combined/combined_fragments.tsv.gz
n=$(zcat combined/combined_fragments.tsv.gz | wc -l)
echo "[$(date)] DONE: $n valid-cell autosomal fragments in combined file"
echo -n "chroms: "; tabix -l combined/combined_fragments.tsv.gz | tr '\n' ' '; echo
rm -rf "$TMP"
} 2>&1 | tee "$LOG"
