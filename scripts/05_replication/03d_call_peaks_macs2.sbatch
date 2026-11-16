#!/bin/bash
# 03d_call_peaks_macs2.sbatch
# One MACS2 peak call on the three pools' ATAC fragments (Signac CallPeaks parameters).
# 04_build_combined_object.R imports peaks/macs2/combined_peaks.narrowPeak and applies the
# autosome and ENCODE-blacklist filters.
# Run    : sbatch scripts/05_replication/03d_call_peaks_macs2.sbatch   (or bash; from the repository root)
#SBATCH --job-name=peaks_comb
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=replication/logs/slurm_peaks_%j.out
#SBATCH --error=replication/logs/slurm_peaks_%j.err
cd "${TNBC_ROOT:?set TNBC_ROOT}/replication"
LOG=logs/03d_call_peaks_macs2.log
MACS2=macs2
mkdir -p peaks/macs2 logs
{
  echo "=== MACS2 peak calling, pools 9/11/12  $(date) ==="
  "$MACS2" callpeak \
     -t GSM7821199_KTBpool9_atac_fragments.tsv.gz \
        data/fragments_hg38/pool11_atac_fragments.hg38.tsv.gz \
        GSM7821201_KTBpool12_atac_fragments.tsv.gz \
     -g 2.7e+09 -f BED --nomodel --extsize 200 --shift -100 \
     -n combined --outdir peaks/macs2
  echo "=== MACS2 exit=$?  $(date) ==="
} 2>&1 | tee "$LOG"
