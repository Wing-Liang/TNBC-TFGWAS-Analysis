#!/bin/bash
# ==============================================================================
# 08_tfgwas_combine.sbatch
# SLURM wrapper for 08_tfgwas_combine.R.
# Produces: Table S7
# Run     : sbatch scripts/05_replication/08_tfgwas_combine.sbatch   (from the repository root)
# ==============================================================================
#SBATCH --job-name=lmer_combine
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=replication/logs/slurm_lmercomb_%j.out
#SBATCH --error=replication/logs/slurm_lmercomb_%j.err

# Propagate the R exit status through the `| tee` pipeline. Without this the pipeline
# reports tee's status (0) and Slurm records COMPLETED even when Rscript failed, which
# silently releases an --dependency=afterok chain onto missing inputs.
set -o pipefail
cd "${TNBC_ROOT:?set TNBC_ROOT}"   # repository root; the R script switches to replication/
mkdir -p replication/logs
LOG=replication/logs/08_tfgwas_combine.log
{ echo "=== 08_tfgwas_combine.R  run $(date)  host=$(hostname)  job=$SLURM_JOB_ID ==="
  Rscript scripts/05_replication/08_tfgwas_combine.R
  rc=$?
  echo "=== exit=$rc end=$(date) ==="
  exit $rc
} 2>&1 | tee "$LOG"
exit "${PIPESTATUS[0]}"
