#!/bin/bash
# ==============================================================================
# 07b_tfgwas_model.sbatch
# SLURM array wrapper for 07b_tfgwas_model.R.
# Produces: Table S7
# Run     : sbatch scripts/05_replication/07b_tfgwas_model.sbatch   (from the repository root)
# ==============================================================================
#SBATCH --job-name=lmer_bymotif
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-37
#SBATCH --output=replication/logs/slurm_lmerM_%a_%j.out
#SBATCH --error=replication/logs/slurm_lmerM_%a_%j.err
# Set --array to match tfgwas/lmer_batches/n_batch.txt (07a prints it).
# Submit:  mkdir -p replication/logs; JA=$(sbatch --parsable scripts/05_replication/07b_tfgwas_model.sbatch)
# Combine: sbatch --dependency=afterok:${JA} scripts/05_replication/08_tfgwas_combine.sbatch
set -euo pipefail
cd "${TNBC_ROOT:?set TNBC_ROOT}"   # repository root; the R script switches to replication/
mkdir -p replication/logs
echo "=== 07b_tfgwas_model.R  BATCH=${SLURM_ARRAY_TASK_ID}  $(date)  host=$(hostname) ==="
Rscript scripts/05_replication/07b_tfgwas_model.R ${SLURM_ARRAY_TASK_ID}
echo "=== exit=$? end=$(date) ==="
