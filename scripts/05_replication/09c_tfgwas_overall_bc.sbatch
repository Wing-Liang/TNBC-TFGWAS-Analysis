#!/bin/bash
# ==============================================================================
# 09c_tfgwas_overall_bc.sbatch
# SLURM array wrapper for 09c_tfgwas_overall_bc.R, replication cohort.
# Produces: Table S7 (overall BC)
# Run     : sbatch scripts/05_replication/09c_tfgwas_overall_bc.sbatch   (from the repository root)
# ==============================================================================
#SBATCH --job-name=lmer_OverallBC
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-37
#SBATCH --output=replication/logs/slurm_OverallBC_%a_%j.out
#SBATCH --error=replication/logs/slurm_OverallBC_%a_%j.err
set -euo pipefail
cd "${TNBC_ROOT:?set TNBC_ROOT}"   # repository root; the R script switches to replication/
mkdir -p replication/logs
echo "=== 09c_tfgwas_overall_bc.R  BATCH=${SLURM_ARRAY_TASK_ID}  $(date)  host=$(hostname) ==="
Rscript scripts/05_replication/09c_tfgwas_overall_bc.R ${SLURM_ARRAY_TASK_ID}
echo "=== exit=$? end=$(date) ==="
