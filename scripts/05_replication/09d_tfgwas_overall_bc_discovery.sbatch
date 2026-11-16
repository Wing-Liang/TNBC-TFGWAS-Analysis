#!/bin/bash
# ==============================================================================
# 09d_tfgwas_overall_bc_discovery.sbatch
# SLURM array wrapper for 09c_tfgwas_overall_bc.R, discovery cohort.
# Produces: Table S4 (overall BC)
# Run     : sbatch scripts/05_replication/09d_tfgwas_overall_bc_discovery.sbatch   (from the repository root)
# ==============================================================================
#SBATCH --job-name=OverallBC_disc
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --array=1-37
#SBATCH --output=replication/logs/slurm_OverallBCdisc_%a_%j.out
#SBATCH --error=replication/logs/slurm_OverallBCdisc_%a_%j.err
set -euo pipefail
cd "${TNBC_ROOT:?set TNBC_ROOT}"   # repository root; the R script switches to replication/
mkdir -p replication/logs
export OVERALLBC_COHORT=discovery
echo "=== 09c_tfgwas_overall_bc.R DISCOVERY BATCH=${SLURM_ARRAY_TASK_ID} $(date) host=$(hostname) ==="
Rscript scripts/05_replication/09c_tfgwas_overall_bc.R ${SLURM_ARRAY_TASK_ID}
echo "=== exit=$? end=$(date) ==="
