#!/bin/bash
# ==============================================================================
# 12_tf_family_models.sbatch
# SLURM wrapper for 12_tf_family_models.R.
# Produces: Table S8
# Run     : sbatch scripts/05_replication/12_tf_family_models.sbatch   (from the repository root)
# ==============================================================================
#SBATCH --job-name=grp_tfgwas
##SBATCH --partition=<your_partition>   # site-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --output=replication/logs/slurm_grptfgwas_%j.out
#SBATCH --error=replication/logs/slurm_grptfgwas_%j.err
#
# Grouped TF-GWAS (motif 727 / TF family 94), observed statistics only.
#   OBS_COHORT=discovery  OBS_CSV=results/tfgwas/variant_peak_motif_subtypes_ALL.csv  OBS_POSCOL=Pos
#   OBS_COHORT=replication OBS_CSV=replication/tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz OBS_POSCOL=Pos38
# POSCOL must be the build the peak/motif overlap was computed in.
cd "${TNBC_ROOT:?set TNBC_ROOT}"
mkdir -p replication/logs
mkdir -p logs replication/tfgwas_grouped
LOG=replication/logs/12_tf_family_models_${OBS_COHORT:-discovery}.log
{
  echo "=== 12_tf_family_models.R  $(date)  host=$(hostname)  job=$SLURM_JOB_ID ==="
  echo "cohort=${OBS_COHORT} csv=${OBS_CSV} poscol=${OBS_POSCOL}"
  Rscript scripts/05_replication/12_tf_family_models.R
  echo "=== exit=$? $(date) ==="
} 2>&1 | tee -a "$LOG"
