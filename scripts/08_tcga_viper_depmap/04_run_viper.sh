#!/bin/bash
# 04_run_viper.sh -- run the msVIPER analysis reported in Fig 4 and Table S10.
# Regulon minimum size 5 (VIPER_MINSIZE=5); intermediate results are cached under
# results/tf_activity/viper_vst/steps/ and reused on re-runs.
# Run from the repository root:  bash scripts/08_tcga_viper_depmap/04_run_viper.sh
set -euo pipefail
export VIPER_MINSIZE=5
mkdir -p results/tf_activity/logs
Rscript scripts/08_tcga_viper_depmap/04_viper_tf_activity.R > results/tf_activity/logs/viper_min5.log 2>&1
echo "done -- log: results/tf_activity/logs/viper_min5.log"
