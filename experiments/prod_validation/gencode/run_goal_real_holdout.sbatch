#!/usr/bin/env bash
#SBATCH --job-name=finrealhold
#SBATCH --partition=mne4
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=04:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrealhold_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/real_holdout_r3r1_v2
SAMPLE=SGNex_H9_directRNA_replicate3_run1
VARIANTS=(baseline mono5 mono5_r2)
VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}

echo "[$(date --iso-8601=seconds)] START $SAMPLE/$VARIANT"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --samples "$SAMPLE" \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $SAMPLE/$VARIANT"
