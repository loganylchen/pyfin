#!/usr/bin/env bash
#SBATCH --job-name=finrealprec
#SBATCH --partition=mne4
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrealprec_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/real_precision
VARIANTS=(mono3 mono5 mono5_r2 mono5_r2_ab15)
VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}

echo "[$(date --iso-8601=seconds)] START $VARIANT"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $VARIANT"
