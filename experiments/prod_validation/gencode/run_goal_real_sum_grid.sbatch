#!/usr/bin/env bash
#SBATCH --job-name=finsumgrid
#SBATCH --partition=mne4
#SBATCH --array=0-8%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finsumgrid_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
SNAP=$REPO/experiments/prod_validation/gencode/_goal_opt/profile_code_20260826
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/real_sum_grid
SCRIPT=$REPO/experiments/prod_validation/profile_sweep.py
VARIANTS=(
  sum_m1_f4 sum_m1_f6 sum_m1_f8
  sum_m2_f4 sum_m2_f6 sum_m2_f8
  sum_m3_f4 sum_m3_f6 sum_m3_f8
)

VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
echo "[$(date --iso-8601=seconds)] START $VARIANT snapshot=$SNAP"
PYFIN_CODE_ROOT=$SNAP PYFIN_REAL_THREADS=8 python3 "$SCRIPT" \
  --domain real-drna \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $VARIANT"
