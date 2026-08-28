#!/usr/bin/env bash
#SBATCH --job-name=finrefit
#SBATCH --partition=mne4
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrefit_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/final_profiles_v9_release/balanced
VARIANTS=(baseline refit_off)
VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}

echo "[$(date --iso-8601=seconds)] START gencode_p00/$VARIANT"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --profile real-drna \
  --samples gencode_p00 \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE gencode_p00/$VARIANT"
