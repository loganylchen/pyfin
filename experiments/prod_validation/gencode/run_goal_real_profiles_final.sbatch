#!/usr/bin/env bash
#SBATCH --job-name=finprofiles
#SBATCH --partition=mne4
#SBATCH --array=0-2%1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finprofiles_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
TASKS=(
  "real-drna balanced baseline"
  "real-drna-precision precision baseline"
  "real-drna-precision precision overlap_locus"
)
read -r PROFILE LABEL VARIANT <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/final_profiles_v8/$LABEL

echo "[$(date --iso-8601=seconds)] START $PROFILE/$LABEL/$VARIANT"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --profile "$PROFILE" \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $PROFILE/$LABEL/$VARIANT"
