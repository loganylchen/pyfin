#!/usr/bin/env bash
#SBATCH --job-name=finrefitho
#SBATCH --partition=mne4
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=04:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrefitho_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
TASKS=(
  "real-drna balanced SGNex_H9_directRNA_replicate3_run1"
  "real-drna-precision precision gencode_p00"
)
read -r PROFILE LABEL SAMPLE <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/final_profiles_v9_release/$LABEL

echo "[$(date --iso-8601=seconds)] START $PROFILE/$SAMPLE/baseline"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --profile "$PROFILE" \
  --samples "$SAMPLE" \
  --variants baseline \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $PROFILE/$SAMPLE/baseline"
