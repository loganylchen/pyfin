#!/usr/bin/env bash
#SBATCH --job-name=sirvrefit
#SBATCH --partition=mne4
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --time=01:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4/slurm_logs/sirvrefit_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/sirv4/_goal_opt/final_auto_v7_release
SAMPLE=SGNex_H9_directRNA_replicate3_run2
VARIANTS=(baseline refit_off)
VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}

echo "[$(date --iso-8601=seconds)] START $SAMPLE/p00/$VARIANT"
python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain sirv \
  --profile sirv \
  --samples "$SAMPLE" \
  --guides p00 \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $SAMPLE/p00/$VARIANT"
