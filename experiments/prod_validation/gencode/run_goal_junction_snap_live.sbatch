#!/usr/bin/env bash
#SBATCH --job-name=finjsnap
#SBATCH --partition=mne4
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=04:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finjsnap_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/junction_snap_live
TASKS=(
  "gencode_p00 mono5_snap6"
  "gencode_p00 mono5_r2_snap6"
  "SGNex_H9_directRNA_replicate3_run1 mono5_snap6"
  "SGNex_H9_directRNA_replicate3_run1 mono5_r2_snap6"
)
read -r SAMPLE VARIANT <<< "${TASKS[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date --iso-8601=seconds)] START $SAMPLE/$VARIANT"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --samples "$SAMPLE" \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $SAMPLE/$VARIANT"
