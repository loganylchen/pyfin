#!/usr/bin/env bash
#SBATCH --job-name=finrefitrep
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrefitrep_%j.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/final_profiles_v9_repeat_release/balanced

echo "[$(date --iso-8601=seconds)] START gencode_p00/baseline repeat"
PYFIN_REAL_THREADS=8 python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain real-drna \
  --profile real-drna \
  --samples gencode_p00 \
  --variants baseline \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE gencode_p00/baseline repeat"
