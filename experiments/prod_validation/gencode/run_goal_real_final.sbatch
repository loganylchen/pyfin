#!/usr/bin/env bash
#SBATCH --job-name=finrealfinal
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=02:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrealfinal_%j.out

set -euo pipefail
REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/final_profile_strict_v2
cd "$REPO"
echo "[$(date --iso-8601=seconds)] START live named real-drna profile"
PYFIN_REAL_THREADS=8 python3 experiments/prod_validation/profile_sweep.py \
  --domain real-drna \
  --variants baseline \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE live named real-drna profile"
