#!/usr/bin/env bash
#SBATCH --job-name=finsirvauto
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000
#SBATCH --time=01:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4/slurm_logs/finsirvauto_%j.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/sirv4/_goal_opt/final_auto_v6
SAMPLE=SGNex_H9_directRNA_replicate3_run2

echo "[$(date --iso-8601=seconds)] START sirv auto routing"
python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain sirv \
  --profile sirv \
  --samples "$SAMPLE" \
  --guides p00,full \
  --variants baseline,overlap_locus \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE sirv auto routing"
