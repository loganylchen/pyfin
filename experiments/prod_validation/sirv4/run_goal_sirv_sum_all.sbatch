#!/usr/bin/env bash
#SBATCH --job-name=finsirvsum
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000
#SBATCH --time=01:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4/slurm_logs/finsirvsum_%j.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
OUT=$REPO/experiments/prod_validation/sirv4/_goal_opt/sum_guide_split_true_p00
SAMPLES=SGNex_H9_directRNA_replicate2_run1,SGNex_H9_directRNA_replicate2_run2,SGNex_H9_directRNA_replicate3_run1,SGNex_H9_directRNA_replicate3_run2,SGNex_H9_directRNA_replicate4_run1,SGNex_H9_directRNA_replicate4_run2

echo "[$(date --iso-8601=seconds)] START six-sample mean-vs-sum validation"
python3 "$REPO/experiments/prod_validation/profile_sweep.py" \
  --domain sirv \
  --samples "$SAMPLES" \
  --guides p00,full \
  --variants best_mean,best_sum_m1_f4 \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE six-sample mean-vs-sum validation"
