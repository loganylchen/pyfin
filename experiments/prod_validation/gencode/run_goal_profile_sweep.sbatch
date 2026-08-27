#!/usr/bin/env bash
#SBATCH --job-name=finrealopt
#SBATCH --partition=mne4
#SBATCH --array=0-14%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=200000
#SBATCH --time=06:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/finrealopt_%A_%a.out

set -euo pipefail

REPO=/SSD/logan/dev/pyfin
SNAP=$REPO/experiments/prod_validation/gencode/_goal_opt/profile_code_20260826
OUT=$REPO/experiments/prod_validation/gencode/_goal_opt/profile_sweep_retry
SCRIPT=$REPO/experiments/prod_validation/profile_sweep.py
VARIANTS=(
  baseline
  m2_sum
  m2_off
  abundance_0
  abundance_2
  abundance_3
  finalize_sirv
  fulllen_on
  polya_on
  junction_1
  containment_off
  recheck_off
  canonical_off
  softmass_off
  isoform_fraction_off
)

VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
echo "[$(date --iso-8601=seconds)] START $VARIANT snapshot=$SNAP"
PYFIN_CODE_ROOT=$SNAP PYFIN_REAL_THREADS=8 python3 "$SCRIPT" \
  --domain real-drna \
  --variants "$VARIANT" \
  --output-root "$OUT" \
  --run-only
echo "[$(date --iso-8601=seconds)] DONE $VARIANT"
