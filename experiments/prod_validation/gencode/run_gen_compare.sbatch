#!/usr/bin/env bash
#SBATCH --job-name=gencmp
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=48000
#SBATCH --time=04:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/gencmp_%j.out
set -euo pipefail
DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CODE=/SSD/logan/dev/pyfin
STAGE="$DSROOT/_p00val/stage"
OUT="$DSROOT/_p00val/_gen"; mkdir -p "$OUT"
echo "[$(date)] gen-compare START"
singularity exec -B /SSD "$SIF" env PYTHONPATH="$CODE" /usr/bin/python3.10 -u \
  "$CODE/experiments/wobble_heya8/measure_gen_compare.py" \
  "$STAGE/input.bam" "$STAGE/genome.fa" "$OUT/read_chains.gtf" "$OUT/families.gtf"
echo "[$(date)] gen-compare DONE"
