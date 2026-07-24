#!/usr/bin/env bash
#SBATCH --job-name=pyfin-gc
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=32
#SBATCH --mem=250000
#SBATCH --time=2-00:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/cell_%j.out
set -euo pipefail
DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CELL=$(sed -n "$((CELL_IDX+1))p" "$DSROOT/cells.txt")
S=$(echo "$CELL"|awk '{print $1}'); R=$(echo "$CELL"|awk '{print $2}')
OUT="$DSROOT/$S/$R/pyfin_work"
[ -f "$DSROOT/$S/$R/pyfin.gtf" ] && { echo "skip $S/$R"; exit 0; }
mkdir -p "$OUT" "$DSROOT/slurm_logs"
echo "[$(date)] START $S/$R (human-scale)"
/usr/bin/time -v singularity exec --nv -B /SSD -B /autofs "$SIF" \
  fin --bam "$DSROOT/$S/stage/input.bam" --genome "$DSROOT/$S/stage/genome.fa" \
      --fastq "$DSROOT/$S/stage/reads.fq.gz" --signal "$DSROOT/$S/stage/signal.blow5" \
      --signal-format slow5 --output-dir "$OUT" \
      --gtf "$DSROOT/_ref/$R/annotation.gtf" --gpu --quant-mode m2_em
cp "$OUT/assembly.gtf" "$DSROOT/$S/$R/pyfin.gtf"
echo "[$(date)] DONE $S/$R -> $(grep -c $'\ttranscript\t' "$DSROOT/$S/$R/pyfin.gtf") tx"
