#!/usr/bin/env bash
#SBATCH --job-name=pyfin-gcdn
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=24
#SBATCH --mem=80000
#SBATCH --time=2-00:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/dn_%j.out
set -euo pipefail
DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
S=$(sed -n "$((CELL_IDX+1))p" "$DSROOT/cells_nogtf.txt")
OUT="$DSROOT/$S/nogtf/pyfin_work"
[ -f "$DSROOT/$S/nogtf/pyfin.gtf" ] && { echo "skip $S"; exit 0; }
mkdir -p "$OUT"
echo "[$(date)] START de novo (no --gtf) $S"
singularity exec --nv -B /SSD -B /autofs "$SIF" \
  fin --bam "$DSROOT/$S/stage/input.bam" --genome "$DSROOT/$S/stage/genome.fa" \
      --fastq "$DSROOT/$S/stage/reads.fq.gz" --signal "$DSROOT/$S/stage/signal.blow5" \
      --signal-format slow5 --output-dir "$OUT" --gpu --quant-mode m2_em
cp "$OUT/assembly.gtf" "$DSROOT/$S/nogtf/pyfin.gtf"
echo "[$(date)] DONE $S -> $(grep -c $'\ttranscript\t' "$DSROOT/$S/nogtf/pyfin.gtf") tx"
