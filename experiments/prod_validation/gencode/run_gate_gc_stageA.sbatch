#!/usr/bin/env bash
#SBATCH --job-name=gcgateA
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=24
#SBATCH --mem=90000
#SBATCH --time=20:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/gcgateA_%j.out
set -euo pipefail
DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CODE=/SSD/logan/dev/pyfin
STAGE="$DSROOT/_p00val/stage"
JD="--junction-dominance-filter"
declare -A CFG=(
  [J0]=""
  [JS]="$JD --junction-dominance-min-reads 2 --junction-dominance-window-bp 2 --junction-dominance-tol-bp 2"
  [JD1]="$JD --junction-dominance-min-reads 1 --junction-dominance-window-bp 20 --junction-dominance-tol-bp 2"
  [J10]="$JD --junction-dominance-min-reads 2 --junction-dominance-window-bp 10 --junction-dominance-tol-bp 2"
  [J20]="$JD --junction-dominance-min-reads 2 --junction-dominance-window-bp 20 --junction-dominance-tol-bp 2"
  [J30]="$JD --junction-dominance-min-reads 2 --junction-dominance-window-bp 30 --junction-dominance-tol-bp 2"
  [COFF]="--no-canonical-gate"
)
for NAME in J0 JS JD1 J10 J20 J30 COFF; do
  OUT="$DSROOT/_gate_gc/$NAME/pyfin_work"
  [ -f "$OUT/assembly.gtf" ] && { echo "skip $NAME"; continue; }
  mkdir -p "$OUT"; echo "[$(date)] RUN $NAME  flags=[${CFG[$NAME]}]"
  singularity exec --nv -B /SSD "$SIF" env PYTHONPATH="$CODE" /usr/bin/python3.10 -m fin.cli \
    --bam "$STAGE/input.bam" --genome "$STAGE/genome.fa" --fastq "$STAGE/reads.fq.gz" \
    --signal "$STAGE/signal.blow5" --signal-format slow5 --output-dir "$OUT" \
    --gpu --quant-mode m2_em ${CFG[$NAME]}
  echo "[$(date)] DONE $NAME -> $(grep -c $'\ttranscript\t' "$OUT/assembly.gtf" 2>/dev/null||echo 0)"
done
echo "[$(date)] STAGE A DONE"
