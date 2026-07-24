#!/usr/bin/env bash
#SBATCH --job-name=pyfin1
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000
#SBATCH --time=02:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4/slurm_logs/cell_%j.out
set -euo pipefail
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CELL=$(sed -n "$((CELL_IDX+1))p" "$ROOT/cells.txt")
S=$(echo "$CELL"|awk '{print $1}'); R=$(echo "$CELL"|awk '{print $2}')
OUT="$ROOT/$S/$R/pyfin_work"
[ -f "$ROOT/$S/$R/pyfin.gtf" ] && { echo "skip $S/$R"; exit 0; }
mkdir -p "$OUT"; ANN="$ROOT/_ref/$R/annotation.gtf"
NTX=$(grep -c $'\ttranscript\t' "$ANN" 2>/dev/null || echo 0)
GTF_FLAG=(--gtf "$ANN"); [ "$NTX" -lt 2 ] && GTF_FLAG=()
echo "[$(date)] $S/$R ntx=$NTX"
singularity exec --nv -B /SSD "$SIF" \
  fin --bam "$ROOT/$S/stage/input.bam" --genome "$ROOT/$S/stage/genome.fa" \
      --fastq "$ROOT/$S/stage/reads.fq.gz" --signal "$ROOT/$S/stage/signal.blow5" \
      --signal-format slow5 --output-dir "$OUT" "${GTF_FLAG[@]}" --gpu --quant-mode m2_em
cp "$OUT/assembly.gtf" "$ROOT/$S/$R/pyfin.gtf"
echo "[$(date)] DONE $S/$R -> $(grep -c $'\ttranscript\t' "$ROOT/$S/$R/pyfin.gtf") tx"
