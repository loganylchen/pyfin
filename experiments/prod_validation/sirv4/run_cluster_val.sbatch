#!/usr/bin/env bash
#SBATCH --job-name=clustval
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4/slurm_logs/clustval_%j.out
set -euo pipefail
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CODE=/SSD/logan/dev/pyfin                       # live code (uncommitted clustering flag)
S=SGNex_H9_directRNA_replicate2_run2
VAL="$ROOT/_cluster_val"
mkdir -p "$VAL"
for MODE in read_chains families; do
  for R in p00 full c_jitter20_10bp; do
    OUT="$VAL/$MODE/$R/pyfin_work"
    if [ -f "$OUT/assembly.gtf" ]; then echo "skip $MODE/$R"; continue; fi
    mkdir -p "$OUT"
    ANN="$ROOT/_ref/$R/annotation.gtf"
    NTX=$(grep -c $'\ttranscript\t' "$ANN" 2>/dev/null || echo 0)
    GTF_FLAG=(--gtf "$ANN"); [ "$NTX" -lt 2 ] && GTF_FLAG=()
    echo "[$(date)] RUN $MODE/$R ntx=$NTX"
    singularity exec --nv -B /SSD "$SIF" env PYTHONPATH="$CODE" \
      /usr/bin/python3.10 -m fin.cli \
        --bam "$ROOT/$S/stage/input.bam" --genome "$ROOT/$S/stage/genome.fa" \
        --fastq "$ROOT/$S/stage/reads.fq.gz" --signal "$ROOT/$S/stage/signal.blow5" \
        --signal-format slow5 --output-dir "$OUT" "${GTF_FLAG[@]}" \
        --gpu --quant-mode m2_em --clustering "$MODE" --mono-resolve-post-em
    echo "[$(date)] DONE $MODE/$R -> $(grep -c $'\ttranscript\t' "$OUT/assembly.gtf" 2>/dev/null || echo 0) tx"
  done
done
echo "[$(date)] ALL DONE"
