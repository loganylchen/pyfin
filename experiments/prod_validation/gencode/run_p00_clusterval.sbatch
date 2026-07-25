#!/usr/bin/env bash
#SBATCH --job-name=gcp00val
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=24
#SBATCH --mem=90000
#SBATCH --time=12:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/p00val_%j.out
set -euo pipefail
NFS=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
SAMPLE=H9directRNAreplicate2run2
DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CODE=/SSD/logan/dev/pyfin
STAGE="$DSROOT/_p00val/stage"
mkdir -p "$STAGE" "$DSROOT/slurm_logs"

# --- stage inputs local (idempotent) ---
BAM_SRC=/SSD/logan/dev/pyfin/experiments/wobble_heya8/cf_measure/genome.bam   # already local, full genome aln
FA_SRC="$NFS/refs/host/GRCh38.primary_assembly.genome.fa"
FQ_SRC="$NFS/data/$SAMPLE/fastq/pass.fq.gz"
SIG_SRC="$NFS/data/$SAMPLE/blow5/nanopore.blow5"
ln -sf "$BAM_SRC" "$STAGE/input.bam"; ln -sf "$BAM_SRC.bai" "$STAGE/input.bam.bai"
stage_cp(){ local src="$1" dst="$2"; if [ ! -s "$dst" ] || [ "$(stat -c%s "$src")" != "$(stat -c%s "$dst")" ]; then echo "[$(date)] copy $(basename "$dst")"; cp -f "$src" "$dst"; else echo "have $(basename "$dst")"; fi; }
stage_cp "$FA_SRC" "$STAGE/genome.fa"
stage_cp "$FA_SRC.fai" "$STAGE/genome.fa.fai"
stage_cp "$FQ_SRC" "$STAGE/reads.fq.gz"
stage_cp "$SIG_SRC" "$STAGE/signal.blow5"
stage_cp "$SIG_SRC.idx" "$STAGE/signal.blow5.idx"
echo "[$(date)] staging done"

# --- run both clustering modes, p00 (de novo, no --gtf) ---
for MODE in read_chains families; do
  OUT="$DSROOT/_p00val/$MODE/pyfin_work"
  if [ -f "$OUT/assembly.gtf" ]; then echo "skip $MODE"; continue; fi
  mkdir -p "$OUT"
  echo "[$(date)] RUN gencode p00 $MODE"
  singularity exec --nv -B /SSD "$SIF" env PYTHONPATH="$CODE" \
    /usr/bin/python3.10 -m fin.cli \
      --bam "$STAGE/input.bam" --genome "$STAGE/genome.fa" \
      --fastq "$STAGE/reads.fq.gz" --signal "$STAGE/signal.blow5" \
      --signal-format slow5 --output-dir "$OUT" \
      --gpu --quant-mode m2_em --clustering "$MODE" --mono-resolve-post-em
  echo "[$(date)] DONE $MODE -> $(grep -c $'\ttranscript\t' "$OUT/assembly.gtf" 2>/dev/null || echo 0) tx"
done
echo "[$(date)] ALL DONE"
