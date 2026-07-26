#!/usr/bin/env bash
#SBATCH --job-name=gcflip
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=24
#SBATCH --mem=90000
#SBATCH --time=10:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/slurm_logs/gcflip_%j.out
set -euo pipefail
DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CODE=/SSD/logan/dev/pyfin
STAGE="$DSROOT/_p00val/stage"
run(){ local NAME="$1"; shift
  local OUT="$DSROOT/_flipval/$NAME/pyfin_work"
  [ -f "$OUT/assembly.gtf" ] && { echo "skip $NAME"; return; }
  mkdir -p "$OUT"; echo "[$(date)] RUN $NAME"
  singularity exec --nv -B /SSD "$SIF" env PYTHONPATH="$CODE" /usr/bin/python3.10 -m fin.cli \
    --bam "$STAGE/input.bam" --genome "$STAGE/genome.fa" --fastq "$STAGE/reads.fq.gz" \
    --signal "$STAGE/signal.blow5" --signal-format slow5 --output-dir "$OUT" \
    --gpu --quant-mode m2_em "$@"
  echo "[$(date)] DONE $NAME -> $(grep -c $'\ttranscript\t' "$OUT/assembly.gtf" 2>/dev/null||echo 0)"
}
run baseline_prod  --clustering read_chains
run clean_families --clustering families --mono-resolve-post-em
echo "[$(date)] ALL DONE"
