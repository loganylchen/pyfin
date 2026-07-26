#!/usr/bin/env bash
#SBATCH --job-name=flipval
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=24000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4/slurm_logs/flipval_%j.out
set -euo pipefail
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4
SIF=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
CODE=/SSD/logan/dev/pyfin
S=SGNex_H9_directRNA_replicate2_run2
VAL="$ROOT/_flip_val"; mkdir -p "$VAL"
# CONFIG: name | extra flags
run(){ local NAME="$1"; shift
  for R in p00 full; do
    OUT="$VAL/$NAME/$R/pyfin_work"
    [ -f "$OUT/assembly.gtf" ] && { echo "skip $NAME/$R"; continue; }
    mkdir -p "$OUT"
    ANN="$ROOT/_ref/$R/annotation.gtf"; NTX=$(grep -c $'\ttranscript\t' "$ANN" 2>/dev/null||echo 0)
    GTF=(--gtf "$ANN"); [ "$NTX" -lt 2 ] && GTF=()
    echo "[$(date)] RUN $NAME/$R"
    singularity exec --nv -B /SSD "$SIF" env PYTHONPATH="$CODE" /usr/bin/python3.10 -m fin.cli \
      --bam "$ROOT/$S/stage/input.bam" --genome "$ROOT/$S/stage/genome.fa" \
      --fastq "$ROOT/$S/stage/reads.fq.gz" --signal "$ROOT/$S/stage/signal.blow5" \
      --signal-format slow5 --output-dir "$OUT" "${GTF[@]}" --gpu --quant-mode m2_em "$@"
    echo "[$(date)] DONE $NAME/$R -> $(grep -c $'\ttranscript\t' "$OUT/assembly.gtf" 2>/dev/null||echo 0)"
  done
}
run baseline_prod   --clustering read_chains
run clean_families  --clustering families --mono-resolve-post-em
echo "[$(date)] ALL DONE"
