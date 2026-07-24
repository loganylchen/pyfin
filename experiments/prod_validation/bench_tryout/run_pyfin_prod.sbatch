#!/usr/bin/env bash
#SBATCH --job-name=pyfin_bench
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%A_%a.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
DS=$REPO/experiments/prod_validation
S=SGNex_H9_directRNA_replicate2_run2
RATIOS=(full p99 p90 p50 p10 p00 c_skip20 c_jitter20_10bp c_spurious20 c_merge20 c_flip20 c_ir20)
R=${RATIOS[$SLURM_ARRAY_TASK_ID]}
OUT=$DS/bench_tryout/prod/$S/$R
[ -s "$OUT/pyfin.gtf" ] && { echo "skip $R"; exit 0; }
rm -rf "$OUT/work"; mkdir -p "$OUT/work"
st=$DS/gencode/$S/stage
ANN=$(readlink -f "$DS/gencode/_ref/$R/annotation.gtf")
echo "[$(date)] $S/$R  ann=$ANN  (production defaults, Lever2 on)"
# PURE production defaults: no lever flags (config default novel_junction_min_reads=2)
singularity exec --nv -B /SSD -B /autofs/mnemosyne3_SSD "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
  --bam "$st/input.bam" --genome "$st/genome.fa" --fastq "$st/reads.fq.gz" \
  --signal "$st/signal.blow5" --signal-format slow5 --output-dir "$OUT/work" \
  --gpu --quant-mode m2_em --gtf "$ANN"
if [ -s "$OUT/work/assembly.gtf" ]; then
  cp "$OUT/work/assembly.gtf" "$OUT/pyfin.gtf"
  echo "[$(date)] DONE $R -> $(grep -c $'\ttranscript\t' "$OUT/pyfin.gtf") tx"
else
  echo "[$(date)] FAIL $R"; exit 1
fi
