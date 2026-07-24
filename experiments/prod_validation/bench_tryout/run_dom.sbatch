#!/usr/bin/env bash
#SBATCH --job-name=pyfin_dom
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=48000
#SBATCH --time=06:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%A_%a.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
BAM=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep/SGNex_H9_directRNA_replicate2_run2/align/SGNex_H9_directRNA_replicate2_run2.sorted.bam
GENOME=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/GRCh38.primary_assembly.genome.fa
FASTQ=/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/fastq/pass.fq.gz
SIGNAL=/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5
OUTROOT=$REPO/experiments/prod_validation/bench_tryout/prodfull_dom
RATIOS=(full p99 p90 p50 p10 p00 c_skip20 c_jitter20_10bp c_spurious20 c_merge20 c_flip20 c_ir20)
R=${RATIOS[$SLURM_ARRAY_TASK_ID]}
GTF=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep/_ref/$R/annotation.gtf
OUT=$OUTROOT/$R
[ -s "$OUT/pyfin.gtf" ] && { echo "skip $R"; exit 0; }
for f in "$BAM" "$GENOME" "$FASTQ" "$SIGNAL" "$SIGNAL.idx" "$GTF"; do ls -L "$f" >/dev/null 2>&1 || { echo "MISSING $f"; exit 1; }; done
rm -rf "$OUT/work"; mkdir -p "$OUT/work"
echo "[$(date)] $R  (made genome bam + FULL fastq/blow5 originals, bind mne3+NAS25, production incl Lever2)"
singularity exec --nv -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
  --bam "$BAM" --genome "$GENOME" --fastq "$FASTQ" \
  --signal "$SIGNAL" --signal-format slow5 --output-dir "$OUT/work" \
  --gpu --quant-mode m2_em --junction-dominance-filter --gtf "$GTF"
if [ -s "$OUT/work/assembly.gtf" ]; then
  cp "$OUT/work/assembly.gtf" "$OUT/pyfin.gtf"
  echo "[$(date)] DONE $R -> $(grep -c $'\ttranscript\t' "$OUT/pyfin.gtf") tx"
else echo "[$(date)] FAIL $R"; exit 1; fi
