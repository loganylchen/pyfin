#!/usr/bin/env bash
#SBATCH --job-name=pyfin_pp
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
# nested add-back on top of smart(min-abundance 5, cluster-recheck+support-gate OFF):
#  pp1 = + polya5p 1
#  pp2 = + polya5p 1 + fulllen 0.1
#  pp3 = + polya5p 1 + fulllen 0.1 + soft-mass 2.0
NAMES=(pp1 pp2 pp3)
declare -a EXTRA
EXTRA[0]="--min-polya5p-reads 1 --min-fulllen-fraction 0 --max-soft-mass-ratio 0"
EXTRA[1]="--min-polya5p-reads 1 --min-fulllen-fraction 0.1 --max-soft-mass-ratio 0"
EXTRA[2]="--min-polya5p-reads 1 --min-fulllen-fraction 0.1 --max-soft-mass-ratio 2.0"
N=${NAMES[$SLURM_ARRAY_TASK_ID]}
E=${EXTRA[$SLURM_ARRAY_TASK_ID]}
OUT=$REPO/experiments/prod_validation/bench_tryout/prodfull_$N/p00
GTF=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep/_ref/p00/annotation.gtf
[ -s "$OUT/pyfin.gtf" ] && { echo "skip"; exit 0; }
for f in "$BAM" "$GENOME" "$FASTQ" "$SIGNAL" "$SIGNAL.idx" "$GTF"; do ls -L "$f" >/dev/null 2>&1 || { echo "MISSING $f"; exit 1; }; done
rm -rf "$OUT/work"; mkdir -p "$OUT/work"
echo "[$(date)] p00 $N : smart(min-ab5, no-cluster-recheck, no-support-gate) $E"
singularity exec --nv -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
  --bam "$BAM" --genome "$GENOME" --fastq "$FASTQ" \
  --signal "$SIGNAL" --signal-format slow5 --output-dir "$OUT/work" \
  --gpu --quant-mode m2_em --gtf "$GTF" \
  --min-abundance 5 --no-m2-cluster-recheck --no-m2-support-gate \
  $E
if [ -s "$OUT/work/assembly.gtf" ]; then
  cp "$OUT/work/assembly.gtf" "$OUT/pyfin.gtf"
  echo "[$(date)] DONE $N -> $(grep -c $'\ttranscript\t' "$OUT/pyfin.gtf") tx"
else echo "[$(date)] FAIL $N"; exit 1; fi
