#!/usr/bin/env bash
#SBATCH --job-name=pyfin_unfilt
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=48000
#SBATCH --time=06:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
OUT=$REPO/experiments/prod_validation/bench_tryout/prodfull_unfilt
BAM=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep/SGNex_H9_directRNA_replicate2_run2/align/SGNex_H9_directRNA_replicate2_run2.sorted.bam
GENOME=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/GRCh38.primary_assembly.genome.fa
FASTQ=/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/fastq/pass.fq.gz; SIGNAL=/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5
GTF=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep/_ref/full/annotation.gtf
[ -s "$OUT/work/assembly.gtf" ] && { echo skip; exit 0; }
rm -rf "$OUT/work"; mkdir -p "$OUT/work"
for f in "$BAM" "$GENOME" "$FASTQ" "$SIGNAL" "$SIGNAL.idx" "$GTF"; do ls -L "$f" >/dev/null 2>&1 || { echo "MISSING $f"; exit 1; }; done
echo "[$(date)] full + --write-unfiltered-scores (production defaults)"
singularity exec --nv -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
  --bam "$BAM" --genome "$GENOME" --fastq "$FASTQ" --signal "$SIGNAL" --signal-format slow5 \
  --output-dir "$OUT/work" --gpu --quant-mode m2_em --write-unfiltered-scores --gtf "$GTF"
echo "[$(date)] done. final tx=$(grep -c $'\ttranscript\t' "$OUT/work/assembly.gtf" 2>/dev/null); unfiltered rows=$(wc -l < "$OUT/work/scores.unfiltered.tsv" 2>/dev/null)"
