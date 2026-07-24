#!/usr/bin/env bash
#SBATCH --job-name=pyfin_align
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -euo pipefail
BIN=/SSD/logan/miniforge3/envs/nanofusion/bin
ST=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/stage_H9r2r2
[ -s "$ST/sorted.bam" ] && [ -s "$ST/sorted.bam.bai" ] && { echo "genome bam exists"; exit 0; }
echo "[$(date)] minimap2 dRNA splice: full reads -> GRCh38 genome"
"$BIN/minimap2" -t 16 -ax splice -uf -k14 --secondary=no \
  "$ST/genome.fa" "$ST/reads.fq.gz" 2> "$ST/minimap2.log" \
  | "$BIN/samtools" sort -@ 16 -m 2G -o "$ST/sorted.bam" -
"$BIN/samtools" index -@ 16 "$ST/sorted.bam"
echo "[$(date)] aligned. mapped reads:"
"$BIN/samtools" flagstat "$ST/sorted.bam" | sed -n '1,7p'
