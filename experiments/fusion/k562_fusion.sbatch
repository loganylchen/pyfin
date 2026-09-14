#!/bin/bash
#SBATCH -J k562_fusion
#SBATCH -p mne4
#SBATCH -c 10
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/k562_fusion_%j.out
# First look at what the owner's independent fusion analysis produces on a
# REAL sample with a known positive control (K562, BCR-ABL1 chr22->chr9,
# 7 SA-tagged chimeric reads in this BAM). Sequence-only (no --signal) first.
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
S=SGNex_K562_directRNA_replicate4_run1
O=/SSD/logan/dev/pyfin/experiments/fusion/_k562_rep4_seq
mkdir -p "$O"
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 /SSD/logan/dev/pyfin/experiments/fusion/run_k562_fusion.py \
  --bam $B/results/gencode_full_sweep_ratio_grid_v2/$S/align/$S.sorted.bam \
  --reference $B/refs/host/GRCh38.primary_assembly.genome.fa \
  --gtf $B/refs/host/gencode.v44.primary_assembly.annotation.gtf \
  --fastq $B/data/K562directRNAreplicate4run1/fastq/pass.fq.gz \
  --out "$O" --threads 8 --realign-workers 8
echo "rc=$? $(date +%T)"
