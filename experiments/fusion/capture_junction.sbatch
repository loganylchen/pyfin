#!/bin/bash
#SBATCH -J juncap
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/juncap_%j.out
#SBATCH -e /SSD/logan/dev/pyfin/logs/slurm/juncap_%j.out
set -u
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
D=/SSD/logan/dev/pyfin/experiments/fusion/_k562_chr9_22
O=/SSD/logan/dev/pyfin/experiments/fusion/_juncap
mkdir -p "$O"
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 \
  /SSD/logan/dev/pyfin/experiments/fusion/capture_junction_calls.py \
  --bam $D/reads.bam \
  --reference $B/refs/host/GRCh38.primary_assembly.genome.fa \
  --gtf $B/refs/host/gencode.v44.primary_assembly.annotation.gtf \
  --fastq $B/data/K562directRNAreplicate4run1/fastq/pass.fq.gz \
  --out "$O/run" --corpus "$O/junction_calls.pkl" --threads 8 --realign-workers 8
echo "rc=$? $(date +%T)"
