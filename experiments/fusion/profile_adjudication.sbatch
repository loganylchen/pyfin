#!/bin/bash
#SBATCH -J fusprof
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/fusprof_%j.out
#SBATCH -e /SSD/logan/dev/pyfin/logs/slurm/fusprof_%j.out
set -u
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
D=/SSD/logan/dev/pyfin/experiments/fusion/_k562_chr9_22
O=/SSD/logan/dev/pyfin/experiments/fusion/_prof
mkdir -p "$O"
# chr9+chr22 molecules (small N) against the WHOLE genome + full GENCODE, so
# each molecule pays exactly the candidate-window cost it pays in a real run.
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 -m cProfile -o "$O/adjudication.prof" \
  /SSD/logan/dev/pyfin/experiments/fusion/run_k562_fusion.py \
  --bam $D/reads.bam \
  --reference $B/refs/host/GRCh38.primary_assembly.genome.fa \
  --gtf $B/refs/host/gencode.v44.primary_assembly.annotation.gtf \
  --fastq $B/data/K562directRNAreplicate4run1/fastq/pass.fq.gz \
  --out "$O" --threads 8 --realign-workers 1
echo "rc=$? $(date +%T)"
