#!/bin/bash
#SBATCH -J k562_c922
#SBATCH -p mne4
#SBATCH -c 10
#SBATCH --mem=150G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/k562_c922_%j.out
# Real-library scale on two whole chromosomes (48,861 records): the same
# nomination density and gene index as a full run, but finishes, so the
# repaired code gets a scale check that is not the 98-record locus fixture.
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
O=/SSD/logan/dev/pyfin/experiments/fusion/_k562_chr9_22/out
rm -rf "$O"; mkdir -p "$O"
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 /SSD/logan/dev/pyfin/experiments/fusion/run_k562_fusion.py \
  --bam /SSD/logan/dev/pyfin/experiments/fusion/_k562_chr9_22/reads.bam \
  --reference /SSD/logan/dev/pyfin/experiments/fusion/_k562_bcrabl/fixture/chr9_chr22.fa \
  --gtf /SSD/logan/dev/pyfin/experiments/fusion/_k562_bcrabl/fixture/chr9_chr22.gtf \
  --fastq $B/data/K562directRNAreplicate4run1/fastq/pass.fq.gz \
  --out "$O" --threads 8 --realign-workers 8
echo "rc=$? $(date +%T)"
