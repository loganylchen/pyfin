#!/bin/bash
#SBATCH -J k562_bcrabl
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 6:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/k562_bcrabl_%j.out
# Positive control on the real BCR-ABL1 locus only (98 records, 17 with SA),
# against the FULL human reference so competing mappings are still searched.
# Runs the current tree (transcript layer included).
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
O=/SSD/logan/dev/pyfin/experiments/fusion/_k562_bcrabl/out_cov${COV:-0.9}
rm -rf "$O"; mkdir -p "$O"
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 /SSD/logan/dev/pyfin/experiments/fusion/run_k562_fusion.py \
  --bam /SSD/logan/dev/pyfin/experiments/fusion/_k562_bcrabl/bcr_abl1.bam \
  --reference $B/refs/host/GRCh38.primary_assembly.genome.fa \
  --gtf $B/refs/host/gencode.v44.primary_assembly.annotation.gtf \
  --fastq $B/data/K562directRNAreplicate4run1/fastq/pass.fq.gz \
  --out "$O" --threads 8 --min-clip 250 --min-read-coverage ${COV:-0.9}
echo "rc=$? $(date +%T)"
