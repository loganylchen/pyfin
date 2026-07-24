#!/bin/bash
#SBATCH --job-name=trace_reads
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=00:30:00
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
W=$REPO/experiments/wobble_heya8/cf_measure
TRUTH=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/gencode.v44.primary_assembly.annotation.gtf
echo "[$(date)] trace missed FL reads by name"
singularity exec -B /SSD -B /autofs "$SIF" env PYTHONPATH=$REPO /usr/bin/python3.10 -u \
  $REPO/experiments/wobble_heya8/trace_missed_reads.py "$W/transcriptome.bam" "$W/genome.bam" "$TRUTH" "$W/gc_cf.cf.gtf.refmap"
echo "[$(date)] DONE"
