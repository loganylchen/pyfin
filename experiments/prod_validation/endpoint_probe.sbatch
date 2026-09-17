#!/bin/bash
#SBATCH -J epprobe
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=150G
#SBATCH -t 1:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/epprobe_%j.out
# Codex's hypothesis, timed: _endpoint_context() refits the endpoint mixtures
# over the whole library on EVERY interval. One span, absolute seconds,
# --threads 1 so the work is in-process, real library, prepass cache warm.
set -u
cd /SSD/logan/dev/pyfin
P=/SSD/logan/dev/pyfin/experiments/prod_validation/_annoprobe
S=SGNex_H9_directRNA_replicate2_run2
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
G=$B/refs/host/GRCh38.primary_assembly.genome.fa
ANNO=$B/results/gencode_full_sweep_ratio_grid_v2/_ref/$S/full/annotation.gtf
ST=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage
W=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/_sample_evidence_cache
timeout 1500 singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin PYFIN_SAMPLE_EVIDENCE_CACHE=$W /usr/bin/python3.10 -m fin.cli \
  --bam $ST/input.bam --genome $G --fastq $ST/reads.fq.gz --signal $ST/signal.blow5 \
  --output-dir $P/ep.out --profile real-drna --no-gpu --threads 1 \
  --gtf $ANNO --m2-support-gate-attest-reads 1 > $P/ep.log 2>&1
echo "STOPPED $(date +%T) intervals=$(grep -c 'Completed interval' $P/ep.log) probes=$(grep -c 'PROBE endpoint_context' $P/ep.log)"
grep "PROBE endpoint_context" $P/ep.log | tail -3
