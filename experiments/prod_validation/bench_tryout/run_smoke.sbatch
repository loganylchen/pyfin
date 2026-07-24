#!/usr/bin/env bash
#SBATCH --job-name=pyfin_smoke
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=48000
#SBATCH --time=00:40:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
B=$REPO/experiments/prod_validation/bench_tryout
OUT=$B/smoke/work_cluster; rm -rf $OUT; mkdir -p $OUT
echo "[$(date)] smoke cluster-mode start"
singularity exec --nv -B /SSD -B /autofs "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam $B/smoke/smoke.bam --genome /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/genome.fa \
    --fastq $B/smoke/reads.fq.gz --signal /autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5 --signal-format slow5 \
    --output-dir $OUT --gpu --quant-mode cluster
echo "[$(date)] EXIT=$? -> $(grep -c $'\ttranscript\t' $OUT/assembly.gtf 2>/dev/null) tx"
