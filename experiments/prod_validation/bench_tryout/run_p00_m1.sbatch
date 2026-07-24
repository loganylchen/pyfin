#!/usr/bin/env bash
#SBATCH --job-name=pyfin_p00m1
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=24
#SBATCH --mem=80000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
B=$REPO/experiments/prod_validation/bench_tryout
OUT=$B/p00val/work_cluster_m1; rm -rf $OUT; mkdir -p $OUT
echo "[$(date)] p00 cluster margin=1 (keep wobble), recovery, M2 off, gates off"
singularity exec --nv -B /SSD -B /autofs "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/input.bam --genome /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/genome.fa \
    --fastq $B/p00val/reads_full.fq.gz --signal /autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5 --signal-format slow5 \
    --output-dir $OUT --gpu --quant-mode cluster --no-cluster-use-m2 \
    --cluster-m1-tie-margin 1 --min-polya5p-reads 0 --min-fulllen-fraction 0
cp $OUT/assembly.gtf $B/p00val/pyfin_cluster_m1.gtf
echo "[$(date)] DONE cluster_m1 EXIT=$? -> $(grep -c $'\ttranscript\t' $B/p00val/pyfin_cluster_m1.gtf 2>/dev/null) tx"
