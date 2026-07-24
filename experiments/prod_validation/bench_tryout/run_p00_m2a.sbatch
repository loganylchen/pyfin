#!/usr/bin/env bash
#SBATCH --job-name=pyfin_m2a
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=24
#SBATCH --mem=80000
#SBATCH --time=06:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
B=$REPO/experiments/prod_validation/bench_tryout
OUT=$B/p00val/work_m2a; rm -rf $OUT; mkdir -p $OUT
echo "[$(date)] p00 m2_em gates off sc250 + monoexon-fold + span-guard"
singularity exec --nv -B /SSD -B /autofs "$SIF" \
  env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/input.bam --genome /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/genome.fa \
    --fastq $B/p00val/reads_full.fq.gz --signal /autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5 --signal-format slow5 \
    --output-dir $OUT --gpu --quant-mode m2_em \
    --min-polya5p-reads 0 --min-fulllen-fraction 0 --chain-cluster-fold-monoexon
cp $OUT/assembly.gtf $B/p00val/pyfin_m2a.gtf
echo "[$(date)] DONE m2a EXIT=$? -> $(grep -c $'\ttranscript\t' $B/p00val/pyfin_m2a.gtf 2>/dev/null) tx"
