#!/bin/bash
#SBATCH --job-name=cf_measure
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:45:00
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
BAM=$REPO/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/input.bam
W=$REPO/experiments/wobble_heya8/cf_measure
echo "[$(date)] start cluster_families emission"
singularity exec -B /SSD -B /autofs "$SIF" env PYTHONPATH=$REPO /usr/bin/python3.10 \
  $REPO/experiments/wobble_heya8/measure_cf.py "$BAM" "$W/cf.gtf" --mono
echo "[$(date)] DONE -> $(grep -c $'\ttranscript\t' $W/cf.gtf) transcripts"
