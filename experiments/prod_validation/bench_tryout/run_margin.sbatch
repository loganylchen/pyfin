#!/usr/bin/env bash
#SBATCH --job-name=pyfin_margin
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=16
#SBATCH --mem=48000
#SBATCH --time=03:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
cd $REPO/experiments/prod_validation/bench_tryout
singularity exec --nv -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
  env PYTHONPATH=$REPO LIMIT=5000 USE_GPU=1 /usr/bin/python3.10 margin_dist.py
