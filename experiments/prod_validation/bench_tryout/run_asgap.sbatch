#!/usr/bin/env bash
#SBATCH --job-name=pyfin_asgap
#SBATCH --partition=mne4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --time=01:00:00
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
cd $REPO/experiments/prod_validation/bench_tryout
singularity exec -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
  env PYTHONPATH=$REPO LIMIT=3000 /usr/bin/python3.10 measure_as_gap.py
