#!/bin/bash
#SBATCH --job-name=stepA_test
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:25:00
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
cd $REPO
echo "=== pytest unit suite (Step A tree, --nv GPU) ==="
singularity exec --nv -B /SSD "$SIF" env PYTHONPATH=$REPO /usr/bin/python3.10 -m pytest \
  tests/unit -q --no-header 2>&1 | tail -45
echo "PYTEST_EXIT=${PIPESTATUS[0]}"
