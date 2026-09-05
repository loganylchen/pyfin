#!/bin/bash
#SBATCH -J pyfin_tests
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/pyfin_tests_%j.out
cd /SSD/logan/dev/pyfin
singularity exec --nv -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin python3 -m pytest -q -p no:cacheprovider tests/unit tests/integration 2>&1 | tail -40
