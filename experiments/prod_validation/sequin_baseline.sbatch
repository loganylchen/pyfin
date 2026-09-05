#!/bin/bash
#SBATCH -J sqbase
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 3:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/sqbase_%j.out
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
echo "=== SEQUIN default (real-drna, no extra flags), CURRENT code $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --tag _sq_base --out $SQ/base_real --tsv $SQ/base_real.tsv 2>&1 | tail -2
echo "=== SEQUIN sirv profile, CURRENT code $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile sirv --tag _sq_base_sirv --out $SQ/base_sirv --tsv $SQ/base_sirv.tsv 2>&1 | tail -2
echo DONE $(date +%T)
