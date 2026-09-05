#!/bin/bash
#SBATCH -J heya8
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/heya8_%j.out
cd /SSD/logan/dev/pyfin
H=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
echo "=== HEYA8 dual, production default $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset heya8 --workers 6 \
  --profile real-drna --tag _h8_final --out $H/final --tsv $H/final.tsv 2>&1 | tail -2
echo "=== expression contract vs competitors $(date +%T)"
python3 experiments/prod_validation/expression_rescore.py --dataset heya8 --workers 8 \
  --pyfin "FINAL=$H/final" 2>&1 | tail -16
echo DONE $(date +%T)
