#!/bin/bash
#SBATCH -J attK1_sirv4
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/attK1_sirv4_%j.out
# Spike-in cost of the same operating point.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/att_land
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile real-drna \
  --extra="--m2-support-gate-attest-reads 1" \
  --tag _attK1 --out $R/attK1 --tsv $R/attK1.tsv 2>&1 | tail -2
echo DONE $(date +%T)
