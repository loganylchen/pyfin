#!/bin/bash
#SBATCH -J attK1_heya8
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/attK1_heya8_%j.out
# The frozen promotion rule needs a positive (or inert) result on EVERY grid
# separately. Human K=1 beats the landed K=3 by +0.97 F1 (5 up / 1 down) and
# SIRV4 is byte-identical across K; HEYA8 is the last grid unmeasured.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/att_land
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset heya8 --workers 6 --profile real-drna \
  --extra="--m2-support-gate-attest-reads 1" \
  --tag _attK1 --out $R/attK1 --tsv $R/attK1.tsv 2>&1 | tail -2
echo DONE $(date +%T)
