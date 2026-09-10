#!/bin/bash
#SBATCH -J m2norm_heya8
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/m2norm_heya8_%j.out
# Independent-grid replication of the per-candidate-base normalised metric.
# SIRV4 gave +0.08 F1 with a CI straddling zero but a consistent direction
# (60 up / 27 down, every family mean >= 0); this is the replication.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/m2scram
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset heya8 --workers 6 --profile real-drna \
  --extra="--m2-metric sqrt_count_mean_llr" \
  --tag _m2norm --out $R/m2norm --tsv $R/m2norm.tsv 2>&1 | tail -2
echo DONE $(date +%T)
