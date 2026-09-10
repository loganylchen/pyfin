#!/bin/bash
#SBATCH -J m2norm_sirv4
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/m2norm_sirv4_%j.out
# The third of the NanoSplicer prescription, never benchmarked on the
# production path: per-candidate-base NORMALISED score over the same tight
# differing-junction window, instead of production's undivided sum (the form
# with the known longer-template bias). One lever vs m2base.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/m2scram
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile real-drna \
  --extra="--m2-metric sqrt_count_mean_llr" \
  --tag _m2norm --out $R/m2norm --tsv $R/m2norm.tsv 2>&1 | tail -2
echo DONE $(date +%T)
