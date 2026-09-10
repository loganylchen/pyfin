#!/bin/bash
#SBATCH -J m2off_sirv4
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/m2off_sirv4_%j.out
# Falsification arm: is the guided-path M2 signal tie-break (auto -> metric="mean",
# tiebreak margin 1e-09) doing anything at all? One lever: --m2-metric off.
# Everything else is the production real-drna profile.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/m2off
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile real-drna \
  --extra="--m2-metric off" \
  --tag _m2off --out $R/m2off --tsv $R/m2off.tsv 2>&1 | tail -2
echo DONE $(date +%T)
