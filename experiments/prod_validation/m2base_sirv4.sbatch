#!/bin/bash
#SBATCH -J m2base_sirv4
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/m2base_sirv4_%j.out
# DEFAULT control on the SAME source snapshot, so the comparison is one lever.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/m2off
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile real-drna \
  --tag _m2base --out $R/m2base --tsv $R/m2base.tsv 2>&1 | tail -2
echo DONE $(date +%T)
