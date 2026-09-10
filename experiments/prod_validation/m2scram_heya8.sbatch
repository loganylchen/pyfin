#!/bin/bash
#SBATCH -J m2scram_heya8
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/m2scram_heya8_%j.out
# Independent-grid replication of the SIRV4 three-arm signal decomposition.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/m2scram
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset heya8 --workers 6 --profile real-drna \
  --extra="--m2-scramble-seed 20260831" \
  --tag _m2scram --out $R/m2scram --tsv $R/m2scram.tsv 2>&1 | tail -2
echo DONE $(date +%T)
