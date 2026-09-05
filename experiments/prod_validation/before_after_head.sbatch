#!/bin/bash
#SBATCH -J beforeaft
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/beforeaft_%j.out
# Unambiguous before/after for the WHOLE uncommitted batch: same runner, same
# flags, only PYFIN_REPO differs (HEAD worktree vs working tree).
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
echo "=== SEQUIN at HEAD (ca82319) $(date +%T)"
PYFIN_REPO=/SSD/logan/dev/pyfin_head python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset sequin --workers 6 --profile real-drna \
  --tag _sq_head --out $SQ/head_real --tsv $SQ/head_real.tsv 2>&1 | tail -2
echo "=== SIRV at HEAD $(date +%T)"
PYFIN_REPO=/SSD/logan/dev/pyfin_head python3 experiments/prod_validation/sirv4_grid_runner.py \
  --workers 6 --profile real-drna \
  --tag _s4_head --out $S4/head_real --tsv $S4/head_real.tsv 2>&1 | tail -2
echo DONE $(date +%T)
