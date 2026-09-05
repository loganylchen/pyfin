#!/bin/bash
#SBATCH -J xds2
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/xds2_%j.out
cd /SSD/logan/dev/pyfin
export PYFIN_PYPATH_PREFIX=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cpu
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
BASE="--min-gtf-abundance 0 --no-floor-gtf-abundance --complete-chain-rescue-reads 2"
BASEG="--require-gtf-positive-mass --junction-election-mode require --junction-election-burden defense"
echo "=== SIRV off $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile sirv \
  --tag _xd_s4_off --out $S4/xd_off --tsv $S4/xd_off.tsv 2>&1 | tail -2
echo "=== SIRV best $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile sirv \
  --extra "$BASE" --extra-with-gtf "$BASEG" \
  --tag _xd_s4_best --out $S4/xd_best --tsv $S4/xd_best.tsv 2>&1 | tail -2
echo "=== SIRV expression rescore $(date +%T)"
python3 experiments/prod_validation/sirv4_expression_rescore.py --workers 8 2>&1 | tail -3
echo DONE $(date +%T)
