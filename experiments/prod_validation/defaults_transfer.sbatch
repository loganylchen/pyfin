#!/bin/bash
#SBATCH -J defaults
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/defaults_%j.out
# The two arms nobody has ever measured:
#  (a) the ACTUAL production default (--profile real-drna, all levers off) on SIRV4
#  (b) the current best stack on SEQUIN, the held-out set
cd /SSD/logan/dev/pyfin
export PYFIN_PYPATH_PREFIX=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cpu
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
BASE="--min-gtf-abundance 0 --no-floor-gtf-abundance --complete-chain-rescue-reads 2"
BASEG="--require-gtf-positive-mass --junction-election-mode require --junction-election-burden defense"
echo "=== (b) SEQUIN best stack, real-drna profile $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --extra "$BASE" --extra-with-gtf "$BASEG" \
  --tag _sq_best --out $SQ/best_real --tsv $SQ/best_real.tsv 2>&1 | tail -2
echo "=== (b2) SEQUIN best stack, sirv profile (isolates profile vs levers) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile sirv --extra "$BASE" --extra-with-gtf "$BASEG" \
  --tag _sq_best_sirv --out $SQ/best_sirv --tsv $SQ/best_sirv.tsv 2>&1 | tail -2
echo "=== (a) SIRV4 with the ACTUAL production default profile $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 \
  --profile real-drna --tag _s4_realdrna --out $S4/off_real --tsv $S4/off_real.tsv 2>&1 | tail -2
echo "=== (a2) SIRV4 best stack under real-drna profile $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 \
  --profile real-drna --extra "$BASE" --extra-with-gtf "$BASEG" \
  --tag _s4_best_real --out $S4/best_real --tsv $S4/best_real.tsv 2>&1 | tail -2
echo "=== rescore both grids $(date +%T)"
python3 experiments/prod_validation/expression_rescore.py --dataset sequin --workers 8 \
  --pyfin "pyfin_sq_off_real=$SQ/off_real-drna" --pyfin "pyfin_sq_off_sirv=$SQ/off_sirv" \
  --pyfin "pyfin_sq_best_real=$SQ/best_real" --pyfin "pyfin_sq_best_sirv=$SQ/best_sirv" 2>&1 | tail -16
python3 experiments/prod_validation/expression_rescore.py --dataset sirv4 --workers 8 \
  --pyfin "pyfin_s4_off_sirv=$S4/xd_off" --pyfin "pyfin_s4_best_sirv=$S4/xd_best" \
  --pyfin "pyfin_s4_off_real=$S4/off_real" --pyfin "pyfin_s4_best_real=$S4/best_real" 2>&1 | tail -16
echo DONE $(date +%T)
