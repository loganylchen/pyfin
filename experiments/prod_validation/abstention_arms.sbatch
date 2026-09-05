#!/bin/bash
#SBATCH -J abstain
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/abstain_%j.out
# Does "abstention never convicts" remove the depth dependence?
# prosecution burden = drop only when ballots actually convict.
cd /SSD/logan/dev/pyfin
export PYFIN_PYPATH_PREFIX=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cpu
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
P0="--min-gtf-abundance 0 --no-floor-gtf-abundance"
RES="--complete-chain-rescue-reads 2"
PROS="--require-gtf-positive-mass --junction-election-mode require --junction-election-burden prosecution"
NOGTFMASS="--junction-election-mode require --junction-election-burden prosecution"
echo "=== SEQUIN P0+rescue+election(prosecution) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --extra "$P0 $RES" --extra-with-gtf "$PROS" \
  --tag _sq_pros --out $SQ/pros --tsv $SQ/pros.tsv 2>&1 | tail -2
echo "=== SEQUIN same but WITHOUT require-gtf-positive-mass $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --extra "$P0 $RES" --extra-with-gtf "$NOGTFMASS" \
  --tag _sq_pros_nm --out $SQ/pros_nm --tsv $SQ/pros_nm.tsv 2>&1 | tail -2
echo "=== SIRV P0+rescue+election(prosecution) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 \
  --profile real-drna --extra "$P0 $RES" --extra-with-gtf "$PROS" \
  --tag _s4_pros --out $S4/pros --tsv $S4/pros.tsv 2>&1 | tail -2
echo "=== SIRV same but WITHOUT require-gtf-positive-mass $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 \
  --profile real-drna --extra "$P0 $RES" --extra-with-gtf "$NOGTFMASS" \
  --tag _s4_pros_nm --out $S4/pros_nm --tsv $S4/pros_nm.tsv 2>&1 | tail -2
echo DONE $(date +%T)
