#!/bin/bash
#SBATCH -J mc2
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH -t 30:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/mc2_%j.out
# The nested junction contrast, now adjudicated only between ANNOTATED chains
# (two transcripts), never between a read-derived chain and its container (one
# transcript's 5' degradation ladder, an extent question). Human vs SIRV4 were
# the two poles: +1.25 and -0.37.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/mc2
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
echo "=== human $(date +%T)"
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 3 --profile real-drna --conds full,p50,c_trunc5p50 \
  --extra="--m2-containment-contrast" --tag _mc2 --out $GC/mc2 --tsv $GC/mc2.tsv 2>&1 | tail -1
for ds in sirv4 sequin heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --extra="--m2-containment-contrast" --tag _mc2 --out $root/mc2 --tsv $root/mc2.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
