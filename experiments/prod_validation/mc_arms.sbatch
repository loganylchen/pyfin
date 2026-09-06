#!/bin/bash
#SBATCH -J mcarms
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/mcarms_%j.out
# The DRS-specific step: score the NESTED junction contrast from raw signal,
# charging a candidate for the events it cannot explain. Human: the signal's
# decision rate doubles (2.6% -> 5.5%) and F1 +1.25 (full) / +0.72 (p50).
# Spike-in validation on all three grids.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/mc
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --extra="--m2-containment-contrast" --tag _mc --out $root/mc --tsv $root/mc.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
