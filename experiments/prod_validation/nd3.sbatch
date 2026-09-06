#!/bin/bash
#SBATCH -J nd3
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 20:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/nd3_%j.out
# The promoted production default (profile real-drna now implies
# assignment_channels=ends + folded_tss_recovery_mode=require), run with NO
# extra flags from frozen snapshot _snap/nd3, on all three grids.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/nd3
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds newdefault $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 10 --profile real-drna \
    --tag _newdefault --out $root/newdefault --tsv $root/newdefault.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
