#!/bin/bash
#SBATCH -J tesarms
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/tesarms_%j.out
# The intermediate that decides the production default: 3'/polyA channel only.
# Human GENCODE says the 5' start term costs 9 points of TSS agreement while the
# 3' repair is what fixes truncated ends. If `tes` keeps the spike-in F1 gain,
# it is the better default than `ends`.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/tesarm
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds tes-only $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --extra="--assignment-channels tes" --tag _tesonly --out $root/tesonly --tsv $root/tesonly.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
