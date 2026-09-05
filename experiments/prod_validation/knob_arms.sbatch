#!/bin/bash
#SBATCH -J knobs
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 16:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/knobs_%j.out
# One-at-a-time: which knob that differs between the two profiles is good on
# BOTH grids? Anything good on only one is a dataset fit, not an algorithm.
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
arm() { # name flags
  local name=$1 flags=$2
  for ds in sequin sirv4; do
    local root=$SQ flag="--dataset sequin"
    [ $ds = sirv4 ] && { root=$S4; flag=""; }
    [ -s $root/k_$name.tsv ] && continue
    echo "=== $ds $name $(date +%T)"
    python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 \
      --profile real-drna --extra="$flags" \
      --tag _k_$name --out $root/k_$name --tsv $root/k_$name.tsv 2>&1 | tail -2
  done
}
arm softmass  "--max-soft-mass-ratio 2.0"
arm fulllen   "--min-fulllen-fraction 0.1"
arm metric    "--m2-metric mean"
arm nostrict  "--no-strict-novel-abundance-floor"
arm ab3       "--min-abundance 3"
echo DONE $(date +%T)
