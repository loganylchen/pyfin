#!/bin/bash
#SBATCH -J ce3
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/ce3_%j.out
# Phase 2 audit (spec D2): 3'/polyA channel. Must be byte-inert; must stratify
# c_merge / c_trunc3p corrupted annotation where the structural channel could not.
cd /SSD/logan/dev/pyfin
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
CONDS=c_merge50,c_trunc3p50,c_trunc5p50,c_jitter50_10bp,c_skip50,c_flip50,full,p00
for ds in sequin sirv4; do
  root=$S4; flag=""; [ $ds = sequin ] && { root=$SQ; flag="--dataset sequin"; }
  echo "=== $ds audit $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --conds $CONDS --extra="--candidate-existence-mode audit --candidate-existence-bootstrap 1023" \
    --tag _ce3 --out $root/ce3 --tsv $root/ce3.tsv 2>&1 | tail -1
  echo "=== $ds reference (no audit) $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --conds $CONDS --tag _ce3ref --out $root/ce3ref --tsv $root/ce3ref.tsv 2>&1 | tail -1
  d=0; n=0
  for f in $root/ce3/*/*/out/assembly.gtf $root/ce3/*/*/out/scores.tsv $root/ce3/*/*/out/abundance_refit.json; do
    b=${f/\/ce3\//\/ce3ref\/}; n=$((n+1)); cmp -s "$f" "$b" || { d=$((d+1)); echo "DIFF $f"; }
  done
  echo "$ds byte-inert: $((n-d))/$n identical"
done
echo DONE $(date +%T)
