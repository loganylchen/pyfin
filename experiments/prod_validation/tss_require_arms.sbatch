#!/bin/bash
#SBATCH -J tssreq
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/tssreq_%j.out
# Spec D5: 5' hazard channel via folded-TSS shadow recovery + positional likelihood.
cd /SSD/logan/dev/pyfin
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
for ds in sequin sirv4; do
  root=$S4; flag=""; [ $ds = sequin ] && { root=$SQ; flag="--dataset sequin"; }
  echo "=== $ds tss require $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --extra="--folded-tss-recovery-mode require" \
    --tag _tssreq --out $root/tssreq --tsv $root/tssreq.tsv 2>&1 | tail -1
  echo "=== $ds tss audit (byte-inertness) $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --conds full,c_merge50,c_trunc5p50,p00 --extra="--folded-tss-recovery-mode audit" \
    --tag _tssaud --out $root/tssaud --tsv $root/tssaud.tsv 2>&1 | tail -1
  d=0; n=0
  for f in $root/tssaud/*/*/out/assembly.gtf $root/tssaud/*/*/out/scores.tsv; do
    b=${f/\/tssaud\//\/final\/}; [ -f "$b" ] || continue; n=$((n+1)); cmp -s "$f" "$b" || { d=$((d+1)); echo "DIFF $f"; }
  done
  echo "$ds audit byte-inert vs final: $((n-d))/$n"
done
echo DONE $(date +%T)
