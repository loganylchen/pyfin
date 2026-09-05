#!/bin/bash
#SBATCH -J final
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/final_%j.out
# Post-review revalidation: the exact tree that will be committed, both grids.
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
echo "=== SEQUIN final $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --tag _final --out $SQ/final --tsv $SQ/final.tsv 2>&1 | tail -2
echo "=== SIRV final $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 \
  --profile real-drna --tag _final --out $S4/final --tsv $S4/final.tsv 2>&1 | tail -2
echo "=== determinism rerun (subset) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 3 \
  --profile real-drna --conds full,c_merge50,c_trunc3p50,p00 \
  --tag _finaldet --out $SQ/final_det --tsv $SQ/final_det.tsv 2>&1 | tail -1
d=0; n=0
for f in $SQ/final_det/*/*/out/assembly.gtf $SQ/final_det/*/*/out/scores.tsv; do
  b=${f/final_det/final}; n=$((n+1)); cmp -s "$f" "$b" || { d=$((d+1)); echo "DIFF $f"; }
done
echo "determinism: $((n-d))/$n identical"
echo DONE $(date +%T)
