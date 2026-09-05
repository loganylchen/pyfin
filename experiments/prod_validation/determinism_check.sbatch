#!/bin/bash
#SBATCH -J det
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 3:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/det_%j.out
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
echo "=== determinism: two independent reruns, different worker counts $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --conds full,c_merge50,c_trunc3p50,p00 \
  --tag _det1 --out $SQ/det1 --tsv $SQ/det1.tsv 2>&1 | tail -1
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 3 \
  --profile real-drna --conds full,c_merge50,c_trunc3p50,p00 \
  --tag _det2 --out $SQ/det2 --tsv $SQ/det2.tsv 2>&1 | tail -1
d=0; n=0
for f in $SQ/det1/*/*/out/assembly.gtf; do
  b=${f/det1/det2}; n=$((n+1)); cmp -s "$f" "$b" || { d=$((d+1)); echo "DIFF $f"; }
done
for f in $SQ/det1/*/*/out/scores.tsv $SQ/det1/*/*/out/abundance_refit.json; do
  b=${f/det1/det2}; n=$((n+1)); cmp -s "$f" "$b" || { d=$((d+1)); echo "DIFF $f"; }
done
echo "determinism: $((n-d))/$n identical"
echo "=== mass conservation"
python3 - <<'PY'
import json, glob
bad=0; n=0
for f in glob.glob("/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid/det1/*/*/out/abundance_refit.json"):
    d=json.load(open(f)); n+=1
    e=abs(float(d.get("mass_balance_error", 0.0)))
    if e > 1e-9: bad+=1; print("mass error", e, f)
print(f"mass balance: {n-bad}/{n} within 1e-9")
PY
echo DONE $(date +%T)
