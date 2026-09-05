#!/bin/bash
#SBATCH -J ce_sirv
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/ce_sirv_%j.out
cd /SSD/logan/dev/pyfin
D=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
echo "=== SIRV4 candidate-existence AUDIT (3 conds x 6 samples) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile sirv \
  --conds full,c_jitter10_10bp,c_merge10 --extra "--candidate-existence-mode audit" \
  --tag _ce_audit --out $D/ce_audit --tsv $D/ce_audit.tsv 2>&1 | tail -2
echo "=== reference without audit $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile sirv \
  --conds full,c_jitter10_10bp,c_merge10 \
  --tag _ce_ref --out $D/ce_ref --tsv $D/ce_ref.tsv 2>&1 | tail -2
for s in $D/ce_audit/SGNex_*; do sn=$(basename $s); for c in $s/*/; do cn=$(basename $c)
  for f in assembly.gtf scores.tsv abundance_refit.json; do
    a=$D/ce_ref/$sn/$cn/out/$f; b=$c/out/$f
    [ -f "$a" ] && [ -f "$b" ] && { cmp -s "$a" "$b" && echo IDENTICAL || echo "DIFF $sn/$cn/$f"; }
  done; done; done | sort | uniq -c | sed 's/^ *//'
wc -l $D/ce_audit/SGNex_*/*/out/candidate_existence.tsv 2>/dev/null | tail -3
echo DONE $(date +%T)
