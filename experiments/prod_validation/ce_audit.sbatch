#!/bin/bash
#SBATCH -J ce_audit
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 3:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/ce_audit_%j.out
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
echo "=== Sequin, candidate-existence AUDIT, 3 conditions x 4 samples $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --conds full,c_jitter10_10bp,c_merge10 \
  --extra "--candidate-existence-mode audit" \
  --tag _ce_audit --out $SQ/ce_audit --tsv $SQ/ce_audit.tsv 2>&1 | tail -2
echo "=== same cells WITHOUT the audit (byte-inertness reference) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --conds full,c_jitter10_10bp,c_merge10 \
  --tag _ce_ref --out $SQ/ce_ref --tsv $SQ/ce_ref.tsv 2>&1 | tail -2
echo "=== byte-inertness"
for s in $SQ/ce_audit/SGNex_*; do sn=$(basename $s); for c in $s/*/; do cn=$(basename $c)
  for f in assembly.gtf scores.tsv abundance_refit.json; do
    a=$SQ/ce_ref/$sn/$cn/out/$f; b=$c/out/$f
    if [ -f "$a" ] && [ -f "$b" ]; then cmp -s "$a" "$b" && echo "IDENTICAL $sn/$cn/$f" || echo "DIFF $sn/$cn/$f"; fi
  done; done; done | sort | uniq -c | sed 's/^ *//'
echo "=== determinism rerun of one cell"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 2 \
  --profile real-drna --samples SGNex_K562_directRNA_replicate4_run1 --conds full \
  --extra "--candidate-existence-mode audit" \
  --tag _ce_det --out $SQ/ce_det --tsv $SQ/ce_det.tsv 2>&1 | tail -1
A=$SQ/ce_audit/SGNex_K562_directRNA_replicate4_run1/full/out/candidate_existence.tsv
B=$SQ/ce_det/SGNex_K562_directRNA_replicate4_run1/full/out/candidate_existence.tsv
cmp -s $A $B && echo "DET_IDENTICAL candidate_existence.tsv" || echo "DET_DIFF candidate_existence.tsv"
echo "=== audit content"
head -1 $A | tr '\t' '\n' | head -3
wc -l $SQ/ce_audit/SGNex_*/*/out/candidate_existence.tsv 2>/dev/null | tail -5
echo DONE $(date +%T)
