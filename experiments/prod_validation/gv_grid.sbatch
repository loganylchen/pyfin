#!/bin/bash
#SBATCH -J gv_grid
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gv_grid_%j.out
cd /SSD/logan/dev/pyfin
export PYFIN_PYPATH_PREFIX=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cpu
D=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
BASE="--min-gtf-abundance 0 --no-floor-gtf-abundance --complete-chain-rescue-reads 2"
BASEG="--require-gtf-positive-mass --junction-election-mode require --junction-election-burden defense"
run() { # tag outdir extra extra_gtf [samples] [conds]
  local tag=$1 out=$2 extra=$3 extrag=$4 samples=${5:-} conds=${6:-}
  local args=""; [ -n "$samples" ] && args="$args --samples $samples"; [ -n "$conds" ] && args="$args --conds $conds"
  echo "=== $tag $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $args --extra "$extra" --extra-with-gtf "$extrag" \
    --tag _$tag --out $D/$out --tsv $D/$out.tsv 2>&1 | tail -3
}
run gv_off gv_off "" ""
run gv_best gv_best "$BASE" "$BASEG"
run gv_tss_audit gv_tss_audit "$BASE --folded-tss-recovery-mode audit" "$BASEG"
run gv_tss_req gv_tss_req "$BASE --folded-tss-recovery-mode require" "$BASEG"
run gv_tss_req_det gv_tss_req_det "$BASE --folded-tss-recovery-mode require" "$BASEG" SGNex_H9_directRNA_replicate2_run1 full,p00,c_trunc5p10
echo "=== audit inertness + determinism cmp"
for s in $D/gv_tss_audit/SGNex_*; do sn=$(basename $s); for c in $s/*/; do cn=$(basename $c); for f in assembly.gtf scores.tsv abundance_refit.json; do
  cmp -s $D/gv_best/$sn/$cn/out/$f $c/out/$f || echo "AUDIT_DIFF $sn $cn $f"; done; done; done | sort | uniq -c | head
for c in full p00 c_trunc5p10; do for f in assembly.gtf scores.tsv abundance_refit.json tss_shadow_recovery.tsv; do
  a=$D/gv_tss_req/SGNex_H9_directRNA_replicate2_run1/$c/out/$f; b=$D/gv_tss_req_det/SGNex_H9_directRNA_replicate2_run1/$c/out/$f
  if cmp -s "$a" "$b"; then echo "DET_IDENTICAL $c/$f"; else echo "DET_DIFF $c/$f"; fi; done; done
echo "=== rescore $(date +%T)"
python3 experiments/prod_validation/sirv4_expression_rescore.py --workers 8 2>&1 | tail -3
echo "DONE $(date +%T)"
