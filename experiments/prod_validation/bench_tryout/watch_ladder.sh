set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  squeue -u "$USER" -h -n pyfin_ladder 2>/dev/null | grep -q . || break
  sleep 300
done
echo "=== M1/M2/M3 DISCRIMINATION LADDER — per-mode precision/recall frontier ==="
for M in argmax m1_em m2_em; do
  # build gffcompare tracking for the broad output first
  python3 eval_honest.py prodfull_ladder_$M gc_ladder_$M >/dev/null 2>&1 || true
  python3 frontier_of.py prodfull_ladder_$M gc_ladder_$M "$M" 2>/dev/null || echo "$M: score failed"
  echo
done
echo "Compare corrRec at structPr~=68 across modes: higher = better read-mass concentration."
echo "ref: isoquant structPr91/corrRec34.4 | smart5(m2_em,ab5) 66/34.3 honF1 36.8"
echo LADDER_DONE
