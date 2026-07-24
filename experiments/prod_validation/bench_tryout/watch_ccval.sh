set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  squeue -u "$USER" -h -n pyfin_ccval 2>/dev/null | grep -q . || break
  sleep 300
done
echo "=== cc_b config MULTI-CONDITION validation vs PRODUCTION (honF1) ==="
echo "condition           cc_b(honPr/corrRec/honF1)      OFF(honPr/corrRec/honF1)"
python3 eval_honest.py prodfull_ccval gc_ccval > eval_honest_CCVAL.txt 2>/dev/null
for c in full p10 c_skip20 c_jitter20_10bp c_merge20; do
  n=$(grep "^$c " eval_honest_CCVAL.txt | awk '{printf "%s/%s/%s",$4,$5,$6}')
  o=$(grep "^$c " eval_honest_OFF.txt | awk '{printf "%s/%s/%s",$4,$5,$6}')
  printf "%-18s  %-28s  %s\n" "$c" "${n:-MISSING}" "${o:-?}"
done
echo "p00 (from cc_b): 42.7/33.8/37.8   vs OFF 37.6/22.8/28.4"
echo "NICE-TO-HOLD: c_jitter is pyfin's niche (must stay strong vs isoquant ~13)."
echo CCVAL_DONE
