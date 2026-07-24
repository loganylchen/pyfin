set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for p in cc_a cc_b; do [ -s "prodfull_$p/p00/pyfin.gtf" ] && n=$((n+1)); done
  squeue -u "$USER" -h -n pyfin_cc2 2>/dev/null | grep -q . || break
  [ "$n" -ge 2 ] && break
  sleep 240
done
echo "=== CONTAINMENT-CLUSTER test (p00) [tx F1std orphan honPr corrRec honF1] ==="
for p in cc_a cc_b; do
  python3 eval_honest.py prodfull_$p gc_$p 2>/dev/null | grep "^p00" | sed "s/^p00/$p /"
done
echo "smart5 ref : 7268  3.5  36.9  39.7  34.3  36.8"
echo "OFF    ref : 5102  2.5  42.7  37.6  22.8  28.4"
echo "competitors: isoquant hPr54.4 rec34.4 F1 42.1 | bambu 38.9/40.0/39.4 | stringtie3 33.8/45.3/38.7"
echo CC2_DONE
