set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  squeue -u "$USER" -h -n pyfin_pp,pyfin_wc 2>/dev/null | grep -q . || break
  sleep 300
done
echo "=== FINAL p00 frontier [tx | honPr | corrRec | honF1] ==="
printf "%-12s %8s %7s %8s %7s\n" config tx honPr corrRec honF1
row() {
  python3 eval_honest.py "$1" "gc_$2" 2>/dev/null | awk -v n="$2" '/^p00/{printf "%-12s %8s %7s %8s %7s\n",n,$2,$5,$6,$7}'
}
row prodfull OFF
row prodfull_smart5 smart5
row prodfull_pp1 pp1
row prodfull_pp2 pp2
row prodfull_pp3 pp3
row prodfull_wc_r10 wc_r10
row prodfull_wc_r05 wc_r05
echo "--- competitors ---"
printf "%-12s %8s %7s %8s %7s\n" isoquant 5309 54.4 34.4 42.1
printf "%-12s %8s %7s %8s %7s\n" bambu 8661 38.9 40.0 39.4
printf "%-12s %8s %7s %8s %7s\n" stringtie3 11275 33.8 45.3 38.7
echo FINAL_DONE
