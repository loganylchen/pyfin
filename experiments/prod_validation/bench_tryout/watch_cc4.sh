set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  squeue -u "$USER" -h -n pyfin_cc4 2>/dev/null | grep -q . || break
  sleep 300
done
python3 eval_honest.py prodfull_cc4 gc_cc4 > eval_honest_CC4.txt 2>/dev/null
echo "=== NEW+read-guarded collapse (gates as-is) vs OLD production (honF1) ==="
printf "%-16s | %-10s | %-10s | %-10s\n" condition "NEW+guard" "OLD" "dF1(vs OLD)"
for c in full p00 c_jitter20_10bp; do
  nf=$(grep "^$c " eval_honest_CC4.txt | awk '{print $7}')
  of=$(grep "^$c " eval_honest_OFF.txt | awk '{print $7}')
  cf=$(grep "^$c " eval_honest_CC3.txt 2>/dev/null | awk '{print $7}')
  d=$(python3 -c "print(f'{($nf)-($of):+.1f}')" 2>/dev/null || echo "?")
  printf "%-16s | %-10s | %-10s | %-10s  (unguarded cc3 F1=%s)\n" "$c" "${nf:-NA}" "${of:-NA}" "$d" "${cf:-NA}"
done
echo CC4_DONE
