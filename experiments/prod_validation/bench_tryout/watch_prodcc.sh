set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  squeue -u "$USER" -h -n pyfin_prodcc 2>/dev/null | grep -q . || break
  sleep 300
done
python3 eval_honest.py prodfull_prodcc gc_prodcc > eval_honest_PRODCC.txt 2>/dev/null
echo "=== DEFAULT-ON containment on PURE PRODUCTION defaults vs OLD production (honF1) ==="
printf "%-16s | %-18s | %-18s | %s\n" condition "new honPr/rec/F1" "old honPr/rec/F1" dF1
for c in full p10 p00 c_skip20 c_jitter20_10bp c_merge20; do
  n=$(grep "^$c " eval_honest_PRODCC.txt | awk '{printf "%s/%s/%s",$5,$6,$7}')
  o=$(grep "^$c " eval_honest_OFF.txt | awk '{printf "%s/%s/%s",$5,$6,$7}')
  nf=$(grep "^$c " eval_honest_PRODCC.txt | awk '{print $7}')
  of=$(grep "^$c " eval_honest_OFF.txt | awk '{print $7}')
  d=$(python3 -c "print(f'{($nf)-($of):+.1f}')" 2>/dev/null || echo "?")
  printf "%-16s | %-18s | %-18s | %s\n" "$c" "${n:-NA}" "${o:-NA}" "$d"
done
echo "(want: dF1 >= 0 everywhere; this is the isolated default-on effect on real production defaults)"
echo PRODCC_DONE
