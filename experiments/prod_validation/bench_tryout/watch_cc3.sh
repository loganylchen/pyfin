set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  squeue -u "$USER" -h -n pyfin_cc3 2>/dev/null | grep -q . || break
  sleep 300
done
python3 eval_honest.py prodfull_cc3 gc_cc3 > eval_honest_CC3.txt 2>/dev/null
echo "=== NEW chain-cluster discovery (full pipeline) vs OLD production (honF1) ==="
printf "%-16s | %-18s | %-18s | %s\n" condition "NEW honPr/rec/F1" "OLD honPr/rec/F1" dF1
for c in full p00 c_jitter20_10bp; do
  n=$(grep "^$c " eval_honest_CC3.txt | awk '{printf "%s/%s/%s",$5,$6,$7}')
  o=$(grep "^$c " eval_honest_OFF.txt | awk '{printf "%s/%s/%s",$5,$6,$7}')
  nf=$(grep "^$c " eval_honest_CC3.txt | awk '{print $7}'); of=$(grep "^$c " eval_honest_OFF.txt | awk '{print $7}')
  d=$(python3 -c "print(f'{($nf)-($of):+.1f}')" 2>/dev/null || echo "?")
  printf "%-16s | %-18s | %-18s | %s\n" "$c" "${n:-NA}" "${o:-NA}" "$d"
done
echo "(NEW = chain-cluster discovery default-on + all existing gates; OLD = production)"
echo CC3_DONE
