set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
RATIOS="full p99 p90 p50 p10 p00 c_skip20 c_jitter20_10bp c_spurious20 c_merge20 c_flip20 c_ir20"
end=$(( $(date +%s) + 8*3600 ))
sleep 60  # let the array register
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in $RATIOS; do [ -s "prodfull/$r/pyfin.gtf" ] && n=$((n+1)); done
  q=$(squeue -u "$USER" -h -n pyfin_full 2>/dev/null | grep -c .)
  echo "[$(date +%H:%M:%S)] gtfs=$n/12 pyfin_full_in_queue=$q"
  { [ "$n" -ge 12 ] || [ "$q" -eq 0 ]; } && break
  sleep 120
done
echo "=== scoring $(date) ==="
python3 scorecard_full.py > scorecard_full_result.txt 2>score_err.txt
echo "SCORE_DONE"; cat scorecard_full_result.txt
