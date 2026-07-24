set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
mkdir -p gc_g2
end=$(( $(date +%s) + 8*3600 )); sleep 60
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in full p99 p90 p50 p10 p00 c_skip20 c_jitter20_10bp c_spurious20 c_merge20 c_flip20 c_ir20; do [ -s "prodfull_g2/$r/pyfin.gtf" ] && n=$((n+1)); done
  q=$(squeue -u "$USER" -h -n pyfin_g2 2>/dev/null | grep -c .)
  echo "[$(date +%H:%M:%S)] g2 gtfs=$n/12 in_queue=$q"
  { [ "$n" -ge 12 ] || [ "$q" -eq 0 ]; } && break
  sleep 120
done
echo "=== scoring ablation $(date) ==="
python3 compare_g2.py > compare_g2_result.txt 2>compare_g2_err.txt
echo COMPARE_DONE; cat compare_g2_result.txt
