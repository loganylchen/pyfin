set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 4*3600 )); sleep 60
while [ "$(date +%s)" -lt "$end" ]; do
  n=0
  [ -s prodfull_g2/full/pyfin.gtf ] && n=$((n+1))
  [ -s prodfull_g2/c_jitter20_10bp/pyfin.gtf ] && n=$((n+1))
  q=$(squeue -u "$USER" -h -n pyfin_g2 2>/dev/null | grep -c .)
  echo "[$(date +%H:%M:%S)] decisive gtfs=$n/2 in_queue=$q"
  { [ "$n" -ge 2 ] || [ "$q" -eq 0 ]; } && break
  sleep 120
done
echo "=== scoring decisive ablation $(date) ==="
python3 compare_g2.py > compare_g2_decisive.txt 2>compare_g2_derr.txt
echo DECISIVE_DONE
# focused view: full + c_jitter rows + orphan
grep -E "condition|^full |^c_jitter20_10bp |orphan" compare_g2_decisive.txt
