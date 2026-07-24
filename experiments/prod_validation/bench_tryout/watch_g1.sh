set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
RATIOS="full p99 p90 p50 p10 p00 c_skip20 c_jitter20_10bp c_spurious20 c_merge20 c_flip20 c_ir20"
end=$(( $(date +%s) + 20*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in $RATIOS; do [ -s "prodfull_g1/$r/pyfin.gtf" ] && n=$((n+1)); done
  q=$(squeue -u "$USER" -h -n pyfin_g1 2>/dev/null | grep -c .)
  echo "[$(date +%H:%M)] g1 gtfs=$n/12 in_queue=$q"
  { [ "$n" -ge 12 ] || [ "$q" -eq 0 ]; } && break
  sleep 300
done
echo "=== eval min_reads=1 $(date) ==="
python3 eval_honest.py prodfull_g1 gc_g1 > eval_honest_G1.txt 2>>g1_err.txt
python3 campaign_rank.py prodfull_g1 gc_g1 > campaign_rank_G1.txt 2>>g1_err.txt
python3 campaign_conditions.py prodfull_g1 gc_g1 > campaign_conditions_G1.txt 2>>g1_err.txt
echo G1_EVAL_DONE
echo "===== OFF vs min_reads=1 : honest F1 per condition ====="
paste <(sed 1d eval_honest_OFF.txt) <(sed 1d eval_honest_G1.txt) 2>/dev/null | awk '{printf "%-16s OFF honF1=%-6s  ON honF1=%-6s\n",$1,$7,$14}'
echo "===== composite rank (min_reads=1) ====="; cat campaign_rank_G1.txt
