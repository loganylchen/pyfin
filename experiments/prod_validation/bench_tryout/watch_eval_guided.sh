set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 5*3600 )); sleep 90
while [ "$(date +%s)" -lt "$end" ]; do
  [ -s prodfull_g2/full/pyfin.gtf ] && [ -s prodfull_g2/c_jitter20_10bp/pyfin.gtf ] && break
  squeue -u "$USER" -h -n pyfin_g2 2>/dev/null | grep -q . || break
  sleep 120
done
echo "=== eval_guided $(date) ==="
python3 eval_guided.py > eval_guided_result.txt 2>eval_guided_err.txt
echo EVAL_GUIDED_DONE; cat eval_guided_result.txt
