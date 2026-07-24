set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 5*3600 )); sleep 90
while [ "$(date +%s)" -lt "$end" ]; do
  [ -s prodfull_unfilt/work/scores.unfiltered.tsv ] && [ -s prodfull_unfilt/work/assembly.gtf ] && break
  squeue -j 211166 -h 2>/dev/null | grep -q . || break
  sleep 120
done
echo "=== unfiltered run done $(date) ==="
ls -la prodfull_unfilt/work/scores.unfiltered.tsv prodfull_unfilt/work/assembly.gtf 2>/dev/null
# sanity: assembly matches prodfull (same config)?
echo "new tx=$(grep -c $'\ttranscript\t' prodfull_unfilt/work/assembly.gtf 2>/dev/null) | prodfull tx=$(grep -c $'\ttranscript\t' prodfull/full/pyfin.gtf 2>/dev/null)"
python3 recall_split.py > recall_split_result.txt 2>recall_split_err.txt
echo RECALL_SPLIT_DONE; cat recall_split_result.txt
