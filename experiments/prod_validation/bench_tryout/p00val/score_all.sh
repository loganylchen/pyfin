#!/usr/bin/env bash
set -uo pipefail
OUT=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val
IMG=quay.io/biocontainers/gffcompare:0.12.10--h9948957_0
# wait for local truth copy to finish
while [ "$(stat -c%s $OUT/truth_full.gtf 2>/dev/null || echo 0)" -lt 1569038390 ]; do sleep 10; done
echo "truth ready $(date)"
for m in cluster_em m2_em cluster; do
  [ -f $OUT/pyfin_$m.gtf ] || continue
  docker run --rm -u $(id -u):$(id -g) -v /SSD:/SSD -w "$OUT" "$IMG" \
    gffcompare -r truth_full.gtf -o gcL_$m pyfin_$m.gtf >/dev/null 2>&1
  eq=$(awk 'NR>1 && $3=="="' $OUT/gcL_$m.pyfin_$m.gtf.tmap 2>/dev/null | wc -l)
  pr=$(grep 'Transcript level' $OUT/gcL_$m.stats 2>/dev/null | grep -oE '[0-9.]+' | tail -1)
  nt=$(grep -c $'\ttranscript\t' $OUT/pyfin_$m.gtf)
  echo "$m: out=$nt = =$eq Pr=$pr"
done
echo "SCORING DONE $(date)"
