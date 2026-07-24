#!/usr/bin/env bash
set -uo pipefail
OUT=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val
IMG=quay.io/biocontainers/gffcompare:0.12.10--h9948957_0
# wait for job 211266 to leave the queue
while squeue -h -j 211266 2>/dev/null | grep -q .; do sleep 30; done
echo "job 211266 finished $(date)"
GTF=$OUT/pyfin_ceil_sc250.gtf
for i in $(seq 1 30); do [ -f "$GTF" ] && break; sleep 10; done
docker run --rm -u $(id -u):$(id -g) -v /SSD:/SSD -w "$OUT" "$IMG" \
  gffcompare -r truth_full.gtf -o gcL_ceil_sc250 pyfin_ceil_sc250.gtf >/dev/null 2>&1
eq=$(awk 'NR>1 && $3=="="' $OUT/gcL_ceil_sc250.pyfin_ceil_sc250.gtf.tmap 2>/dev/null | wc -l)
pr=$(grep 'Transcript level' $OUT/gcL_ceil_sc250.stats 2>/dev/null | grep -oE '[0-9.]+' | tail -1)
nt=$(grep -c $'\ttranscript\t' "$GTF")
echo "SC250 CEILING: out=$nt  ==$eq  Pr=$pr   (baseline: out=37078 ==9803 Pr=26.4)"
echo "SCORING DONE $(date)"
