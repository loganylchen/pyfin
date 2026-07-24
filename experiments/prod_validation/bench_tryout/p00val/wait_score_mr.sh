#!/usr/bin/env bash
set -uo pipefail
OUT=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val
IMG=quay.io/biocontainers/gffcompare:0.12.10--h9948957_0
while squeue -h -j 211275 2>/dev/null | grep -q .; do sleep 30; done
echo "job 211275 finished $(date)"
GTF=$OUT/pyfin_mr.gtf
for i in $(seq 1 30); do [ -f "$GTF" ] && break; sleep 10; done
docker run --rm -u $(id -u):$(id -g) -v /SSD:/SSD -w "$OUT" "$IMG" \
  gffcompare -r truth_full.gtf -o gcL_mr pyfin_mr.gtf >/dev/null 2>&1
eq=$(awk 'NR>1 && $3=="="' $OUT/gcL_mr.pyfin_mr.gtf.tmap 2>/dev/null | wc -l)
pr=$(grep 'Transcript level' $OUT/gcL_mr.stats 2>/dev/null | grep -oE '[0-9.]+' | tail -1)
nt=$(grep -c $'\ttranscript\t' "$GTF")
echo "MONO-RESOLVE: out=$nt  ==$eq  Pr=$pr   (mono-split baseline: out=9146 ==6610 Pr=72.3)"
echo "SCORING DONE $(date)"
