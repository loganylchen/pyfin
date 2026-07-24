#!/usr/bin/env bash
set -uo pipefail
OUT=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val
IMG=quay.io/biocontainers/gffcompare:0.12.10--h9948957_0
while squeue -h -j 211268 2>/dev/null | grep -q .; do sleep 30; done
echo "job 211268 finished $(date)"
GTF=$OUT/pyfin_mono.gtf
for i in $(seq 1 30); do [ -f "$GTF" ] && break; sleep 10; done
docker run --rm -u $(id -u):$(id -g) -v /SSD:/SSD -w "$OUT" "$IMG" \
  gffcompare -r truth_full.gtf -o gcL_mono pyfin_mono.gtf >/dev/null 2>&1
eq=$(awk 'NR>1 && $3=="="' $OUT/gcL_mono.pyfin_mono.gtf.tmap 2>/dev/null | wc -l)
pr=$(grep 'Transcript level' $OUT/gcL_mono.stats 2>/dev/null | grep -oE '[0-9.]+' | tail -1)
nt=$(grep -c $'\ttranscript\t' "$GTF")
# single-exon FP count (class not =, single exon) for before/after
echo "MONO-FOLD: out=$nt  ==$eq  Pr=$pr   (baseline sc250: out=9581 ==6323 Pr=66.0)"
echo "SCORING DONE $(date)"
