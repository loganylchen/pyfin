#!/usr/bin/env bash
set -uo pipefail
OUT=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val
while squeue -h -j 211277 2>/dev/null | grep -q .; do sleep 30; done
echo "job 211277 finished $(date)"
for i in $(seq 1 30); do [ -f "$OUT/pyfin_bid.gtf" ] && break; sleep 10; done
BASE=$OUT/pyfin_ms.gtf; NEW=$OUT/pyfin_bid.gtf
echo "baseline lines: $(wc -l < $BASE)   new lines: $(wc -l < $NEW)"
if diff -q "$BASE" "$NEW" >/dev/null 2>&1; then
  echo "BYTE-IDENTICAL: M0+M1 reproduce the baseline exactly ✓"
else
  echo "DIFFERS — M0/M1 changed output. diff summary:"
  diff "$BASE" "$NEW" | head -20
  echo "... transcript counts: base=$(grep -c $'\ttranscript\t' $BASE) new=$(grep -c $'\ttranscript\t' $NEW)"
fi
echo "DONE $(date)"
