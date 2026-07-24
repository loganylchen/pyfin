#!/bin/bash
set -u
REPO=/SSD/logan/dev/pyfin
B=$REPO/experiments/prod_validation/bench_tryout
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
# wait for job 211278
while squeue -j 211278 -h 2>/dev/null | grep -q .; do sleep 15; done
echo "[$(date)] job 211278 gone"
M2A=$B/p00val/pyfin_m2a.gtf
MS=$B/p00val/pyfin_ms.gtf
if [ ! -s "$M2A" ]; then echo "MISSING $M2A"; exit 2; fi
echo "=== struct_diff m2a vs ms(baseline) ==="
singularity exec -B /SSD "$SIF" /usr/bin/python3.10 "$B/struct_diff.py" "$M2A" "$MS"
echo "EXIT=$?"
