#!/usr/bin/env bash
# Submit sirv4 cells to SLURM keeping <=MAX jobs in the queue at any time.
set -uo pipefail
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4
MAX=5
N=$(wc -l < "$ROOT/cells.txt")
me=$(whoami)
for i in $(seq 0 $((N-1))); do
  cell=$(sed -n "$((i+1))p" "$ROOT/cells.txt")
  s=$(echo "$cell"|awk '{print $1}'); r=$(echo "$cell"|awk '{print $2}')
  [ -f "$ROOT/$s/$r/pyfin.gtf" ] && { echo "[skip] $s/$r"; continue; }
  # wait for a free slot
  while [ "$(squeue -u "$me" -h 2>/dev/null | wc -l)" -ge "$MAX" ]; do sleep 20; done
  sbatch --export=ALL,CELL_IDX=$i "$ROOT/run_one.sbatch" >/dev/null 2>&1 \
    && echo "[submit] idx=$i $s/$r (inq=$(squeue -u "$me" -h 2>/dev/null|wc -l))" \
    || echo "[FAIL submit] idx=$i $s/$r"
  sleep 3
done
# wait for the last ones to drain
while [ "$(squeue -u "$me" -h 2>/dev/null | wc -l)" -gt 0 ]; do sleep 20; done
echo "DRIVER_DONE done=$(ls $ROOT/*/*/pyfin.gtf 2>/dev/null|wc -l)/$N"
