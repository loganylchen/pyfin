#!/usr/bin/env bash
set -uo pipefail
DSROOT=$1; PREREQ=${2:-}
me=$(whoami); MAX=5
mkdir -p "$DSROOT/slurm_logs"
[ -n "$PREREQ" ] && { echo "[driver] waiting for prereq: $PREREQ"; while [ ! -f "$PREREQ" ]; do sleep 30; done; }
N=$(wc -l < "$DSROOT/cells.txt")
for i in $(seq 0 $((N-1))); do
  cell=$(sed -n "$((i+1))p" "$DSROOT/cells.txt"); s=$(echo "$cell"|awk '{print $1}'); r=$(echo "$cell"|awk '{print $2}')
  [ -f "$DSROOT/$s/$r/pyfin.gtf" ] && continue
  while [ "$(squeue -u "$me" -h 2>/dev/null | wc -l)" -ge "$MAX" ]; do sleep 20; done
  ( cd "$DSROOT/slurm_logs" && sbatch --export=ALL,DSROOT="$DSROOT",CELL_IDX=$i /SSD/logan/dev/pyfin/experiments/prod_validation/run_one_generic.sbatch >/dev/null 2>&1 ) \
    && echo "[submit] $s/$r" || echo "[FAIL] $s/$r"
  sleep 3
done
while [ "$(squeue -u "$me" -h 2>/dev/null | wc -l)" -gt 0 ]; do sleep 20; done
echo "DRIVER_DONE $(ls $DSROOT/*/*/pyfin.gtf 2>/dev/null|wc -l)/$N"
