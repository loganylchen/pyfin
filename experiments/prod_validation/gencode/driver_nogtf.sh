#!/usr/bin/env bash
me=$(whoami); MAX=5; DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
N=$(wc -l < $DSROOT/cells_nogtf.txt)
for i in $(seq 0 $((N-1))); do
  s=$(sed -n "$((i+1))p" $DSROOT/cells_nogtf.txt)
  [ -f "$DSROOT/$s/nogtf/pyfin.gtf" ] && continue
  while [ "$(squeue -u $me -h 2>/dev/null|wc -l)" -ge "$MAX" ]; do sleep 30; done
  ( cd $DSROOT/slurm_logs && sbatch --export=ALL,CELL_IDX=$i $DSROOT/run_one_nogtf.sbatch >/dev/null 2>&1 ) && echo "[submit dn] $s"
  sleep 3
done
while [ "$(squeue -u $me -h 2>/dev/null|wc -l)" -gt 0 ]; do sleep 30; done
echo "NOGTF_DONE $(ls $DSROOT/*/nogtf/pyfin.gtf 2>/dev/null|wc -l)/$N"
