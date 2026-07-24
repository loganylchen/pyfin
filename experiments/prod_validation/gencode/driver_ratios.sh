#!/usr/bin/env bash
me=$(whoami); MAX=5; DSROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
N=$(wc -l < $DSROOT/cells_ratios.txt)
for i in $(seq 0 $((N-1))); do
  c=$(sed -n "$((i+1))p" $DSROOT/cells_ratios.txt); s=$(echo "$c"|awk '{print $1}'); r=$(echo "$c"|awk '{print $2}')
  [ -f "$DSROOT/$s/$r/pyfin.gtf" ] && continue
  while [ "$(squeue -u $me -h 2>/dev/null|wc -l)" -ge "$MAX" ]; do sleep 30; done
  ( cd $DSROOT/slurm_logs && sbatch --export=ALL,CELL_IDX=$i $DSROOT/run_one_ratio.sbatch >/dev/null 2>&1 ) && echo "[sub] $s/$r"
  sleep 3
done
while [ "$(squeue -u $me -h 2>/dev/null|wc -l)" -gt 0 ]; do sleep 30; done
echo "RATIOS_DONE $(ls $DSROOT/*/{p99,p90,p50,p10,c_skip20,c_jitter20_10bp,c_spurious20,c_merge20,c_flip20,c_ir20}/pyfin.gtf 2>/dev/null|wc -l)/$N"
