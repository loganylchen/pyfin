#!/bin/bash
#SBATCH -J endsarm6
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 16:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/endsarm6_%j.out
# endsarm6: as endsarm4 + only the contradicted end is re-estimated (5-degraded reads keep the annotated TSS), config guard, id disambiguation, confidence pooling.
# contradiction with annotated alternatives) + extent-pair resolution with
# pooled R/ties. Runs from a FROZEN source snapshot (PYFIN_REPO) so edits to
# the working tree cannot leak into a running grid.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/endsarm6
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds assignment-channels ends (snapshot) $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --extra="--assignment-channels ends" --tag _endsarm6 --out $root/endsarm6 --tsv $root/endsarm6.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
