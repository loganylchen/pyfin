#!/bin/bash
#SBATCH -J arm8
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 20:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/arm8_%j.out
# arm8: frozen snapshot _snap/arm8 = ec3d7f9 + structure-level M2 gate (group
# share conserved, fusion identity). Two arms per grid: default8 (channels none;
# measures the default-path change of the gate) and endsarm8 (channels ends).
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/arm8
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds default8 $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 5 --profile real-drna \
    --tag _default8 --out $root/default8 --tsv $root/default8.tsv 2>&1 | tail -1 &
  echo "=== $ds endsarm8 $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 5 --profile real-drna \
    --extra="--assignment-channels ends" --tag _endsarm8 --out $root/endsarm8 --tsv $root/endsarm8.tsv 2>&1 | tail -1 &
  wait
done
echo DONE $(date +%T)
