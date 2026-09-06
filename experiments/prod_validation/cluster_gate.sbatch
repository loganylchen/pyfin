#!/bin/bash
#SBATCH -J clustgate
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/clustgate_%j.out
# The 5' contradiction now needs a START CLUSTER (>=3 reads within 10 bp), the
# analogue of the polyA-PASS specificity on the 3' side. Two questions:
#   human  -- does TSS agreement recover (the -9 point regression)?
#   spikes -- does the +1.5 F1 the 5' term contributed survive?
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/clust
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
echo "=== human $(date +%T)"
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 3 --profile real-drna \
  --conds full,p50,p00,c_trunc3p50,c_trunc5p50,c_jitter50_10bp,c_merge50 \
  --tag _gcclust --out $GC/gcclust --tsv $GC/gcclust.tsv 2>&1 | tail -1
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --tag _clust --out $root/clust --tsv $root/clust.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
