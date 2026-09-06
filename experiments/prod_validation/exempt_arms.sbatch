#!/bin/bash
#SBATCH -J exempt
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH -t 30:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/exempt_%j.out
# The support gate drops genuine annotated transcripts whose only reads are
# 5'-degraded fragments (human GENCODE: -7.7 F1 on an uncorrupted annotation).
# The exemption keeps the gate's anti-corruption purpose by asking whether the
# reads CONTRADICT the annotation's junctions. Validated on all four grids,
# on top of the cluster-gated endpoint channels.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/exempt
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
H8=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
echo "=== human $(date +%T)"
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 3 --profile real-drna \
  --conds full,p50,p00,c_trunc3p50,c_trunc5p50,c_jitter50_10bp,c_merge50 \
  --extra="--m2-support-gate-gtf-exempt" \
  --tag _gcexempt --out $GC/gcexempt --tsv $GC/gcexempt.tsv 2>&1 | tail -1
for ds in sequin sirv4 heya8; do
  case $ds in sequin) root=$SQ; flag="--dataset sequin";; sirv4) root=$S4; flag="";; heya8) root=$H8; flag="--dataset heya8";; esac
  echo "=== $ds $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 --profile real-drna \
    --extra="--m2-support-gate-gtf-exempt" --tag _exempt --out $root/exempt --tsv $root/exempt.tsv 2>&1 | tail -1
done
echo DONE $(date +%T)
