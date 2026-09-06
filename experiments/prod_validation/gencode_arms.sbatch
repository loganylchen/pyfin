#!/bin/bash
#SBATCH -J gcarms
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH -t 48:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gcarms_%j.out
# Does the promoted default generalise from spike-ins to a REAL human
# transcriptome? Same two arms as the spike-in grids, on GENCODE v44 / GRCh38:
#   gcold  = the pre-promotion default (endpoint channels off)
#   gcnew  = the promoted default (no extra flags)
# Conditions: the two endpoint corruptions the channels target, a junction
# corruption as a control, the uncorrupted annotation, a half annotation, and
# de novo.
cd /SSD/logan/dev/pyfin
export PYFIN_GRID_DATASET=gencode_sub
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/gcarms
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
CONDS=full,p50,p00,c_trunc3p50,c_trunc5p50,c_jitter50_10bp,c_merge50
echo "=== gcnew (promoted default) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode_sub --workers 3 \
  --profile real-drna --conds $CONDS \
  --tag _gcnew --out $GC/gcnew --tsv $GC/gcnew.tsv 2>&1 | tail -2
echo "=== gcold (pre-promotion default) $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode_sub --workers 3 \
  --profile real-drna --conds $CONDS \
  --extra="--assignment-channels none --folded-tss-recovery-mode off" \
  --tag _gcold --out $GC/gcold --tsv $GC/gcold.tsv 2>&1 | tail -2
echo DONE $(date +%T)
