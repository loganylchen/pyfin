#!/bin/bash
#SBATCH -J aatt_gc
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/aatt_gc_%j.out
# The recall lever, on the grid that defines the target: +2.1 pp Sn at constant
# Pr closes every cell where stringtie3 is ahead. Default control is m2base
# (same behaviour -- verified byte-identical assembly at the default).
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/att_dec
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 2 --profile real-drna \
  --conds full,p50,p00,c_merge50,c_trunc3p50,c_trunc5p50,c_jitter50_10bp \
  --extra="--m2-support-gate-attest-reads 3" \
  --tag _aatt --out $GC/aatt --tsv $GC/aatt.tsv 2>&1 | tail -2
echo DONE $(date +%T)
