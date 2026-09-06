#!/bin/bash
#SBATCH -J gctes
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gctes_%j.out
# The intermediate that was never benchmarked head-to-head: 3'/polyA channel
# only (no 5' start term). The human endpoint-accuracy measurement says the 3'
# repair works and the 5' term is what costs TSS agreement.
cd /SSD/logan/dev/pyfin
export PYFIN_GRID_DATASET=gencode_sub
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/gcarms
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
CONDS=full,p50,p00,c_trunc3p50,c_trunc5p50,c_jitter50_10bp,c_merge50
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode_sub --workers 3 \
  --profile real-drna --conds $CONDS --extra="--assignment-channels tes" \
  --tag _gctes --out $GC/gctes --tsv $GC/gctes.tsv 2>&1 | tail -2
echo DONE $(date +%T)
