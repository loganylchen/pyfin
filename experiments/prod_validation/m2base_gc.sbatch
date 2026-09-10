#!/bin/bash
#SBATCH -J m2base_gc
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/m2base_gc_%j.out
# REAL HUMAN (gencode v44, chr1/19/21, H9 rep2 run2) replication of the
# three-arm signal decomposition. This is the grid the spike-ins do not stand
# in for: .ext read-extent siblings are far more common here, so the
# structure-level grouping the control preserves carries much more weight.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/m2scram
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 2 --profile real-drna \
  --conds full,p50,p00,c_merge50,c_trunc3p50,c_trunc5p50,c_jitter50_10bp \
  --extra="" \
  --tag _m2base --out $GC/m2base --tsv $GC/m2base.tsv 2>&1 | tail -2
echo DONE $(date +%T)
