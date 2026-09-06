#!/bin/bash
#SBATCH -J gcprobe
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gcprobe_%j.out
# Timing/feasibility probe: one human GENCODE cell with the promoted default.
cd /SSD/logan/dev/pyfin
export PYFIN_GRID_DATASET=gencode
export PYFIN_REPO=/SSD/logan/dev/pyfin
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
date +%T
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode --workers 1 \
  --profile real-drna --conds full --tag _probe --out $GC/probe --tsv $GC/probe.tsv 2>&1 | tail -3
date +%T
