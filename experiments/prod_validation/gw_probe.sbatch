#!/bin/bash
#SBATCH -J gw_probe
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=120G
#SBATCH -t 8:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gw_probe_%j.out
# COST PROBE for the genome-wide campaign: one cell, current production
# default, to replace a 6x extrapolation with a measurement.
cd /SSD/logan/dev/pyfin
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
echo "START $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode --workers 4 \
  --profile real-drna --extra="--m2-support-gate-attest-reads 1" \
  --samples SGNex_H9_directRNA_replicate2_run2 --conds full \
  --tag _gwK1 --out $R/gwK1 --tsv $R/gwK1.tsv 2>&1 | tail -4
echo "DONE $(date +%T)"
