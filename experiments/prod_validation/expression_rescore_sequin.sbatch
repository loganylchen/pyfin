#!/bin/bash
#SBATCH -J expr_sq
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/expr_sq_%j.out
cd /SSD/logan/dev/pyfin
python3 experiments/prod_validation/expression_rescore.py --dataset sequin --workers 8 \
  --pyfin "pyfin_sirvprof=experiments/prod_validation/_sequingrid/off_sirv" \
  --pyfin "pyfin_realdrna=experiments/prod_validation/_sequingrid/off_real-drna"
