#!/bin/bash
#SBATCH -J ce_disc
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/ce_disc_%j.out
cd /SSD/logan/dev/pyfin
D=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
echo "=== SIRV4 discriminative-power audit: 50% corruption conditions $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 --profile sirv \
  --conds c_jitter50_10bp,c_flip50,c_skip50,c_merge50,c_trunc5p50,c_trunc3p50 \
  --extra "--candidate-existence-mode audit" \
  --tag _ce_disc --out $D/ce_disc --tsv $D/ce_disc.tsv 2>&1 | tail -2
echo "=== Sequin same conditions $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna \
  --conds c_jitter50_10bp,c_flip50,c_skip50,c_merge50,c_trunc5p50,c_trunc3p50 \
  --extra "--candidate-existence-mode audit" \
  --tag _ce_disc --out $SQ/ce_disc --tsv $SQ/ce_disc.tsv 2>&1 | tail -2
echo DONE $(date +%T)
