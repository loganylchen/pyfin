#!/bin/bash
#SBATCH -J bridge
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/bridge_%j.out
# Evidence-gated annotation bridging: does it capture the recall of unconditional
# bridging WITHOUT the merge/jitter precision damage? Both grids, default flags.
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
echo "=== SEQUIN evidence-gated bridging $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
  --profile real-drna --tag _sq_gb --out $SQ/gated_bridge --tsv $SQ/gated_bridge.tsv 2>&1 | tail -2
echo "=== SIRV evidence-gated bridging $(date +%T)"
python3 experiments/prod_validation/sirv4_grid_runner.py --workers 6 \
  --profile real-drna --tag _s4_gb --out $S4/gated_bridge --tsv $S4/gated_bridge.tsv 2>&1 | tail -2
echo DONE $(date +%T)
