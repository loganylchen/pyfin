#!/bin/bash
#SBATCH -J finalchk
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 8:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/finalchk_%j.out
# Result-neutrality of the sample-evidence rework (Codex rounds 4-11) against
# the nd3 arm that produced the promotion numbers: same profile, same inputs,
# current code. Expect cell-identical scores.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/final
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 --profile real-drna \
  --tag _finalchk --out $SQ/finalchk --tsv $SQ/finalchk.tsv 2>&1 | tail -1
echo DONE $(date +%T)
