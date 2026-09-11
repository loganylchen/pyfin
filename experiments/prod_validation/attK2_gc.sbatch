#!/bin/bash
#SBATCH -J attK2_gc
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/attK2_gc_%j.out
# Operating CURVE for the attestation lever, not a point. Even after K=3
# landed, m2_support_gate still owns 81.8% (3945/4823) of the remaining
# recall loss, and the tier it has not reached (suffix 1-2 reads) is 59%
# truth against an F1 break-even of 21.4%.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/att_land
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 2 --profile real-drna \
  --conds full,p50,p00,c_merge50,c_trunc3p50,c_trunc5p50,c_jitter50_10bp \
  --extra="--m2-support-gate-attest-reads 2" \
  --tag _attK2 --out $GC/attK2 --tsv $GC/attK2.tsv 2>&1 | tail -2
echo DONE $(date +%T)
