#!/bin/bash
#SBATCH -J aatt_heya8
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=80G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/aatt_heya8_%j.out
# Decision experiment: can the support-gate attestation exemption become a
# DEFAULT? It closes most of the human recall gap as an opt-in; the blocking
# question is what it costs on the spike-ins. One lever, one snapshot (HEAD
# 735cbf8), default arm re-run on the same snapshot.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/att_dec
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset heya8 --workers 6 --profile real-drna \
  --extra="--m2-support-gate-attest-reads 3" \
  --tag _aatt --out $R/aatt --tsv $R/aatt.tsv 2>&1 | tail -2
echo DONE $(date +%T)
