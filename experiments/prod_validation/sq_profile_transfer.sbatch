#!/bin/bash
#SBATCH -J sq_prof
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/sq_prof_%j.out
cd /SSD/logan/dev/pyfin
D=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
for P in sirv real-drna; do
  echo "=== profile=$P $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py --dataset sequin --workers 6 \
    --profile $P --tag _sq_${P//-/} --out $D/off_$P --tsv $D/off_$P.tsv 2>&1 | tail -2
done
echo "=== quant axis"
for P in sirv real-drna; do
  python3 experiments/prod_validation/sequin_quant_score.py --runs $D/off_$P \
    --tag pyfin_off_$P --out $D/quant_off_$P.tsv
done
echo DONE $(date +%T)
