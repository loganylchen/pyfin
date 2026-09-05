#!/bin/bash
#SBATCH -J floors
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/floors_%j.out
# The fixed read-count floor is the suspected source of the depth dependence.
# Replace it with "must have evidence (mass > 0)" and measure on BOTH grids.
cd /SSD/logan/dev/pyfin
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
EVID="--min-abundance 0 --min-gtf-abundance 0 --no-floor-gtf-abundance"
for ds in sequin sirv4; do
  root=$SQ; flag="--dataset sequin"
  [ $ds = sirv4 ] && { root=$S4; flag=""; }
  echo "=== $ds: evidence-only floors (mass>0, no read-count floor) $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 \
    --profile real-drna --extra="$EVID" --extra-with-gtf="--require-gtf-positive-mass" \
    --tag _${ds}_evid --out $root/evid --tsv $root/evid.tsv 2>&1 | tail -2
  echo "=== $ds: novel floor 0 only $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 \
    --profile real-drna --extra="--min-abundance 0" \
    --tag _${ds}_nf0 --out $root/nf0 --tsv $root/nf0.tsv 2>&1 | tail -2
  echo "=== $ds: no relative gates (isoform-fraction 0) $(date +%T)"
  python3 experiments/prod_validation/sirv4_grid_runner.py $flag --workers 6 \
    --profile real-drna --extra="--min-isoform-fraction 0" \
    --tag _${ds}_nif --out $root/nif --tsv $root/nif.tsv 2>&1 | tail -2
done
echo DONE $(date +%T)
