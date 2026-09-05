#!/bin/bash
#SBATCH -J levers
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/levers_%j.out
# Which lever costs recall on the held-out set? Add them one at a time.
cd /SSD/logan/dev/pyfin
export PYFIN_PYPATH_PREFIX=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cpu
SQ=/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid
S4=/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid
run() { # tag out extra extra_gtf dataset profile
  echo "=== $1 $(date +%T)"
  local ds=$5 prof=$6 root=$SQ
  [ "$ds" = sirv4 ] && root=$S4
  local dsflag=""; [ "$ds" = sequin ] && dsflag="--dataset sequin"
  python3 experiments/prod_validation/sirv4_grid_runner.py $dsflag --workers 6 \
    --profile $prof --extra="$3" --extra-with-gtf="$4" \
    --tag _$1 --out $root/$2 --tsv $root/$2.tsv 2>&1 | tail -2
}
P0="--min-gtf-abundance 0 --no-floor-gtf-abundance"
P0G="--require-gtf-positive-mass"
RES="--complete-chain-rescue-reads 2"
ELE="--junction-election-mode require --junction-election-burden defense"
run lv_p0    lv_p0    "$P0"       "$P0G"          sequin real-drna
run lv_p0r   lv_p0r   "$P0 $RES"  "$P0G"          sequin real-drna
run lv_ele   lv_ele   ""          "$ELE"          sequin real-drna
run lv_res   lv_res   "$RES"      ""              sequin real-drna
echo "=== SIRV mirror of the same decomposition $(date +%T)"
run s4_p0    lv_p0    "$P0"       "$P0G"          sirv4 sirv
run s4_p0r   lv_p0r   "$P0 $RES"  "$P0G"          sirv4 sirv
echo DONE $(date +%T)
