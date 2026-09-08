#!/bin/bash
#SBATCH -J sigaudit
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=180G
#SBATCH -t 24:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/sigaudit_%j.out
# Does raw signal decide junction pairs that the count-based snap cannot?
# Both audit modes are byte-inert: they write junction_election.tsv and
# junction_correction.tsv and pass the original chains through, so the output
# must equal the production default (checked below). full = annotation is the
# truth proxy for "which junction is real"; p00 = the same reads unguided.
cd /SSD/logan/dev/pyfin
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/sigaudit
export PYFIN_PYPATH_PREFIX=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cu122
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
PYFIN_GRID_DATASET=gencode_sub python3 experiments/prod_validation/sirv4_grid_runner.py \
  --dataset gencode_sub --workers 2 --profile real-drna --conds full,p00 \
  --extra="--junction-election-mode audit --junction-correction-mode audit" \
  --tag _gcsigaudit --out $GC/gcsigaudit --tsv $GC/gcsigaudit.tsv 2>&1 | tail -1
echo "=== byte-inertness vs the production default (gcsmc)"
for c in full p00; do
  a=$GC/gcsigaudit/SGNex_H9_directRNA_replicate2_run2/$c/out/assembly.gtf
  b=$GC/gcsmc/SGNex_H9_directRNA_replicate2_run2/$c/out/assembly.gtf
  cmp -s "$a" "$b" && echo "SAME $c" || echo "DIFF $c"
  ls -la $GC/gcsigaudit/SGNex_H9_directRNA_replicate2_run2/$c/out/junction_*.tsv
done
echo DONE $(date +%T)
