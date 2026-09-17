#!/bin/bash
#SBATCH -J epgate
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 8:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/epgate_%j.out
# Codex's layer-3 gate for the reader-reuse patch: byte-identity of the three
# published artifacts, at --threads 1 AND --threads 2 (the spawned-worker path),
# on the LARGE annotation. Baseline is the pre-patch tree via PYFIN_REPO.
set -u
cd /SSD/logan/dev/pyfin
P=/SSD/logan/dev/pyfin/experiments/prod_validation/_epgate
mkdir -p $P
S=SGNex_H9_directRNA_replicate2_run2
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
G=$B/refs/host/GRCh38.primary_assembly.genome.fa
ANNO=$B/results/gencode_full_sweep_ratio_grid_v2/_ref/$S/full/annotation.gtf
BAM=/SSD/logan/dev/pyfin/experiments/prod_validation/_annoprobe/small.bam
FQ=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/reads.fq.gz
SIG=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/signal.blow5
run () {
  local tag=$1 repo=$2 threads=$3
  echo "=== $tag repo=$repo threads=$threads START $(date +%T) ==="
  singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
    env PYTHONPATH=$repo /usr/bin/python3.10 -m fin.cli \
    --bam $BAM --genome $G --fastq $FQ --signal $SIG \
    --output-dir $P/$tag --profile real-drna --no-gpu --threads $threads \
    --gtf $ANNO --m2-support-gate-attest-reads 1 > $P/$tag.log 2>&1
  echo "=== $tag DONE rc=$? $(date +%T) parses=$(grep -c 'Parsed .* genes and' $P/$tag.log) ==="
}
BASE=$P/_baseline_src
rm -rf $BASE; mkdir -p $BASE
git archive HEAD | tar -x -C $BASE
run base_t1 $BASE 1
run new_t1  /SSD/logan/dev/pyfin 1
run base_t2 $BASE 2
run new_t2  /SSD/logan/dev/pyfin 2
echo "=== BYTE IDENTITY ==="
for pair in "base_t1 new_t1" "base_t2 new_t2"; do
  set -- $pair
  for f in assembly.gtf scores.tsv abundance_refit.json; do
    if cmp -s $P/$1/$f $P/$2/$f; then echo "IDENTICAL $1 vs $2 : $f"; else echo "DIFFERS   $1 vs $2 : $f"; fi
  done
done
echo ALLDONE $(date +%T)
