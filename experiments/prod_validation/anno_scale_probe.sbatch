#!/bin/bash
#SBATCH -J annoprobe
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 6:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/annoprobe_%j.out
# Is the 5.4x per-interval slowdown paid for by ANNOTATION SIZE? One small BAM,
# identical in both arms; the only variable is the annotation the cell is given.
# Profiled, so the answer comes with attribution rather than a yes/no.
set -u
cd /SSD/logan/dev/pyfin
P=/SSD/logan/dev/pyfin/experiments/prod_validation/_annoprobe
S=SGNex_H9_directRNA_replicate2_run2
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
G=$B/refs/host/GRCh38.primary_assembly.genome.fa
SMALL=$B/results/gencode_full_sweep_ratio_grid_v2/_ref/$S/full/annotation.gtf
SUB=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/_sub/_ref/$S/full/annotation.gtf
FQ=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/reads.fq.gz
SIG=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/signal.blow5
run () {
  local tag=$1 anno=$2
  echo "=== $tag START $(date +%T) anno=$(stat -c%s $anno) ==="
  singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
    env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 -m cProfile -o $P/$tag.prof \
    -m fin.cli --bam $P/small.bam --genome $G --fastq $FQ --signal $SIG \
    --output-dir $P/$tag.out --profile real-drna --no-gpu --threads 2 \
    --gtf $anno --m2-support-gate-attest-reads 1 > $P/$tag.log 2>&1
  echo "=== $tag DONE rc=$? $(date +%T) ==="
}
run subanno  $SUB
run fullanno $SMALL
echo ALLDONE $(date +%T)
