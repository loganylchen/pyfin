#!/bin/bash
#SBATCH -J annoraw
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=120G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/annoraw_%j.out
# The annotation-scale contrast WITHOUT cProfile. The profiled run inflated
# Python-heavy loops (parse looked like 137 s/call; unprofiled it is ~28 s), so
# the attribution has to be redone on real wall time.
set -u
cd /SSD/logan/dev/pyfin
P=/SSD/logan/dev/pyfin/experiments/prod_validation/_annoprobe
S=SGNex_H9_directRNA_replicate2_run2
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
G=$B/refs/host/GRCh38.primary_assembly.genome.fa
FULL=$B/results/gencode_full_sweep_ratio_grid_v2/_ref/$S/full/annotation.gtf
SUB=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/_sub/_ref/$S/full/annotation.gtf
FQ=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/reads.fq.gz
SIG=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/signal.blow5
run () {
  local tag=$1 anno=$2
  echo "=== $tag START $(date +%s) $(date +%T) ==="
  singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
    env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 -m fin.cli \
    --bam $P/small.bam --genome $G --fastq $FQ --signal $SIG \
    --output-dir $P/raw_$tag --profile real-drna --no-gpu --threads 2 \
    --gtf $anno --m2-support-gate-attest-reads 1 > $P/raw_$tag.log 2>&1
  echo "=== $tag DONE $(date +%s) $(date +%T) intervals=$(grep -oE 'Generated [0-9]+ intervals' $P/raw_$tag.log | tail -1) ==="
}
run subanno  $SUB
run fullanno $FULL
echo ALLDONE $(date +%T)
