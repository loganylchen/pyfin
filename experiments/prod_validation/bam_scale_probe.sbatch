#!/bin/bash
#SBATCH -J bamscale
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=150G
#SBATCH -t 2:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/bamscale_%j.out
# Same annotation (the full 253k one), same threads, ONLY the BAM differs.
# Tests whether per-interval cost scales with the library's TOTAL read count
# rather than with the reads in the interval -- the flat 4.0 s median held
# across a 50x range of per-interval read counts, so the driver is global.
# Each arm is killed after a fixed wall time; only the run.log timings matter.
set -u
cd /SSD/logan/dev/pyfin
P=/SSD/logan/dev/pyfin/experiments/prod_validation/_annoprobe
S=SGNex_H9_directRNA_replicate2_run2
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
G=$B/refs/host/GRCh38.primary_assembly.genome.fa
ANNO=$B/results/gencode_full_sweep_ratio_grid_v2/_ref/$S/full/annotation.gtf
BIG=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/input.bam
FQ=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/reads.fq.gz
SIG=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/signal.blow5
arm () {
  local tag=$1 bam=$2 secs=$3
  echo "=== $tag bam_reads=$(samtools idxstats $bam | awk '{s+=$3}END{print s}') START $(date +%T) ==="
  timeout $secs singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
    env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10 -m fin.cli \
    --bam $bam --genome $G --fastq $FQ --signal $SIG \
    --output-dir $P/bs_$tag.out --profile real-drna --no-gpu --threads 2 \
    --gtf $ANNO --m2-support-gate-attest-reads 1 > $P/bs_$tag.log 2>&1
  echo "=== $tag STOPPED $(date +%T) intervals=$(grep -c 'Completed interval' $P/bs_$tag.log) ==="
}
arm smallbam $P/small.bam 1500
arm fullbam  $BIG        1500
echo ALLDONE $(date +%T)
