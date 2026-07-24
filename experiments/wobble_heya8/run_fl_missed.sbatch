#!/bin/bash
#SBATCH --job-name=fl_missed
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
W=$REPO/experiments/wobble_heya8/cf_measure
TBAM=$W/transcriptome.bam
GSRC=$REPO/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/input.bam
GBAM=$W/genome.bam
TRUTH=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/gencode.v44.primary_assembly.annotation.gtf
echo "[$(date)] copy genome BAM local (resolve NFS symlink)"
cp -L "$GSRC" "$GBAM"; cp -L "$GSRC.bai" "$GBAM.bai"
echo "[$(date)] diagnosing FL>=3 missed"
singularity exec -B /SSD -B /autofs "$SIF" env PYTHONPATH=$REPO /usr/bin/python3.10 \
  $REPO/experiments/wobble_heya8/fl_missed_diagnose.py "$TBAM" "$GBAM" "$TRUTH" "$W/gc_cf.cf.gtf.refmap"
echo "[$(date)] DONE"
