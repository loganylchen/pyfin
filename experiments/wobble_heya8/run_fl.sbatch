#!/bin/bash
#SBATCH --job-name=fl_rescore
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
W=$REPO/experiments/wobble_heya8/cf_measure
SRC=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep/SGNex_H9_directRNA_replicate2_run2/full/align/SGNex_H9_directRNA_replicate2_run2.transcriptome.bam
TRUTH=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/gencode.v44.primary_assembly.annotation.gtf
echo "[$(date)] copy transcriptome BAM local"
cp "$SRC" "$W/transcriptome.bam"
[ -f "$SRC.bai" ] && cp "$SRC.bai" "$W/transcriptome.bam.bai" || \
  singularity exec -B /SSD "$SIF" samtools index "$W/transcriptome.bam"
CFNOUT=$(grep -c $'\ttranscript\t' $W/cf.gtf)
PRODNOUT=$(grep -c $'\ttranscript\t' $W/prod.gtf)
echo "[$(date)] compute FL support + rescore (cf_nout=$CFNOUT prod_nout=$PRODNOUT)"
singularity exec -B /SSD -B /autofs "$SIF" env PYTHONPATH=$REPO /usr/bin/python3.10 \
  $REPO/experiments/wobble_heya8/compute_fl.py \
  "$W/transcriptome.bam" "$TRUTH" "$W/gc_cf.cf.gtf.refmap" "$W/gc_prod.prod.gtf.refmap" "$CFNOUT" "$PRODNOUT"
echo "[$(date)] DONE"
