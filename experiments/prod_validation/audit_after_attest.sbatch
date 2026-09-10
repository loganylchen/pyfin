#!/bin/bash
#SBATCH -J audit_att
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/audit_att_%j.out
# With attestation now on by default, WHO owns the remaining recall loss?
# Same selection audit that identified m2_support_gate as the previous owner.
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep_ratio_grid_v2
GC=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
S=SGNex_H9_directRNA_replicate2_run2
O=$GC/_audit_after_attest
rm -rf "$O"; mkdir -p "$O"
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/att_land \
  /usr/bin/python3.10 -m fin.cli \
  --bam $GC/_sub/input.bam --genome $GC/_ref/GRCh38.primary_assembly.genome.fa \
  --fastq /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/reads.fq.gz \
  --signal /SSD/logan/dev/pyfin/experiments/prod_validation/gencode/_p00val/stage/signal.blow5 \
  --gtf $GC/_sub/_ref/$S/full/annotation.gtf \
  --output-dir "$O" --profile real-drna --no-gpu --threads 6 --selection-audit \
  > "$O/run.log" 2>&1
echo "rc=$? $(date +%T)"
