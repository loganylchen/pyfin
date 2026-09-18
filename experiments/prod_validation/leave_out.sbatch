#!/bin/bash
#SBATCH -J leaveout
#SBATCH -p mne4
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH -t 2:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/leaveout_%j.out
set -uo pipefail
cd /SSD/logan/dev/pyfin
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
S=SGNex_H9_directRNA_replicate2_run2
singularity exec --pwd /tmp -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
  env PYTHONPATH=/SSD/logan/dev/pyfin PYFIN_SAMPLE_EVIDENCE_CACHE=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/_sample_evidence_cache \
  /usr/bin/python3.10 /SSD/logan/dev/pyfin/experiments/prod_validation/measure_leave_out_effect.py \
  --stage /SSD/logan/dev/pyfin/experiments/prod_validation/_gwstage/$S \
  --gtf $B/results/gencode_full_sweep_ratio_grid_v2/_ref/$S/full/annotation.gtf \
  --genome $B/refs/host/GRCh38.primary_assembly.genome.fa \
  --evidence /SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign/gw/SGNex_H9_directRNA_replicate2_run2/c_spurious50/out/_sample_evidence_bdffe5a7960eee80.json --n 60 --out /SSD/logan/dev/pyfin/experiments/prod_validation/_leaveout 2>&1 | { grep -v "INFO" || true; }
RC=${PIPESTATUS[0]}
# A stale leave_out_effect.json from an earlier run would otherwise sit under a
# job SLURM reports as COMPLETED. The measurement's own status decides.
if [ "$RC" -ne 0 ]; then echo "FAILED rc=$RC $(date +%T)"; exit 1; fi
echo DONE $(date +%T)
