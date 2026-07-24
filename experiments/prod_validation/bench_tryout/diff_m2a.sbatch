#!/bin/bash
#SBATCH --job-name=m2a_diff
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
REPO=/SSD/logan/dev/pyfin
B=$REPO/experiments/prod_validation/bench_tryout
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
echo "=== struct_diff pyfin_m2a vs pyfin_ms(baseline) ==="
singularity exec -B /SSD "$SIF" /usr/bin/python3.10 "$B/struct_diff.py" \
  "$B/p00val/pyfin_m2a.gtf" "$B/p00val/pyfin_ms.gtf"
echo "DIFF_EXIT=$?"
