#!/bin/bash
#SBATCH -J jc_cand3
#SBATCH -p mne4
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/jc_cand2_%j.out
cd /SSD/logan/dev/pyfin
for CFG in "$JC_CFG"; do set -- $CFG; M=$1; R=$2; K=$3
OUT=experiments/prod_validation/_sirv4grid/jc_additive_m${M}_r${R}_k${K}; mkdir -p $OUT
for S in $(ls /autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2 | grep '^SGNex_H9'); do for C in p00 c_jitter10_10bp c_jitter50_10bp full; do
  [ -s $OUT/$S.$C.summary.json ] && continue
  G=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2/_ref/$S/$C/annotation.gtf; GA=""; [ "$C" != "p00" ] && GA="--gtf $G"
  singularity exec --nv -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif \
    env PYTHONPATH=/SSD/logan/dev/pyfin/experiments/prod_validation/_img/krill016_cu122:/SSD/logan/dev/pyfin \
    python3 experiments/prod_validation/junction_correction_candidate_probe.py --sample $S --cond $C \
      --margin $M --runner-up-margin $R --min-reads $K \
      --bam /autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2/$S/align/$S.sorted.bam --fq /autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2/$S/subset/mapped.fq.gz --blow5 /autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2/$S/subset/mapped.blow5 $GA \
      --genome experiments/prod_validation/_sirv4grid/_ref/SIRV_Set4.fa \
      --truth experiments/prod_validation/_sirv4grid/_ref/SIRV_Set4.gtf \
      --expressed experiments/prod_validation/_sirv4grid/_ref/expressed/$S.gtf \
      --out $OUT > /dev/null 2> $OUT.$S.$C.log || echo "FAIL $CFG $S $C"
  echo "done m$M r$R k$K $S $C $(date +%T)"
done; done; done
