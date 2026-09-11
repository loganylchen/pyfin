#!/bin/bash
#SBATCH -J flank_art
#SBATCH -p mne4
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH -t 6:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/flank_art_%j.out
cd /SSD/logan/dev/pyfin
SIF="singularity exec -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif env PYTHONPATH=/SSD/logan/dev/pyfin /usr/bin/python3.10"
G=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2
H=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/heya8_dual_full_sweep_ratio_grid_v2
S=SGNex_H9_directRNA_replicate4_run2
$SIF experiments/prod_validation/splice_flank_artifact.py \
  --gtf $G/expression_contract/_truth/SGNex_H9_directRNA_replicate2_run2.gtf \
  --genome /autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/GRCh38.primary_assembly.genome.fa \
  --bam $G/_sub/input.bam --chroms chr1,chr19,chr21 --label human
echo
HS=SGNex_HEYA8_directRNA_replicate1_run1
$SIF experiments/prod_validation/splice_flank_artifact.py \
  --gtf $H/_ref/$HS/full/annotation.gtf \
  --genome /SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid/_ref/spikein_sequin_sirv1.fa \
  --bam $H/$HS/align/$HS.sorted.bam --label HEYA8
echo DONE $(date +%T)
