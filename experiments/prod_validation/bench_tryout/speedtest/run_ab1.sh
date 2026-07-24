set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
KGPU=$REPO/experiments/prod_validation/_img/krill_cu122
T=$REPO/experiments/prod_validation/bench_tryout/speedtest
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
NAS=/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2
run () { local out=$T/$1_1t; rm -rf "$out"; mkdir -p "$out"; local t0=$SECONDS
  singularity exec --nv -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
    env PYTHONPATH=$2 /usr/bin/python3.10 -m fin.cli \
    --bam "$T/chr21sub.bam" --genome "$B/refs/host/GRCh38.primary_assembly.genome.fa" \
    --fastq "$NAS/fastq/pass.fq.gz" --signal "$NAS/blow5/nanopore.blow5" --signal-format slow5 \
    --output-dir "$out" --gpu --quant-mode m2_em --threads 1 \
    --gtf "$B/results/gencode_full_sweep/_ref/full/annotation.gtf" > "$T/$1_1t.log" 2>&1
  echo "$1_1t wall=$((SECONDS-t0))s"; }
run gpu "$KGPU:$REPO"; run cpu "$REPO"
echo "AB1_DONE"
