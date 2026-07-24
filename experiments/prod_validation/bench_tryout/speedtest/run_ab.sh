set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
KGPU=$REPO/experiments/prod_validation/_img/krill_cu122
T=$REPO/experiments/prod_validation/bench_tryout/speedtest
B=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex
NAS=/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2
GENOME=$B/refs/host/GRCh38.primary_assembly.genome.fa
FASTQ=$NAS/fastq/pass.fq.gz; SIGNAL=$NAS/blow5/nanopore.blow5
GTF=$B/results/gencode_full_sweep/_ref/full/annotation.gtf
run () {  # $1=label $2=pythonpath
  local out=$T/$1_work; rm -rf "$out"; mkdir -p "$out"
  local t0=$SECONDS
  singularity exec --nv -B /autofs/mnemosyne3_SSD -B /autofs/NAS25_Shared -B /SSD "$SIF" \
    env PYTHONPATH=$2 /usr/bin/python3.10 -m fin.cli \
    --bam "$T/chr21.bam" --genome "$GENOME" --fastq "$FASTQ" --signal "$SIGNAL" \
    --signal-format slow5 --output-dir "$out" --gpu --quant-mode m2_em --threads 8 \
    --gtf "$GTF" > "$T/$1.log" 2>&1
  echo "$1 wall=$((SECONDS-t0))s tx=$(grep -c $'\ttranscript\t' "$out/assembly.gtf" 2>/dev/null || echo NA)"
}
# GPU-krill run with a background nvidia-smi peak sampler
( for i in $(seq 1 600); do nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits; sleep 0.5; done > "$T/gpu_util.samples" 2>/dev/null ) &
SMPID=$!
echo "[$(date +%T)] === A: GPU krill ==="; run gpu "$KGPU:$REPO"
echo "[$(date +%T)] === B: CPU krill ==="; run cpu "$REPO"
kill $SMPID 2>/dev/null
echo "peak GPU util during A+B window: $(sort -n "$T/gpu_util.samples" 2>/dev/null | tail -1)%"
# confirm which backend each used
echo "-- gpu run krill backend --"; grep -iE "krill|gpu|cpu.*fallback|dtw" "$T/gpu.log" | grep -iE "backend|gpu|cuda|fallback" | head -3
echo "AB_DONE"
