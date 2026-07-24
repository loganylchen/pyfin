#!/usr/bin/env bash
# Two configs per cell:
#   out_bp20cas : bp20 + cassette, displace OFF  (finish to 40, explicit --no-displace)
#   out_v2      : bp20 + cassette + displace ON  (new full default; 40 cells)
# 5-way parallel over samples, serial inside. Idempotent (skip existing assembly.gtf).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
SAMPLES=(
  SGNex_HEYA8_directRNA_replicate1_run1 SGNex_HEYA8_directRNA_replicate1_run2
  SGNex_HEYA8_directRNA_replicate2_run1 SGNex_HEYA8_directRNA_replicate2_run2
  SGNex_HEYA8_directRNA_replicate3_run1
)
mkdir -p run_logs
worker() {
  local RATIOS=(full p00 c_jitter10bp c_skip10 c_spurious5 c_merge5 c_flip5 c_ir10)
  local s=$1
  for r in "${RATIOS[@]}"; do
    # config bp20cas (displace OFF)
    local OUT="matrix/$s/$r/out_bp20cas"
    if [[ ! -f "$OUT/assembly.gtf" ]]; then
      mkdir -p "$OUT"
      echo "  ▶ [$(date +%H:%M:%S)] $s/$r/bp20cas"
      docker run --rm --gpus all -u "$(id -u):$(id -g)" -v "$ROOT":/data -w /data pyfin-gpu:dev \
        fin --bam "matrix/$s/stage/input.bam" --gtf "matrix/_ref/$r/annotation.gtf" \
            --genome "matrix/$s/stage/genome.fa" --fastq "matrix/$s/stage/reads.fq.gz" \
            --signal "matrix/$s/stage/signal.blow5" --output-dir "$OUT" --gpu --quant-mode m2_em \
            --no-m2-cluster-recheck-novel-displaces-gtf \
        >"run_logs/${s}_${r}_bp20cas.log" 2>&1 && echo "    ✓ $s/$r/bp20cas" || echo "  ❌ $s/$r/bp20cas"
    fi
    # config v2 (displace ON = all current defaults)
    OUT="matrix/$s/$r/out_v2"
    if [[ ! -f "$OUT/assembly.gtf" ]]; then
      mkdir -p "$OUT"
      echo "  ▶ [$(date +%H:%M:%S)] $s/$r/v2"
      docker run --rm --gpus all -u "$(id -u):$(id -g)" -v "$ROOT":/data -w /data pyfin-gpu:dev \
        fin --bam "matrix/$s/stage/input.bam" --gtf "matrix/_ref/$r/annotation.gtf" \
            --genome "matrix/$s/stage/genome.fa" --fastq "matrix/$s/stage/reads.fq.gz" \
            --signal "matrix/$s/stage/signal.blow5" --output-dir "$OUT" --gpu --quant-mode m2_em \
        >"run_logs/${s}_${r}_v2.log" 2>&1 && echo "    ✓ $s/$r/v2" || echo "  ❌ $s/$r/v2"
    fi
  done
  echo "  [$(date +%H:%M:%S)] sample $s DONE"
}
export -f worker; export ROOT
printf '%s\n' "${SAMPLES[@]}" | parallel -j5 --line-buffer worker
bp20cas=$(ls -1 matrix/SGNex_HEYA8_directRNA_replicate*/*/out_bp20cas/assembly.gtf 2>/dev/null | wc -l)
v2=$(ls -1 matrix/SGNex_HEYA8_directRNA_replicate*/*/out_v2/assembly.gtf 2>/dev/null | wc -l)
echo "=== FINAL DONE: bp20cas=$bp20cas/40  v2=$v2/40 ==="
