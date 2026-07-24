#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
SAMPLES=(SGNex_HEYA8_directRNA_replicate1_run1 SGNex_HEYA8_directRNA_replicate1_run2 SGNex_HEYA8_directRNA_replicate2_run1 SGNex_HEYA8_directRNA_replicate2_run2 SGNex_HEYA8_directRNA_replicate3_run1)
mkdir -p run_logs
worker() {
  local RATIOS=(full p00 c_jitter10bp c_skip10 c_spurious5 c_merge5 c_flip5 c_ir10); local s=$1
  for r in "${RATIOS[@]}"; do
    local OUT="matrix/$s/$r/out_v3"
    [[ -f "$OUT/assembly.gtf" ]] && continue
    mkdir -p "$OUT"; echo "  ▶ [$(date +%H:%M:%S)] $s/$r/v3"
    docker run --rm --gpus all -u "$(id -u):$(id -g)" -v "$ROOT":/data -w /data pyfin-gpu:dev \
      fin --bam "matrix/$s/stage/input.bam" --gtf "matrix/_ref/$r/annotation.gtf" \
          --genome "matrix/$s/stage/genome.fa" --fastq "matrix/$s/stage/reads.fq.gz" \
          --signal "matrix/$s/stage/signal.blow5" --output-dir "$OUT" --gpu --quant-mode m2_em \
      >"run_logs/${s}_${r}_v3.log" 2>&1 && echo "    ✓ $s/$r/v3" || echo "  ❌ $s/$r/v3"
  done
  echo "  [$(date +%H:%M:%S)] $s DONE"
}
export -f worker; export ROOT
printf '%s\n' "${SAMPLES[@]}" | parallel -j5 --line-buffer worker
echo "=== V3 DONE: $(ls -1 matrix/SGNex_HEYA8_directRNA_replicate*/*/out_v3/assembly.gtf 2>/dev/null|wc -l)/40 ==="
