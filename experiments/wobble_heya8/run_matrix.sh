#!/usr/bin/env bash
# Orchestrate pyfin matrix: 5 samples × 8 ratios × 2 configs (bp10, bp20).
# Parallel over samples (one sample = one worker chain, 16 cells serial inside).
# Each cell idempotent: skip if assembly.gtf exists.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

SAMPLES=(
  SGNex_HEYA8_directRNA_replicate1_run1
  SGNex_HEYA8_directRNA_replicate1_run2
  SGNex_HEYA8_directRNA_replicate2_run1
  SGNex_HEYA8_directRNA_replicate2_run2
  SGNex_HEYA8_directRNA_replicate3_run1
)

mkdir -p run_logs

worker() {
  # Arrays defined inside (bash export -f doesn't carry arrays).
  local RATIOS=(full p00 c_jitter10bp c_skip10 c_spurious5 c_merge5 c_flip5 c_ir10)
  local CONFIGS=(bp10 bp20)
  local s=$1
  for r in "${RATIOS[@]}"; do
    for c in "${CONFIGS[@]}"; do
      local OUT="matrix/$s/$r/out_$c"
      local LOG="run_logs/${s}_${r}_${c}.log"
      local BP
      case "$c" in bp10) BP=10 ;; bp20) BP=20 ;; esac
      if [[ -f "$OUT/assembly.gtf" ]]; then continue; fi
      mkdir -p "$OUT"
      echo "  ▶ [$(date +%H:%M:%S)] $s / $r / $c"
      if docker run --rm --gpus all -u "$(id -u):$(id -g)" \
          -v "$ROOT":/data -w /data pyfin-gpu:dev \
          fin --bam "matrix/$s/stage/input.bam" \
              --gtf "matrix/_ref/$r/annotation.gtf" \
              --genome "matrix/$s/stage/genome.fa" \
              --fastq "matrix/$s/stage/reads.fq.gz" \
              --signal "matrix/$s/stage/signal.blow5" \
              --output-dir "$OUT" --gpu --quant-mode m2_em \
              --m2-cluster-recheck-bp "$BP" \
          >"$LOG" 2>&1; then
        echo "    ✓ [$(date +%H:%M:%S)] $s / $r / $c"
      else
        echo "  ❌ FAIL: $s/$r/$c (see $LOG)" | tee -a run_logs/_failures.log
      fi
    done
  done
  echo "  [$(date +%H:%M:%S)] sample $s DONE"
}
export -f worker
export ROOT

started=$(date +%s)
printf '%s\n' "${SAMPLES[@]}" | parallel -j5 --line-buffer worker
elapsed=$(( $(date +%s) - started ))
total=$(ls -1 matrix/SGNex_HEYA8_directRNA_replicate*/*/out_bp*/assembly.gtf 2>/dev/null | wc -l)
echo "=== MATRIX DONE: cells_with_gtf=$total/80  elapsed=${elapsed}s ==="
