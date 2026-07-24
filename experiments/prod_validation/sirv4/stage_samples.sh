#!/usr/bin/env bash
set -euo pipefail
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/sirv4
SW=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep
stage_one() {
  local s=$1; local D=$ROOT/$s/stage; mkdir -p "$D"
  [[ -f $D/input.bam ]]     || cp "$SW/$s/align/${s}.sorted.bam" "$D/input.bam"
  [[ -f $D/input.bam.bai ]] || { cp "$SW/$s/align/${s}.sorted.bam.bai" "$D/input.bam.bai" 2>/dev/null || samtools index "$D/input.bam"; }
  [[ -f $D/signal.blow5 ]]     || cp "$SW/$s/subset/mapped.blow5" "$D/signal.blow5"
  [[ -f $D/signal.blow5.idx ]] || cp "$SW/$s/subset/mapped.blow5.idx" "$D/signal.blow5.idx" 2>/dev/null || true
  local nc=$(find "$SW/$s" -name nanocount.tsv 2>/dev/null | head -1)
  [[ -f $D/nanocount.tsv ]] || cp "$nc" "$D/nanocount.tsv"
  ln -sf ../../_ref/genome.fa "$D/genome.fa"; ln -sf ../../_ref/genome.fa.fai "$D/genome.fa.fai"
  [[ -f $D/reads.fq.gz && $(stat -c%s "$D/reads.fq.gz") -gt 100000 ]] || python3 "$ROOT/extract_reads.py" "$D/input.bam" "$D/reads.fq.gz"
  echo "[$(date +%H:%M:%S)] staged $s"
}
export -f stage_one; export ROOT SW
parallel -j6 stage_one ::: \
  SGNex_H9_directRNA_replicate2_run1 SGNex_H9_directRNA_replicate2_run2 \
  SGNex_H9_directRNA_replicate3_run1 SGNex_H9_directRNA_replicate3_run2 \
  SGNex_H9_directRNA_replicate4_run1 SGNex_H9_directRNA_replicate4_run2
echo STAGE_DONE
