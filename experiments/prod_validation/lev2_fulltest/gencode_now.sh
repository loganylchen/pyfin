#!/usr/bin/env bash
# Run the gencode (human) Lever-2 cells NOW (priority-ordered), in parallel with
# the main driver's heya8 phase. ON live; OFF reused (e268c9b). Idempotent ->
# the main driver skips these when it reaches gencode. GPU, 3-wide.
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
OUT=$REPO/experiments/prod_validation/lev2_fulltest
mkdir -p "$OUT/logs"
PAR=${PAR:-3}

run_pyfin() {
  local work="$1" gtf="$2"; shift 2
  local gflag=(); [ "$gtf" != "-" ] && gflag=(--gtf "$gtf")
  rm -rf "$work"; mkdir -p "$work"
  singularity exec --nv -B /SSD -B /autofs/mnemosyne3_SSD "$SIF" \
    env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam "$ST/input.bam" --genome "$ST/genome.fa" --fastq "$ST/reads.fq.gz" \
    --signal "$ST/signal.blow5" --signal-format slow5 --output-dir "$work" \
    --gpu --quant-mode m2_em "${gflag[@]}" "$@" >/dev/null 2>&1
}
run_cell() {
  IFS=$'\t' read -r ds samp ratio gtf off gpu run_off <<< "$1"
  local dst="$OUT/$ds/$samp/$ratio"; mkdir -p "$dst"
  local ST="$REPO/experiments/prod_validation/gencode/$samp/stage"
  if [ ! -s "$dst/on.gtf" ]; then
    run_pyfin "$dst/work_on" "$gtf" --novel-junction-min-reads 2
    [ -s "$dst/work_on/assembly.gtf" ] && cp "$dst/work_on/assembly.gtf" "$dst/on.gtf"
  fi
  if [ ! -s "$dst/off.gtf" ] && [ "$off" != "-" ] && [ -s "$off" ]; then cp "$off" "$dst/off.gtf"; fi
  echo "[gc-done $(date +%H:%M:%S)] $samp/$ratio on=$(grep -c $'\ttranscript\t' "$dst/on.gtf" 2>/dev/null||echo X) off=$(grep -c $'\ttranscript\t' "$dst/off.gtf" 2>/dev/null||echo X)"
}
export -f run_cell run_pyfin; export OUT SIF REPO

while IFS= read -r line; do
  run_cell "$line" &
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 8; done
done < "$OUT/gc_priority.tsv"
wait
echo "GENCODE_NOW_DONE $(date)"
