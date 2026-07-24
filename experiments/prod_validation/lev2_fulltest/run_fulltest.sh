#!/usr/bin/env bash
# Full Lever-2 (--novel-junction-min-reads 2) test across ALL datasets/ratios.
# ON arm always run live; OFF reused (sirv/gencode=e268c9b) or run live (heya8).
# Processed dataset-by-dataset (sirv -> heya8 -> gencode) so each finishes and
# can be scored progressively. CPU for sirv/heya8, GPU for gencode. Idempotent.
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
OUT=$REPO/experiments/prod_validation/lev2_fulltest
CELLS=$OUT/cells.tsv
mkdir -p "$OUT/logs"
CPUPAR=${CPUPAR:-6}
GPUPAR=${GPUPAR:-3}

stage_dir() { case "$1" in
  sirv) echo "$REPO/experiments/prod_validation/sirv4/$2/stage";;
  gencode) echo "$REPO/experiments/prod_validation/gencode/$2/stage";;
  heya8) echo "$REPO/experiments/wobble_heya8/matrix/$2/stage";;
esac; }

run_pyfin() {  # outdir gpu gtf extra... ; uses stage from globals ST
  local work="$1" gpu="$2" gtf="$3"; shift 3
  local gflag=(); [ "$gtf" != "-" ] && gflag=(--gtf "$gtf")
  local gpuflag=(--no-gpu); [ "$gpu" = "1" ] && gpuflag=(--gpu)
  rm -rf "$work"; mkdir -p "$work"
  singularity exec --nv -B /SSD -B /autofs/mnemosyne3_SSD "$SIF" \
    env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam "$ST/input.bam" --genome "$ST/genome.fa" --fastq "$ST/reads.fq.gz" \
    --signal "$ST/signal.blow5" --signal-format slow5 --output-dir "$work" \
    "${gpuflag[@]}" --quant-mode m2_em "${gflag[@]}" "$@" >/dev/null 2>&1
}

run_cell() {
  IFS=$'\t' read -r ds samp ratio gtf off gpu run_off <<< "$1"
  local dst="$OUT/$ds/$samp/$ratio"; mkdir -p "$dst"
  local ST; ST=$(stage_dir "$ds" "$samp")
  local lg="$OUT/logs/${ds}__${samp}__${ratio}"
  # ON arm
  if [ ! -s "$dst/on.gtf" ]; then
    run_pyfin "$dst/work_on" "$gpu" "$gtf" --novel-junction-min-reads 2 2>"$lg.on.err"
    [ -s "$dst/work_on/assembly.gtf" ] && cp "$dst/work_on/assembly.gtf" "$dst/on.gtf"
  fi
  # OFF arm
  if [ ! -s "$dst/off.gtf" ]; then
    if [ "$run_off" = "1" ]; then
      run_pyfin "$dst/work_off" "$gpu" "$gtf" 2>"$lg.off.err"
      [ -s "$dst/work_off/assembly.gtf" ] && cp "$dst/work_off/assembly.gtf" "$dst/off.gtf"
    elif [ "$off" != "-" ] && [ -s "$off" ]; then
      cp "$off" "$dst/off.gtf"
    fi
  fi
  echo "[done $(date +%H:%M:%S)] $ds/$samp/$ratio on=$(grep -c $'\ttranscript\t' "$dst/on.gtf" 2>/dev/null||echo X) off=$(grep -c $'\ttranscript\t' "$dst/off.gtf" 2>/dev/null||echo X)"
}
export -f run_cell run_pyfin stage_dir; export OUT SIF REPO

for ds in sirv heya8 gencode; do
  par=$CPUPAR; [ "$ds" = "gencode" ] && par=$GPUPAR
  echo "=== dataset $ds (par=$par) $(date) ==="
  while IFS= read -r line; do
    run_cell "$line" &
    while [ "$(jobs -rp | wc -l)" -ge "$par" ]; do sleep 8; done
  done < <(grep -P "^$ds\t" "$CELLS")
  wait
  echo "DATASET_DONE $ds $(date)"
done
echo "FULLTEST_DONE $(date)"
