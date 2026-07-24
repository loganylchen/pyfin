#!/usr/bin/env bash
# Lever-1 ablation on HUMAN GENCODE (real transcriptome).
# OFF baseline already exists: gencode/<S>/{full,nogtf}/pyfin.gtf were produced by
# the e268c9b SIF == the direct parent of the Lever-1 commit, so they are
# byte-identical to live-code OFF. Here we run only the ON arm (live code,
# --containment-collapse) for guided (full GTF) and de novo (no GTF).
# GPU; samples run with bounded parallelism. Idempotent.
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
GC=$REPO/experiments/prod_validation/gencode
ANN=$GC/_ref/full/annotation.gtf
OUT=$REPO/experiments/prod_validation/lev1_ablation_gencode
mkdir -p "$OUT/logs"
MAXPAR=${MAXPAR:-3}

# Samples to run (default: all 5 that have OFF baselines).
SAMPLES=(${SAMPLES:-$(ls -d "$GC"/SGNex_*/stage 2>/dev/null | xargs -n1 dirname | xargs -n1 basename)})

run_on() {  # sample mode(full|nogtf)
  local S="$1" mode="$2"
  local dst="$OUT/$S/$mode"; mkdir -p "$dst"
  [ -s "$dst/on.gtf" ] && { echo "[skip] $S/$mode"; return; }
  local work="$dst/work_on"; rm -rf "$work"; mkdir -p "$work"
  local st="$GC/$S/stage"
  local gtf_flag=(); [ "$mode" = "full" ] && gtf_flag=(--gtf "$ANN")
  echo "[run $(date +%H:%M:%S)] $S/$mode/on"
  singularity exec --nv -B /SSD -B /autofs/mnemosyne3_SSD "$SIF" \
    env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam "$st/input.bam" --genome "$st/genome.fa" --fastq "$st/reads.fq.gz" \
    --signal "$st/signal.blow5" --signal-format slow5 --output-dir "$work" \
    --gpu --quant-mode m2_em "${gtf_flag[@]}" --containment-collapse \
    > "$OUT/logs/${S}__${mode}__on.log" 2>&1
  if [ -s "$work/assembly.gtf" ]; then
    cp "$work/assembly.gtf" "$dst/on.gtf"
    # stage the existing OFF baseline next to it for scoring convenience
    [ -s "$GC/$S/$mode/pyfin.gtf" ] && cp "$GC/$S/$mode/pyfin.gtf" "$dst/off.gtf"
    echo "[done $(date +%H:%M:%S)] $S/$mode/on -> $(grep -c $'\ttranscript\t' "$dst/on.gtf") tx (off=$(grep -c $'\ttranscript\t' "$dst/off.gtf" 2>/dev/null || echo NA))"
  else
    echo "[FAIL] $S/$mode/on"
  fi
}

n=0
for S in "${SAMPLES[@]}"; do
  for mode in full nogtf; do
    run_on "$S" "$mode" &
    n=$((n+1))
    while [ "$(jobs -rp | wc -l)" -ge "$MAXPAR" ]; do sleep 15; done
  done
done
wait
echo "GENCODE_ON_DONE $(date)"
