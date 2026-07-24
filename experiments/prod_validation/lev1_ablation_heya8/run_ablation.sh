#!/usr/bin/env bash
# Phase-1 Lever-1 (containment) ablation on heya8 (dense-locus dRNA spike-in).
# LIVE repo code via SIF runtime (PYTHONPATH inject). 5 samples x {denovo,guided}
# x {off,on}. PARALLEL across samples (GPU idle on this workload -> CPU krill is
# the bottleneck; 5 workers x ~16 threads fits the 128-core box). Idempotent.
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
DS=$REPO/experiments/wobble_heya8/matrix
ANN=$DS/_ref/full/annotation.gtf
OUTROOT=$REPO/experiments/prod_validation/lev1_ablation_heya8
mkdir -p "$OUTROOT/logs"

run_one() {
  local S="$1" mode="$2" arm="$3"
  local dst="$OUTROOT/$S/$mode/$arm"; mkdir -p "$dst"
  [ -s "$dst/pyfin.gtf" ] && { echo "[skip] $S/$mode/$arm"; return; }
  local work="$dst/work"; rm -rf "$work"; mkdir -p "$work"
  local st="$DS/$S/stage"
  local gtf_flag=(); [ "$mode" = "guided" ] && gtf_flag=(--gtf "$ANN")
  local arm_flag=(--no-containment-collapse); [ "$arm" = "on" ] && arm_flag=(--containment-collapse)
  echo "[run $(date +%H:%M:%S)] $S/$mode/$arm"
  singularity exec --nv -B /SSD "$SIF" \
    env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam "$st/input.bam" --genome "$st/genome.fa" --fastq "$st/reads.fq.gz" \
    --signal "$st/signal.blow5" --signal-format slow5 --output-dir "$work" \
    --no-gpu --quant-mode m2_em "${gtf_flag[@]}" "${arm_flag[@]}" \
    > "$OUTROOT/logs/${S}__${mode}__${arm}.log" 2>&1
  if [ -s "$work/assembly.gtf" ]; then
    cp "$work/assembly.gtf" "$dst/pyfin.gtf"
    echo "[done $(date +%H:%M:%S)] $S/$mode/$arm -> $(grep -c $'\ttranscript\t' "$dst/pyfin.gtf") tx"
  else
    echo "[FAIL] $S/$mode/$arm"
  fi
}

worker() {  # one sample, its 4 cells serial
  local S="$1"
  for mode in denovo guided; do
    for arm in off on; do run_one "$S" "$mode" "$arm"; done
  done
  echo "[WORKER DONE $(date +%H:%M:%S)] $S"
}

for S in $(ls -d "$DS"/SGNex_HEYA8*/ 2>/dev/null | xargs -n1 basename); do
  worker "$S" &
done
wait
echo "ABLATION_DONE $(date)"
