#!/usr/bin/env bash
# Phase-1 Lever-1 (containment collapse) ablation on SIRV4.
# Runs the LIVE repo code (containment_shadow_drops, commit 75c01f1) via the
# production SIF's compiled runtime (krill etc.) using PYTHONPATH injection.
# Matrix: each SIRV4 sample x {denovo (no GTF), guided (full GTF)} x {off, on}.
# CPU (--no-gpu) for determinism. Idempotent (skips existing pyfin.gtf).
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
DS=$REPO/experiments/prod_validation/sirv4
OUTROOT=$REPO/experiments/prod_validation/lev1_ablation
ANN=$DS/_ref/full/annotation.gtf
mkdir -p "$OUTROOT/logs"

run_one() {  # sample mode arm
  local S="$1" mode="$2" arm="$3"
  local dst="$OUTROOT/$S/$mode/$arm"; mkdir -p "$dst"
  [ -s "$dst/pyfin.gtf" ] && { echo "[skip] $S/$mode/$arm"; return; }
  local work="$dst/work"; rm -rf "$work"; mkdir -p "$work"
  local st="$DS/$S/stage"
  local gtf_flag=(); [ "$mode" = "guided" ] && gtf_flag=(--gtf "$ANN")
  local arm_flag=(--no-containment-collapse); [ "$arm" = "on" ] && arm_flag=(--containment-collapse)
  echo "[run] $S/$mode/$arm"
  singularity exec --nv -B /SSD "$SIF" \
    env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam "$st/input.bam" --genome "$st/genome.fa" --fastq "$st/reads.fq.gz" \
    --signal "$st/signal.blow5" --signal-format slow5 --output-dir "$work" \
    --no-gpu --quant-mode m2_em "${gtf_flag[@]}" "${arm_flag[@]}" \
    > "$OUTROOT/logs/${S}__${mode}__${arm}.log" 2>&1
  if [ -s "$work/assembly.gtf" ]; then
    cp "$work/assembly.gtf" "$dst/pyfin.gtf"
    echo "[done] $S/$mode/$arm -> $(grep -c $'\ttranscript\t' "$dst/pyfin.gtf") tx"
  else
    echo "[FAIL] $S/$mode/$arm (see log)"
  fi
}

SAMPLES=$(ls -d "$DS"/SGNex_*/ 2>/dev/null | xargs -n1 basename)
for S in $SAMPLES; do
  for mode in denovo guided; do
    for arm in off on; do
      run_one "$S" "$mode" "$arm"
    done
  done
done
echo "ABLATION_DONE $(date)"
