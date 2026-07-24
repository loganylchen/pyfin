#!/usr/bin/env bash
# Lever-3 (mono-exon gate) ablation. ON arm only; OFF reused from prior runs
# (all levers default-off -> byte-identical OFF baseline).
# ON flags: --drop-mono-exon-novel --min-mono-exon-reads 3 --min-mono-exon-length 200
# Cells: SIRV4 (denovo+guided, CPU), heya8 (denovo, CPU), gencode (nogtf, GPU).
set -uo pipefail
REPO=/SSD/logan/dev/pyfin
SIF=$REPO/experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif
OUT=$REPO/experiments/prod_validation/lev2_ablation
mkdir -p "$OUT/logs"
MONO=(--novel-junction-min-reads 2 --novel-junction-reads-tol 2)
MAXPAR=${MAXPAR:-4}

# cell = dataset|sample|mode|stage_dir|gtf(or -)|off_gtf|gpu(0/1)
CELLS=()
SIRV=$REPO/experiments/prod_validation/sirv4
for S in $(ls -d "$SIRV"/SGNex_*/ | xargs -n1 basename); do
  st="$SIRV/$S/stage"
  CELLS+=("sirv|$S|denovo|$st|-|$REPO/experiments/prod_validation/lev1_ablation/$S/denovo/off/pyfin.gtf|0")
  CELLS+=("sirv|$S|guided|$st|$SIRV/_ref/full/annotation.gtf|$REPO/experiments/prod_validation/lev1_ablation/$S/guided/off/pyfin.gtf|0")
done
HEYA=$REPO/experiments/wobble_heya8/matrix
for S in $(ls -d "$HEYA"/SGNex_HEYA8*/ | xargs -n1 basename); do
  st="$HEYA/$S/stage"
  CELLS+=("heya8|$S|denovo|$st|-|$REPO/experiments/prod_validation/lev1_ablation_heya8/$S/denovo/off/pyfin.gtf|0")
done
GENC=$REPO/experiments/prod_validation/gencode
for S in SGNex_H9_directRNA_replicate2_run2 SGNex_H9_directRNA_replicate4_run2; do
  st="$GENC/$S/stage"
  CELLS+=("gencode|$S|nogtf|$st|-|$GENC/$S/nogtf/pyfin.gtf|1")
done

run_cell() {
  IFS='|' read -r ds S mode st gtf off gpu <<< "$1"
  local dst="$OUT/$ds/$S/$mode"; mkdir -p "$dst"
  [ -s "$dst/on.gtf" ] && { echo "[skip] $ds/$S/$mode"; return; }
  local work="$dst/work_on"; rm -rf "$work"; mkdir -p "$work"
  local gtf_flag=(); [ "$gtf" != "-" ] && gtf_flag=(--gtf "$gtf")
  local gpu_flag=(--no-gpu); [ "$gpu" = "1" ] && gpu_flag=(--gpu)
  echo "[run $(date +%H:%M:%S)] $ds/$S/$mode"
  singularity exec --nv -B /SSD -B /autofs/mnemosyne3_SSD "$SIF" \
    env PYTHONPATH=$REPO /usr/bin/python3.10 -m fin.cli \
    --bam "$st/input.bam" --genome "$st/genome.fa" --fastq "$st/reads.fq.gz" \
    --signal "$st/signal.blow5" --signal-format slow5 --output-dir "$work" \
    "${gpu_flag[@]}" --quant-mode m2_em "${gtf_flag[@]}" "${MONO[@]}" \
    > "$OUT/logs/${ds}__${S}__${mode}.log" 2>&1
  if [ -s "$work/assembly.gtf" ]; then
    cp "$work/assembly.gtf" "$dst/on.gtf"
    [ -s "$off" ] && cp "$off" "$dst/off.gtf"
    echo "[done $(date +%H:%M:%S)] $ds/$S/$mode -> on=$(grep -c $'\ttranscript\t' "$dst/on.gtf") off=$(grep -c $'\ttranscript\t' "$dst/off.gtf" 2>/dev/null||echo NA)"
  else
    echo "[FAIL] $ds/$S/$mode"
  fi
}
export -f run_cell; export OUT SIF REPO; export MONO_STR="${MONO[*]}"

for c in "${CELLS[@]}"; do
  run_cell "$c" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXPAR" ]; do sleep 10; done
done
wait
echo "LEV2_DONE $(date)"
