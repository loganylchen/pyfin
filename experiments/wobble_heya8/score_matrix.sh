#!/usr/bin/env bash
# Score every (sample, ratio, tool, config) cell with gffcompare vs truth=full annotation.
# Tools = pyfin@bp10, pyfin@bp20, + 8 baseline tools (GTFs symlinked from NFS).
# Idempotent: skips when <prefix>.stats exists.
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
RATIOS=(full p00 c_jitter10bp c_skip10 c_spurious5 c_merge5 c_flip5 c_ir10)
BASELINE_TOOLS=(bambu espresso flair isoquant isotools lafite stringtie3 talon)
NFS_BASE=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/heya8_dual_full_sweep
TRUTH="matrix/_ref/full/annotation.gtf"
GFFC_IMG="quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
MNT="$(readlink -f "$ROOT")"
NFS_MNT="/autofs/mnemosyne3_SSD"

gffc() {
  local outdir=$1 prefix=$2 query=$3
  mkdir -p "$outdir"
  local abs_outdir; abs_outdir="$(readlink -f "$outdir")"
  local abs_truth;  abs_truth="$(readlink -f "$TRUTH")"
  local abs_query;  abs_query="$(readlink -f "$query")"
  [[ -f "$abs_outdir/${prefix}.stats" ]] && return 0
  docker run --rm -u "$(id -u):$(id -g)" \
    -v "$MNT":"$MNT" -v "$NFS_MNT":"$NFS_MNT":ro \
    -w "$abs_outdir" "$GFFC_IMG" \
    gffcompare -r "$abs_truth" -o "$prefix" "$abs_query" >/dev/null 2>&1 || return 1
}

# Stage baseline tool GTFs as symlinks (don't copy)
for s in "${SAMPLES[@]}"; do
  for r in "${RATIOS[@]}"; do
    for t in "${BASELINE_TOOLS[@]}"; do
      src="$NFS_BASE/$s/$r/assembly/$t.gtf"
      dst="matrix/_other_tools/$s/$r/$t.gtf"
      mkdir -p "$(dirname "$dst")"
      [[ -e "$src" && ! -e "$dst" ]] && ln -sf "$src" "$dst" || true
    done
  done
done

# Score
total=0; scored=0; skipped=0; missing=0
for s in "${SAMPLES[@]}"; do
  for r in "${RATIOS[@]}"; do
    SDIR="matrix/$s/$r/scoring"; mkdir -p "$SDIR"
    # pyfin configs: bp10 baseline, bp20 wobble-fix, bp20cas (+cassette), v2 (+displace)
    for c in bp10 bp20cas v3 prodfinal; do
      q="matrix/$s/$r/out_$c/assembly.gtf"
      total=$((total+1))
      if [[ -f "$q" ]]; then
        if [[ -f "$SDIR/pyfin_${c}.stats" ]]; then skipped=$((skipped+1))
        else gffc "$SDIR" "pyfin_${c}" "$(readlink -f "$q")" && scored=$((scored+1)) || true
        fi
      else missing=$((missing+1))
      fi
    done
    # baselines
    for t in "${BASELINE_TOOLS[@]}"; do
      q="matrix/_other_tools/$s/$r/$t.gtf"
      total=$((total+1))
      if [[ -e "$q" ]]; then
        if [[ -f "$SDIR/${t}.stats" ]]; then skipped=$((skipped+1))
        else gffc "$SDIR" "$t" "$q" && scored=$((scored+1)) || true
        fi
      else missing=$((missing+1))
      fi
    done
  done
done
echo "scoring done: total=$total scored=$scored skipped=$skipped missing=$missing"
