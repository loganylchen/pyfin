#!/usr/bin/env bash
# map_drna.sh — map nanopore (d)RNA reads to genome with minimap2 (Docker)
#
# Uses the btrspg/minimap2 image (bundles minimap2 + samtools) to produce a
# sorted, indexed BAM suitable as pyfin input.
#
# Usage:
#   bash benchmarks/map_drna.sh \
#       --fastq sample.fastq.gz \
#       --ref   genome.fa \
#       --output-dir out/ \
#       [--sample-name H9_rep3] \
#       [--threads 8] \
#       [--mode drna|cdna]   (default: drna)
#       [--image btrspg/minimap2]
#
# Multiple fastq inputs (e.g. multiple replicates concatenated) can be passed
# by repeating --fastq, or by providing a single space/comma-separated list.
#
# After creation: chmod +x benchmarks/map_drna.sh
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
FASTQS=()
REF=""
OUTPUT_DIR=""
SAMPLE_NAME=""
THREADS=8
MODE="drna"
IMAGE="btrspg/minimap2"

usage() {
    cat <<EOF
Usage: map_drna.sh [OPTIONS]

Map nanopore reads to a reference genome with minimap2 inside Docker
(btrspg/minimap2), producing sorted+indexed BAM for pyfin.

Required:
  --fastq PATH         Input fastq[.gz]. Repeat for multiple files.
  --ref   PATH         Reference genome FASTA (.fa or .fa.gz)
  --output-dir PATH    Output directory (created if missing)

Optional:
  --sample-name NAME   BAM basename (default: derived from first fastq)
  --threads N          Threads for minimap2 + samtools sort (default: 8)
  --mode {drna,cdna}   Splice preset (default: drna)
                       drna: -ax splice -uf -k14   (forward-strand only)
                       cdna: -ax splice            (both strands)
  --image NAME         Docker image (default: btrspg/minimap2)
  --help, -h           Show this help

Output:
  <output-dir>/<sample-name>.bam
  <output-dir>/<sample-name>.bam.bai
  <output-dir>/<sample-name>.minimap2.log

Example:
  bash benchmarks/map_drna.sh \\
      --fastq H9_rep3_run1.fastq.gz \\
      --ref   GRCh38_plus_SIRV4.fa \\
      --output-dir bam/ \\
      --sample-name H9_rep3_run1 \\
      --threads 16
EOF
    exit 0
}

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq)
            # Accept comma- or space-separated; expand into FASTQS
            IFS=', ' read -ra _items <<< "$2"
            for f in "${_items[@]}"; do
                [[ -n "$f" ]] && FASTQS+=("$f")
            done
            shift 2 ;;
        --ref)         REF="$2";          shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2";   shift 2 ;;
        --sample-name) SAMPLE_NAME="$2";  shift 2 ;;
        --threads)     THREADS="$2";      shift 2 ;;
        --mode)        MODE="$2";         shift 2 ;;
        --image)       IMAGE="$2";        shift 2 ;;
        --help|-h)     usage ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if [[ ${#FASTQS[@]} -eq 0 ]]; then
    echo "ERROR: --fastq is required" >&2; exit 2
fi
if [[ -z "$REF" ]]; then
    echo "ERROR: --ref is required" >&2; exit 2
fi
if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --output-dir is required" >&2; exit 2
fi
if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: docker not found in PATH" >&2; exit 2
fi
if [[ ! -f "$REF" ]]; then
    echo "ERROR: reference not found: $REF" >&2; exit 2
fi
for fq in "${FASTQS[@]}"; do
    [[ -f "$fq" ]] || { echo "ERROR: fastq not found: $fq" >&2; exit 2; }
done

case "$MODE" in
    drna) MM2_PRESET=(-ax splice -uf -k14) ;;
    cdna) MM2_PRESET=(-ax splice) ;;
    *) echo "ERROR: --mode must be 'drna' or 'cdna' (got: $MODE)" >&2; exit 2 ;;
esac

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR_ABS="$(cd "$OUTPUT_DIR" && pwd)"
REF_ABS="$(cd "$(dirname "$REF")" && pwd)/$(basename "$REF")"

FASTQS_ABS=()
for fq in "${FASTQS[@]}"; do
    FASTQS_ABS+=("$(cd "$(dirname "$fq")" && pwd)/$(basename "$fq")")
done

# Derive sample name from first fastq if not given
if [[ -z "$SAMPLE_NAME" ]]; then
    base="$(basename "${FASTQS_ABS[0]}")"
    # Strip .fastq.gz / .fq.gz / .fastq / .fq
    SAMPLE_NAME="${base%.fastq.gz}"
    SAMPLE_NAME="${SAMPLE_NAME%.fq.gz}"
    SAMPLE_NAME="${SAMPLE_NAME%.fastq}"
    SAMPLE_NAME="${SAMPLE_NAME%.fq}"
fi

OUT_BAM="$OUTPUT_DIR_ABS/$SAMPLE_NAME.bam"
OUT_BAI="$OUT_BAM.bai"
OUT_LOG="$OUTPUT_DIR_ABS/$SAMPLE_NAME.minimap2.log"

# ---------------------------------------------------------------------------
# Build docker mount list (each input's parent dir, mapped 1:1)
# ---------------------------------------------------------------------------
declare -A SEEN_DIRS=()
MOUNT_ARGS=()
add_mount_ro() {
    local d="$1"
    if [[ -z "${SEEN_DIRS[$d]:-}" ]]; then
        MOUNT_ARGS+=(-v "$d:$d:ro")
        SEEN_DIRS[$d]=1
    fi
}
add_mount_rw() {
    local d="$1"
    if [[ -z "${SEEN_DIRS[$d]:-}" ]]; then
        MOUNT_ARGS+=(-v "$d:$d")
        SEEN_DIRS[$d]=1
    fi
}

add_mount_ro "$(dirname "$REF_ABS")"
for fq in "${FASTQS_ABS[@]}"; do
    add_mount_ro "$(dirname "$fq")"
done
add_mount_rw "$OUTPUT_DIR_ABS"

# ---------------------------------------------------------------------------
# Run mapping
# ---------------------------------------------------------------------------
echo "[map_drna] image:       $IMAGE"
echo "[map_drna] sample:      $SAMPLE_NAME"
echo "[map_drna] mode:        $MODE  (minimap2 ${MM2_PRESET[*]})"
echo "[map_drna] threads:     $THREADS"
echo "[map_drna] reference:   $REF_ABS"
echo "[map_drna] fastq(s):    ${FASTQS_ABS[*]}"
echo "[map_drna] output BAM:  $OUT_BAM"

# Build the in-container shell command.
# Concatenate fastqs into minimap2 via cat (handles plain + gz transparently
# since minimap2 reads gzipped fastq directly; we use 'cat' so multiple files
# stream as one input).
MM2_CMD="minimap2 ${MM2_PRESET[*]} -t $THREADS \
    --secondary=no \
    -L --MD \
    \"$REF_ABS\" ${FASTQS_ABS[*]} 2>\"$OUT_LOG\" \
  | samtools sort -@ $THREADS -m 1G -o \"$OUT_BAM\" - \
  && samtools index -@ $THREADS \"$OUT_BAM\""

# Use current uid:gid so outputs aren't root-owned.
USER_ID="$(id -u)"
GROUP_ID="$(id -g)"

docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    "${MOUNT_ARGS[@]}" \
    "$IMAGE" \
    bash -c "$MM2_CMD"

# ---------------------------------------------------------------------------
# Verify outputs
# ---------------------------------------------------------------------------
if [[ ! -s "$OUT_BAM" ]]; then
    echo "ERROR: BAM not produced or empty: $OUT_BAM" >&2; exit 3
fi
if [[ ! -s "$OUT_BAI" ]]; then
    echo "ERROR: BAM index not produced: $OUT_BAI" >&2; exit 3
fi

# Quick stats summary (also via docker so samtools is available)
echo "[map_drna] flagstat:"
docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "$OUTPUT_DIR_ABS:$OUTPUT_DIR_ABS:ro" \
    "$IMAGE" \
    samtools flagstat "$OUT_BAM" | sed 's/^/    /'

echo
echo "[map_drna] Done."
echo "  BAM:   $OUT_BAM"
echo "  BAI:   $OUT_BAI"
echo "  LOG:   $OUT_LOG"
echo
echo "Use as pyfin input:"
echo "  pyfin --bam $OUT_BAM --genome <ref.fa> --gtf <annotation.gtf> ..."
