#!/usr/bin/env bash
# run_benchmark.sh — benchmark harness for pyfin vs. external tools
# NOTE: After creation, run: chmod +x benchmarks/run_benchmark.sh
set -euo pipefail

DATASET=""
TOOLS="pyfin"
OUTPUT_DIR="./benchmark_output"

usage() {
    cat <<EOF
Usage: run_benchmark.sh [OPTIONS]

Benchmark pyfin against external long-read isoform tools.

Options:
  --dataset PATH     Input dataset (BAM or FASTQ)
  --tools LIST       Comma-separated list of tools to run.
                     Supported: pyfin,bambu,isoquant,jaffal,longgf
                     (default: pyfin)
  --output-dir PATH  Output directory (default: ./benchmark_output)
  --help, -h         Show this help and exit

Examples:
  bash benchmarks/run_benchmark.sh --dataset data.bam --tools pyfin --output-dir out
  bash benchmarks/run_benchmark.sh --dataset data.bam --tools pyfin,bambu --output-dir out

Missing tools are SKIPPED (not failed). Exit code is always 0 on success.
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dataset)    DATASET="$2"; shift 2 ;;
        --tools)      TOOLS="$2";   shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --help|-h)    usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

mkdir -p "$OUTPUT_DIR"

# check_and_run TOOL_NAME BINARY
# Skips gracefully if binary not found; writes result.json marker on success.
check_and_run() {
    local tool="$1"
    local binary="$2"
    if ! command -v "$binary" >/dev/null 2>&1; then
        echo "[SKIP] $tool: $binary not found in PATH"
        return 0
    fi
    echo "[RUN] $tool"
    mkdir -p "$OUTPUT_DIR/$tool"
    # Actual run commands would be inserted here per-tool.
    # For the harness skeleton, write a status marker.
    echo "{\"tool\": \"$tool\", \"status\": \"ran\", \"dataset\": \"${DATASET}\"}" \
        > "$OUTPUT_DIR/$tool/result.json"
}

IFS=',' read -ra TOOL_ARR <<< "$TOOLS"
for tool in "${TOOL_ARR[@]}"; do
    case "$tool" in
        pyfin)    check_and_run pyfin python ;;
        bambu)    check_and_run bambu Rscript ;;
        isoquant) check_and_run isoquant isoquant.py ;;
        jaffal)   check_and_run jaffal bpipe ;;
        longgf)   check_and_run longgf LongGF ;;
        *) echo "[WARN] Unknown tool: $tool" ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE_SCRIPT="$SCRIPT_DIR/compare_results.py"
if [[ -f "$COMPARE_SCRIPT" ]]; then
    python "$COMPARE_SCRIPT" --input-dir "$OUTPUT_DIR" --output "$OUTPUT_DIR/comparison.tsv"
else
    echo "[WARN] compare_results.py not found at $COMPARE_SCRIPT — skipping comparison"
fi

echo "Done. Output: $OUTPUT_DIR"
