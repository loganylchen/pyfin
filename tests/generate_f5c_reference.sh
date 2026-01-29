#!/bin/bash
# Generate f5c reference eventalign outputs for test data
#
# Prerequisites:
#   - f5c installed (https://github.com/hasindu2008/f5c)
#   - samtools installed
#   - minimap2 installed (for indexing if needed)
#
# Usage:
#   ./generate_f5c_reference.sh
#
# This script generates f5c eventalign TSV outputs that serve as the
# reference ("ground truth") for comparing PyFIN's implementation.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TEST_DATA_DIR="$PROJECT_ROOT/examples/test_data"

echo "=============================================="
echo "Generating f5c Reference Outputs"
echo "=============================================="
echo ""
echo "Test data directory: $TEST_DATA_DIR"
echo ""

# Check if f5c is installed
if ! command -v f5c &> /dev/null; then
    echo "ERROR: f5c not found in PATH"
    echo "Please install f5c: https://github.com/hasindu2008/f5c"
    exit 1
fi

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in PATH"
    exit 1
fi

cd "$TEST_DATA_DIR"

# Function to generate f5c eventalign for a dataset
generate_eventalign() {
    local PREFIX=$1
    local MODEL_TYPE=$2  # rna002 or rna004
    
    local POD5_FILE="${PREFIX}.test.pod5"
    local FASTQ_FILE="${PREFIX}.test.fastq.gz"
    local BAM_FILE="${PREFIX}.test.bam"
    local OUTPUT_FILE="${PREFIX}.test.eventalign.tsv.gz"
    
    echo "Processing ${PREFIX} (${MODEL_TYPE})..."
    
    # Check if input files exist
    if [[ ! -f "$POD5_FILE" ]]; then
        echo "  WARNING: $POD5_FILE not found, skipping"
        return
    fi
    
    if [[ ! -f "$FASTQ_FILE" ]]; then
        echo "  WARNING: $FASTQ_FILE not found, skipping"
        return
    fi
    
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "  WARNING: $BAM_FILE not found, skipping"
        return
    fi
    
    if [[ ! -f "reference.fa" ]]; then
        echo "  WARNING: reference.fa not found, skipping"
        return
    fi
    
    # Index reference if needed
    if [[ ! -f "reference.fa.fai" ]]; then
        echo "  Indexing reference.fa..."
        samtools faidx reference.fa
    fi
    
    # Index BAM if needed
    if [[ ! -f "${BAM_FILE}.bai" ]]; then
        echo "  Indexing ${BAM_FILE}..."
        samtools index "$BAM_FILE"
    fi
    
    # Create f5c index if needed
    local INDEX_FILE="${PREFIX}.test.index"
    if [[ ! -f "${INDEX_FILE}.readdb" ]]; then
        echo "  Creating f5c index..."
        f5c index \
            --slow5 "$POD5_FILE" \
            "$FASTQ_FILE"
    fi
    
    # Run f5c eventalign
    echo "  Running f5c eventalign..."
    
    # Determine model flag based on type
    local MODEL_FLAG=""
    if [[ "$MODEL_TYPE" == "rna004" ]]; then
        MODEL_FLAG="--rna"  # f5c auto-detects RNA004 from model
    else
        MODEL_FLAG="--rna"  # RNA002 is default RNA model
    fi
    
    f5c eventalign \
        $MODEL_FLAG \
        --reads "$FASTQ_FILE" \
        --bam "$BAM_FILE" \
        --genome reference.fa \
        --slow5 "$POD5_FILE" \
        --signal-index \
        --scale-events \
        --print-read-names \
        --samples \
        2>/dev/null | gzip > "$OUTPUT_FILE"
    
    # Count output lines
    local N_LINES=$(zcat "$OUTPUT_FILE" | wc -l)
    echo "  Generated $N_LINES event alignments -> $OUTPUT_FILE"
}

# Generate for RNA004
if [[ -f "RNA004.test.pod5" ]]; then
    generate_eventalign "RNA004" "rna004"
else
    echo "RNA004 test data not found, skipping"
fi

echo ""

# Generate for RNA002
if [[ -f "RNA002.test.pod5" ]]; then
    generate_eventalign "RNA002" "rna002"
else
    echo "RNA002 test data not found, skipping"
fi

echo ""
echo "=============================================="
echo "Done!"
echo "=============================================="
echo ""
echo "Generated reference outputs:"
ls -la *.eventalign.tsv.gz 2>/dev/null || echo "  (none)"
