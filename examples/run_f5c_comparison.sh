#!/bin/bash
# Helper script to run f5c eventalign for comparison
# 
# Usage: ./run_f5c_comparison.sh <fast5_file> <fastq_file> <reference.fasta>

set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fast5_file> <fastq_file> <reference.fasta>"
    echo ""
    echo "This script runs f5c eventalign and saves output for comparison with PyFin"
    exit 1
fi

FAST5=$1
FASTQ=$2
REFERENCE=$3

echo "=================================================="
echo "F5C EVENTALIGN COMPARISON"
echo "=================================================="
echo ""
echo "Input files:"
echo "  FAST5: $FAST5"
echo "  FASTQ: $FASTQ"
echo "  Reference: $REFERENCE"
echo ""

# Check if f5c is available
if ! command -v f5c &> /dev/null; then
    echo "❌ Error: f5c not found in PATH"
    echo ""
    echo "Install f5c from: https://github.com/hasindu2008/f5c"
    echo ""
    echo "Quick install:"
    echo "  git clone https://github.com/hasindu2008/f5c"
    echo "  cd f5c && make"
    echo "  sudo cp f5c /usr/local/bin/"
    exit 1
fi

echo "✓ f5c found: $(which f5c)"
echo ""

# Create output directory
OUTDIR="f5c_comparison_output"
mkdir -p $OUTDIR

# Index reference
echo "[1] Indexing reference..."
if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx $REFERENCE
    echo "    ✓ Reference indexed"
else
    echo "    ✓ Reference already indexed"
fi

# Run minimap2 alignment (f5c needs this)
echo ""
echo "[2] Running minimap2 alignment..."
READS_SAM="$OUTDIR/reads.sam"
minimap2 -ax map-ont -t 4 $REFERENCE $FASTQ > $READS_SAM
echo "    ✓ Alignment complete: $READS_SAM"

# Convert to BAM and sort
echo ""
echo "[3] Converting to sorted BAM..."
READS_BAM="$OUTDIR/reads.sorted.bam"
samtools view -bS $READS_SAM | samtools sort -o $READS_BAM
samtools index $READS_BAM
echo "    ✓ Sorted BAM: $READS_BAM"

# Run f5c index
echo ""
echo "[4] Running f5c index..."
f5c index -d $FAST5 $FASTQ
echo "    ✓ Index complete"

# Run f5c eventalign
echo ""
echo "[5] Running f5c eventalign..."
F5C_OUTPUT="$OUTDIR/f5c_eventalign.tsv"
f5c eventalign \
    -b $READS_BAM \
    -g $REFERENCE \
    -r $FASTQ \
    --signal-index \
    --scale-events \
    --summary $OUTDIR/f5c_summary.txt \
    > $F5C_OUTPUT

echo "    ✓ Eventalign complete: $F5C_OUTPUT"

# Show summary
echo ""
echo "=================================================="
echo "F5C OUTPUT SUMMARY"
echo "=================================================="
echo ""
wc -l $F5C_OUTPUT
head -20 $F5C_OUTPUT

echo ""
echo "=================================================="
echo "FILES CREATED"
echo "=================================================="
ls -lh $OUTDIR/
echo ""
echo "Next steps:"
echo "  1. Compare with PyFin output using compare_f5c_output.py"
echo "  2. Check $OUTDIR/f5c_summary.txt for statistics"

