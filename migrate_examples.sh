#!/bin/bash
# Script to migrate examples from fin._f5c to fin._align

echo "Migrating examples from fin._f5c to fin._align..."

# List of files to update
FILES=(
    "examples/test_profile_hmm.py"
    "examples/compare_profile_hmm_with_f5c.py"
    "examples/test.py"
    "examples/diagnose_eventalign_diff.py"
    "examples/debug_eventalign.py"
    "examples/raw_signal_alignment_example.py"
    "examples/test_event.py"
    "examples/benchmark_cpu_vs_gpu.py"
    "examples/region_transcript_analysis_workflow.py"
    "examples/test_trimming_fix.py"
    "examples/test_eventalign.py"
    "examples/compare_with_f5c.py"
    "examples/eventalign_example.py"
)

for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "Updating $file..."
        # Replace imports
        sed -i 's/from fin\._f5c/from fin._align/g' "$file"
        sed -i 's/fin\._f5c\./fin._align./g' "$file"
    else
        echo "Warning: $file not found"
    fi
done

echo "Migration complete!"
echo ""
echo "Next steps:"
echo "1. Complete the C wrapper in fin/_align/align_wrapper.c"  
echo "2. Update setup.py to use the new wrapper"
echo "3. Rebuild: pip install -e ."
echo "4. Test: python examples/test.py"
