# Migration Guide: fin/_f5c → fin/_align

## Overview

The package is being restructured to use the original f5c source code directly from `fin/_align` instead of the simplified implementation in `fin/_f5c`.

## Status

### What's Complete
- ✅ `fin/_align/__init__.py` - Python wrapper with GPU/CPU backend selection
- ✅ Original f5c C/C++ source code in `fin/_align/`
- ✅ Build configuration in `setup.py` for both CPU and CUDA versions

### What's Needed
- ⏳ Complete C wrapper implementation in `fin/_align/align_wrapper.c`
  - Currently has placeholder implementations
  - Needs to call actual f5c functions from `eventalign.c`
  - Reference: `fin/_f5c/eventalign.c` for wrapper pattern

### Migration Path

1. **Complete the C wrapper** (`fin/_align/align_wrapper.c`):
   - Copy implementation from `fin/_f5c/eventalign.c`
   - Adapt to call original f5c functions
   - Handle both `eventalign` and `profile_hmm_eventalign`

2. **Update setup.py** to use new wrapper:
   ```python
   align_extension = Extension(
       name="fin._align._align",
       sources=[
           os.path.join(ALIGN_DIR, "align_wrapper.c"),  # New wrapper
           os.path.join(ALIGN_DIR, "eventalign.c"),     # Original f5c
           os.path.join(ALIGN_DIR, "align.c"),
           os.path.join(ALIGN_DIR, "hmm.c"),
           os.path.join(ALIGN_DIR, "model.c"),
           # ... other sources
       ],
       # ... configuration
   )
   ```

3. **Update all examples** to import from `fin._align`:
   ```python
   # OLD
   from fin._f5c import eventalign, profile_hmm_eventalign
   
   # NEW  
   from fin._align import eventalign, profile_hmm_eventalign
   ```

4. **Test and verify**:
   - Build with: `pip install -e .`
   - Run examples to ensure compatibility
   - Compare output with original f5c

5. **Remove fin/_f5c** (final step):
   - After all tests pass
   - Update documentation
   - Remove from setup.py

## Examples to Update

The following files import from `fin._f5c` and need updating:

```
examples/test_profile_hmm.py
examples/compare_profile_hmm_with_f5c.py  
examples/test.py
examples/diagnose_eventalign_diff.py
examples/debug_eventalign.py
examples/raw_signal_alignment_example.py
examples/test_event.py
examples/benchmark_cpu_vs_gpu.py
examples/region_transcript_analysis_workflow.py
examples/test_trimming_fix.py
examples/test_eventalign.py
examples/compare_with_f5c.py
examples/eventalign_example.py
```

## Simple Find & Replace

For all Python files in `examples/`:

```bash
# Replace module imports
sed -i 's/from fin\._f5c/from fin._align/g' examples/*.py
sed -i 's/fin\._f5c\./fin._align./g' examples/*.py
```

## Key Differences

### API Compatibility
The `fin._align` API is designed to be a drop-in replacement:
- Same function names: `eventalign()`, `profile_hmm_eventalign()`
- Same arguments and return values
- Additional `use_gpu` parameter for backend selection

### Backend Selection
```python
# Auto-detect (uses CUDA if available)
result = eventalign(signal, sequence)

# Force CPU
result = eventalign(signal, sequence, use_gpu=False)

# Force GPU
result = eventalign(signal, sequence, use_gpu=True)
```

## Testing

After migration, verify with:
```bash
# Build
pip install -e .

# Test basic functionality
python examples/test.py

# Test Profile HMM
python examples/test_profile_hmm.py

# Compare with f5c output
python examples/compare_with_f5c.py
```

## Notes

- The original f5c code in `fin/_align` is C++ (uses `<map>`, `<vector>`, etc.)
- Build system now compiles `.c` files as C++ using g++
- Both CPU and CUDA versions share same source files
- CUDA version enabled with `-DCUDA_ENABLED` flag
