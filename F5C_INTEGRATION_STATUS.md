# F5C Original Source Integration - Status Report

## What Was Done

### 1. Created Python Wrapper API (`fin/_align/__init__.py`)
✅ Complete high-level Python API matching `fin/_f5c` interface:
- `eventalign()` - Basic event-to-kmer alignment
- `profile_hmm_eventalign()` - Full Profile HMM with detailed output
- Auto-detection of CPU vs GPU backends
- Drop-in replacement for `fin._f5c`

### 2. Created Migration Documentation
✅ `MIGRATION_F5C.md` - Complete migration guide
✅ `migrate_examples.sh` - Automated script to update examples

### 3. Prepared C Wrapper Template
✅ `fin/_align/align_wrapper.c` - Template structure for C bindings
- Module initialization for both CPU and CUDA
- Method definitions
- Placeholder implementations

## What Needs To Be Done

### Critical: Complete C Wrapper Implementation

The file `fin/_align/align_wrapper.c` currently has placeholder implementations. You need to:

1. **Copy the implementation from `fin/_f5c/eventalign.c`**:
   - `py_eventalign()` function (lines 118-310)
   - `py_profile_hmm_eventalign()` function (lines 313-500)
   - `estimate_scalings()` helper function
   
2. **Adapt to use original f5c functions**:
   The original f5c code in `fin/_align/` uses different structures:
   - Check `fin/_align/eventalign.c` for available functions
   - May need to call `realign_read()` or similar
   - Ensure compatibility with f5c's data structures

3. **Key differences to handle**:
   ```c
   // fin/_f5c uses simplified structures:
   simple_aligned_pair_t, simple_model_t, simple_scalings_t
   
   // Original f5c uses full structures:
   AlignedPair, model_t, scalings_t, event_alignment_t
   ```

### Update setup.py

Replace `fin/_align/align_python.c` with `fin/_align/align_wrapper.c` in the Extension sources:

```python
align_extension = Extension(
    name="fin._align._align",
    sources=[
        os.path.join(ALIGN_DIR, "align_wrapper.c"),  # ← Use new wrapper
        os.path.join(ALIGN_DIR, "eventalign.c"),
        os.path.join(ALIGN_DIR, "align.c"),
        os.path.join(ALIGN_DIR, "hmm.c"),
        os.path.join(ALIGN_DIR, "model.c"),
        os.path.join(ALIGN_DIR, "nanopolish_read_db.c"),
    ],
    # ... rest of configuration
)
```

### Migration Steps

1. **Complete C wrapper**:
   ```bash
   # Edit fin/_align/align_wrapper.c
   # Copy implementations from fin/_f5c/eventalign.c
   # Adapt to original f5c function calls
   ```

2. **Update setup.py**:
   ```bash
   # Change align_python.c → align_wrapper.c in Extension sources
   ```

3. **Rebuild**:
   ```bash
   pip install -e .
   ```

4. **Test basic functionality**:
   ```bash
   python -c "from fin._align import eventalign; print('Import successful!')"
   ```

5. **Migrate examples**:
   ```bash
   chmod +x migrate_examples.sh
   ./migrate_examples.sh
   ```

6. **Test examples**:
   ```bash
   python examples/test.py
   python examples/test_profile_hmm.py
   ```

7. **Verify output matches f5c**:
   ```bash
   python examples/compare_with_f5c.py
   ```

8. **Remove fin/_f5c** (final step after all tests pass)

## API Compatibility

The new `fin._align` API is designed as a drop-in replacement:

```python
# OLD (fin._f5c)
from fin._f5c import eventalign, profile_hmm_eventalign

result = eventalign(signal, sequence, kmer_size=5)
result = profile_hmm_eventalign(signal, sequence, kmer_size=5)

# NEW (fin._align) - Same interface!
from fin._align import eventalign, profile_hmm_eventalign

result = eventalign(signal, sequence, kmer_size=5)
result = profile_hmm_eventalign(signal, sequence, kmer_size=5)

# NEW - Optional GPU selection
result = eventalign(signal, sequence, use_gpu=True)   # Force GPU
result = eventalign(signal, sequence, use_gpu=False)  # Force CPU
result = eventalign(signal, sequence)                 # Auto-detect
```

## File Structure

```
fin/_align/
├── __init__.py              ← ✅ Complete Python API
├── align_wrapper.c          ← ⏳ Needs implementation
├── align.c                  ← ✅ Original f5c source
├── align.cu                 ← ✅ Original f5c CUDA
├── eventalign.c             ← ✅ Original f5c source
├── hmm.c                    ← ✅ Original f5c source
├── model.c                  ← ✅ Original f5c source
├── model.h                  ← ✅ Model data
├── nanopolish_read_db.c     ← ✅ Original f5c source
└── ...                      ← ✅ Other f5c headers
```

## Testing Strategy

1. **Unit tests**: Ensure wrapper functions return correct data types
2. **Integration tests**: Compare output with `fin._f5c` implementation
3. **F5C comparison**: Verify results match original f5c command-line tool
4. **Performance tests**: Benchmark CPU vs GPU versions

## Benefits of This Approach

1. **Uses original f5c code** - Maximum compatibility
2. **Drop-in replacement** - No API changes needed
3. **GPU acceleration** - Automatic CUDA support
4. **Clean separation** - Original source in `_align`, simplified in `_f5c`
5. **Future-proof** - Easy to sync with f5c updates

## Next Action Required

**Implement the C wrapper in `fin/_align/align_wrapper.c`** by copying and adapting the code from `fin/_f5c/eventalign.c`. This is the only blocking task for the migration.
