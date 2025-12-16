# Examples Fixed - Summary

## Issues Found and Fixed

### 1. Import Path Issues

**Problem**: Many examples used incorrect import paths:
```python
from fin import eventalign, detect_events  # WRONG - doesn't exist
```

**Fixed**: Updated to correct module structure:
```python
from fin._f5c._event import detect_events
from fin._f5c._eventalign import eventalign
```

**Files Fixed**:
- test_eventalign.py
- debug_eventalign.py
- test_trimming_fix.py
- compare_with_f5c.py
- test_event.py
- raw_signal_alignment_example.py

### 2. Missing Error Handling

**Problem**: Examples would crash with unclear errors if extensions not built

**Fixed**: Added try/except blocks with helpful error messages:
```python
try:
    from fin._f5c._eventalign import eventalign
except ImportError:
    print("Error: eventalign extension not available")
    print("Build with: python setup.py build_ext --inplace")
    sys.exit(1)
```

### 3. Deprecated Examples

**Problem**: Some examples referenced old API that no longer exists

**Fixed**:
- eventalign_example.py - Updated to redirect users to working examples
- Added deprecation notices

### 4. Missing Documentation

**Problem**: No central guide for examples

**Fixed**:
- Created examples/README.md with full documentation
- Categorized examples by functionality
- Added troubleshooting section
- Added workflow guide

### 5. No Testing Infrastructure

**Problem**: No way to verify all examples work

**Fixed**:
- Created test_all_examples.py script
- Automatically tests all simple examples
- Lists advanced/skipped examples
- Provides summary of results

## Files Modified

### Updated Examples (7 files)
1. test_eventalign.py - Fixed imports, removed pod5 dependency
2. debug_eventalign.py - Fixed imports
3. test_trimming_fix.py - Fixed imports
4. compare_with_f5c.py - Fixed imports
5. test_event.py - Fixed imports
6. raw_signal_alignment_example.py - Fixed imports
7. eventalign_example.py - Added deprecation notice

### New Files (2 files)
1. test_all_examples.py - Test harness for all examples
2. README.md - Comprehensive examples documentation

## Status of All Examples

### ✓ Working (8 files)
- test_event.py
- test_trimming_fix.py
- test_profile_hmm.py
- raw_signal_alignment_example.py
- debug_eventalign.py
- test_adapter_trimming.py
- dtw_nanopore_example.py
- test_package.py

### ⚠ Needs Special Setup (3 files)
- test_eventalign.py - Needs matplotlib
- compare_with_f5c.py - Needs f5c installed
- eventalign_example.py - Deprecated

### ⊘ Requires Data/Hardware (6 files)
- compare_profile_hmm_with_f5c.py - Needs f5c output
- test_cuda_dtw.py - Needs CUDA
- benchmark_dtw_comparison.py - Needs data
- benchmark_dtw_pairwise.py - Needs data
- dtw_pairwise_example.py - Needs data
- interval_workflow_example.py - Needs BAM files

## Testing

Run the test suite:
```bash
python examples/test_all_examples.py
```

Test specific example:
```bash
python examples/test_profile_hmm.py
```

## Next Steps for Users

1. Build extensions:
   ```bash
   cd /home/logan/Projects/pyfin
   python setup.py build_ext --inplace
   ```

2. Run test suite:
   ```bash
   python examples/test_all_examples.py
   ```

3. Try recommended example:
   ```bash
   python examples/test_profile_hmm.py
   ```

4. Read examples/README.md for detailed documentation

## Breaking Changes

Users upgrading from old code need to update imports:

**Before:**
```python
from fin import detect_events, eventalign
```

**After:**
```python
from fin._f5c._event import detect_events
from fin._f5c._eventalign import eventalign
```

Or use the new Profile HMM API:
```python
from fin._f5c._eventalign import profile_hmm_eventalign
```
