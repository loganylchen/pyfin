# Comparison: PyFin vs F5C Eventalign Implementation

## Critical Differences Found

### 1. **Algorithm Type**

**F5C (Original):**
- Uses **Profile HMM with Viterbi** algorithm
- Implements `profile_hmm_align()` function
- Full forward-backward HMM with 3 states per kmer
- Output: `std::vector<event_alignment_t>` with detailed alignment states

**PyFin (Current):**
- Uses **Adaptive Banded Event Alignment** (ABEA)
- Implements `align_with_flanking_cpu()` function
- Simplified 3-state DP with banding
- Output: `simple_aligned_pair_t` array (event_idx, kmer_idx pairs)

### 2. **File Structure Comparison**

**F5C:**
```
src/eventalign.c     - Main eventalign logic with HMM
src/hmm.c            - Profile HMM implementation
src/align.c          - Adaptive banded alignment (different algorithm!)
src/align.cu         - GPU version of ABEA
```

**PyFin:**
```
fin/_f5c/eventalign.c  - Python wrapper only
fin/_f5c/align.c       - 3-state DP alignment
fin/_f5c/align.cu      - GPU version
```

### 3. **Key Code Differences**

#### F5C Eventalign Flow:
```c
// In eventalign.c: realign_read()
EventAlignmentParameters params;
params.et = events;
params.model = model;
params.kmer_size = kmer_size;

// Calls HMM alignment
std::vector<event_alignment_t> alignment = align_read_to_ref(params, ref);

// align_read_to_ref() uses:
std::vector<HMMAlignmentState> event_alignment = profile_hmm_align(
    fwd_subseq, rc_subseq, event, scaling, model,
    events_per_base, strand, rc, k, e_start, e_end, event_stride
);

// profile_hmm_align() uses:
profile_hmm_fill_generic_r9(...)  // Full HMM forward pass
// Then backtracks to get alignment states
```

#### PyFin Current Flow:
```c
// In eventalign.c: py_eventalign()
// Detects events
event_table et = getevents_simple(nsample, rawptr, is_rna);

// Estimates scaling
simple_scalings_t scaling = estimate_scalings(sequence, seq_len, model, kmer_size, et);

// Calls simplified alignment
int n_pairs = align_with_flanking_cpu(&aligned_pairs, sequence, seq_len, et, model, kmer_size, scaling);

// align_with_flanking_cpu() is NOT the same as f5c's eventalign!
// It's more like f5c's align() function (ABEA)
```

### 4. **Alignment Output Format**

**F5C event_alignment_t:**
```c
struct event_alignment_t {
    uint32_t ref_position;      // Reference position
    char ref_kmer[MAX_KMER_SIZE+1];
    int event_idx;              // Event index
    char hmm_state;             // 'M', 'K', or 'B'
    bool rc;                    // Reverse complement flag
    char model_kmer[MAX_KMER_SIZE+1];
    // Plus more fields...
};
```

**PyFin simple_aligned_pair_t:**
```c
typedef struct {
    int32_t ref_pos;    // Kmer index
    int32_t read_pos;   // Event index
    // That's it - much simpler!
} simple_aligned_pair_t;
```

### 5. **HMM Implementation**

**F5C:**
```c
// Has full 3-state HMM per kmer:
enum ProfileStateR9 {
    PSR9_KMER_SKIP = 0,   // Kmer with no event
    PSR9_BAD_EVENT,       // Bad event to ignore
    PSR9_MATCH,           // Event matches kmer
    PSR9_NUM_STATES = 3
};

// Complex transition model:
BlockTransitions transitions[n_kmers];
// Each has 10 transition probabilities:
// lp_mm_self, lp_mm_next, lp_mb, lp_mk,
// lp_bb, lp_bk, lp_bm_next, lp_bm_self,
// lp_kk, lp_km

// DP table size: [n_events + 1] x [n_kmers * 3 + 6]
// (3 states per kmer + terminal states)
```

**PyFin:**
```c
// Similar structure but implementation differs
// DP table: [n_events + 1] x [n_kmers * 3]
// But the filling logic is different from f5c's profile_hmm_fill_generic_r9()
```

### 6. **What You're Actually Running**

Your current code is running **f5c's "align" function** (Adaptive Banded Event Alignment), NOT f5c's "eventalign" function (Profile HMM alignment).

In f5c:
- `align()` in `src/align.c` = Fast banded DP for getting rough alignment
- `profile_hmm_align()` in `src/eventalign.c` = Detailed HMM for exact alignment

### 7. **To Reproduce F5C Eventalign Results**

You need to implement the **full Profile HMM** from `src/eventalign.c` and `src/hmm.c`, not just the banded alignment.

Key functions missing from your implementation:
1. `profile_hmm_align()` - Main HMM alignment
2. `profile_hmm_fill_generic_r9()` - HMM forward pass
3. `make_transition_probs_r9()` - Transition probability calculation
4. Full traceback through HMM states
5. `event_alignment_t` output structure

### 8. **CUDA Version Comparison**

**F5C CUDA (`src/align.cu`):**
- Implements **Adaptive Banded Event Alignment** on GPU
- NOT the full HMM eventalign
- Uses banded DP with kernels:
  - `align_kernel_pre_2d` - Pre-compute kmer ranks
  - `align_kernel_core_2d_shm` - Main DP filling with shared memory
  - `align_kernel_post` - Traceback

**PyFin CUDA (`fin/_f5c/align.cu`):**
- Similar structure to f5c's ABEA
- But you want eventalign, not ABEA!

## Recommendations

### Option 1: Keep Current Implementation (ABEA)
- Rename to `abea()` or `adaptive_align()`
- Document that it's NOT eventalign
- Use for fast rough alignment

### Option 2: Implement Full F5C Eventalign
Port these files from f5c:
1. `src/eventalign.c` → Profile HMM alignment
2. `src/hmm.c` → HMM functions  
3. Update to use full `event_alignment_t` output
4. Implement proper traceback

### Option 3: Hybrid Approach
1. Use current `align.c` for initial ABEA
2. Add Profile HMM refinement on top
3. Output detailed event alignment

## Testing Strategy

To verify you match f5c:

```bash
# Run f5c eventalign
f5c eventalign -b reads.bam -g ref.fa -r reads.fastq > f5c_output.tsv

# Run your eventalign
python -c "
from fin import eventalign
result = eventalign(signal, sequence)
# Compare outputs
"

# Compare:
# - Number of aligned events
# - Event-to-kmer mappings
# - HMM states ('M', 'K', 'B')
# - Reference positions
```

## Summary

❌ **Current PyFin ≠ F5C Eventalign**
- PyFin implements ABEA (adaptive banded alignment)
- F5C eventalign uses Profile HMM
- Different algorithms, different outputs
- Cannot directly compare results

✅ **To match F5C eventalign:**
- Need to port full HMM implementation
- Use `event_alignment_t` structure
- Implement `profile_hmm_align()`
- Add proper state traceback ('M', 'K', 'B')
