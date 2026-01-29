# PyFIN Eventalign Implementation Plan
## Achieving Result Parity with f5c

This document outlines the detailed plan to make PyFIN's eventalign produce identical results to f5c.

---

## Executive Summary

PyFIN's current eventalign implementation only performs Stage 1 (ABEA - Adaptive Banded Event Alignment). To match f5c's output, we need to implement the complete 2-stage alignment pipeline including Profile HMM realignment and proper alignment QC handling.

---

## Current Status

### What We Have
- ✅ Event detection (from Scrappie, RNA parameters)
- ✅ ABEA alignment (Stage 1)
- ✅ Model of the Mean (MoM) scaling estimation
- ✅ RNA002/RNA004 pore model loading
- ✅ Basic alignment export

### What's Missing
- ❌ Profile HMM realignment (Stage 2)
- ❌ Iterative recalibration loop
- ❌ Complete HMM states (K/B/S missing)
- ❌ `realign_read()` function
- ❌ Proper alignment QC with recovery

---

## Testing Strategy (Phase 1 - PRIORITY)

### Test Framework Goals
1. **Quantitative Comparison**: Exact numeric comparison of event alignments
2. **Position-Level Accuracy**: Compare event-to-kmer mappings position by position
3. **Signal Coordinate Validation**: Verify start_idx/end_idx matches
4. **Scaling Parameter Comparison**: Compare shift/scale/var parameters
5. **Regression Testing**: Ensure changes don't break existing functionality

### Test Metrics
| Metric | Description | Target |
|--------|-------------|--------|
| Position Match Rate | % of positions with same event assignments | >95% |
| Event Index Correlation | Correlation of event indices | >0.99 |
| Mean Value RMSE | Root mean square error of event means | <0.1 pA |
| Scaling Parameter Match | Difference in scale/shift | <0.01 |
| Coverage Match | % of reference covered by events | >98% |

### Test Data Requirements

**RNA004 Test Data:**
- `RNA004.test.pod5`: Raw signal file (10 reads)
- `RNA004.test.fastq.gz`: Basecalled sequences
- `RNA004.test.bam`: Alignments to reference
- `RNA004.test.eventalign.tsv.gz`: f5c reference output (to be generated)

**RNA002 Test Data:**
- `RNA002.test.pod5`: Raw signal file (10 reads)
- `RNA002.test.fastq.gz`: Basecalled sequences
- `RNA002.test.bam`: Alignments to reference
- `RNA002.test.eventalign.tsv.gz`: f5c reference output (to be generated)

**Shared:**
- `reference.fa`: Reference sequence that reads map to

### Test Levels
1. **Unit Tests**: Individual function testing
   - Event detection comparison
   - Scaling parameter estimation
   - K-mer ranking
   
2. **Integration Tests**: Pipeline comparison
   - Full eventalign pipeline vs f5c
   - Event-by-event comparison
   
3. **Regression Tests**: Prevent breaking changes
   - Save baseline results
   - Compare against baseline after changes

---

## Implementation Phases

### Phase 1: Testing Framework (Week 1)
**Goal**: Establish comprehensive testing before any code changes

**Tasks**:
1. Create `test_eventalign_vs_f5c.py` with detailed comparison metrics
2. Create baseline result snapshots
3. Document current discrepancies quantitatively
4. Set up CI-compatible test runner

**Deliverables**:
- `tests/test_eventalign_vs_f5c.py` - Main comparison test
- `tests/test_event_detection.py` - Event detection unit tests
- `tests/test_scaling.py` - Scaling estimation tests
- `tests/fixtures/` - Test data and expected results
- `tests/generate_f5c_reference.sh` - Script to generate f5c reference outputs

---

### Phase 2: Profile HMM Implementation (Weeks 2-3)
**Goal**: Implement Stage 2 alignment (Profile HMM)

**Background**:
f5c uses a 2-stage alignment:
1. **Stage 1 (ABEA)**: Fast banded alignment to get initial event-to-kmer mapping
2. **Stage 2 (Profile HMM)**: Refined alignment using full HMM with all states

**Tasks**:
1. Add missing HMM states to `AlignmentState` enum:
   ```cpp
   enum AlignmentState {
       MATCH = 0,      // Event matches model kmer (existing)
       EVENT_SPLIT = 1, // Event split (existing) 
       KMER = 2,       // New: Skip kmer (no event)
       BAD = 3,        // New: Bad event (outlier)
       SHORTCUT = 4    // New: Large skip
   };
   ```

2. Implement `profile_hmm_align()` function:
   ```cpp
   event_alignment_t profile_hmm_align(
       const model_t* cpg_model,
       const event_table* events,
       const uint32_t* kmer_ranks,
       uint32_t n_kmers,
       const scalings_t& scalings,
       uint32_t start_event,
       uint32_t end_event
   );
   ```

3. Implement transition probabilities:
   - Match → Match: p_stay
   - Match → Event_Split: p_skip_event
   - Match → Kmer: p_skip_kmer
   - etc.

**Key f5c References**:
- `hmm.cpp`: Profile HMM forward/backward algorithms
- `align.cpp:realign_read()`: Combines Stage 1 + Stage 2

---

### Phase 3: Realign Read Function (Week 3)
**Goal**: Implement the complete `realign_read()` pipeline

**Current Flow** (PyFIN):
```
events → ABEA align → output
```

**Target Flow** (f5c):
```
events → ABEA align → Profile HMM realign → merge results → output
```

**Tasks**:
1. Implement `realign_read()` wrapper:
   ```cpp
   void realign_read(
       db_t* db,
       const model_t* model,
       const index_pair_t* ref_range
   ) {
       // Stage 1: ABEA alignment
       event_alignment_t initial = align(...);
       
       // Stage 2: Profile HMM realignment on sub-regions
       for each kmer region:
           event_alignment_t refined = profile_hmm_align(...);
           merge(results, refined);
       
       // Post-processing
       postalign(...);
   }
   ```

2. Handle alignment failure recovery:
   - If ABEA fails QC, attempt recovery with Profile HMM
   - Log statistics for debugging

---

### Phase 4: Iterative Recalibration (Week 4)
**Goal**: Match f5c's scaling parameter estimation exactly

**Current Implementation**:
```cpp
scalings_t estimate_scalings_using_mom(events, model, kmer_ranks);
// Single pass, no iteration
```

**Target Implementation** (f5c):
```cpp
for (int iter = 0; iter < 3; iter++) {
    scalings = estimate_scalings_using_mom(events, model, alignment);
    alignment = align_with_scalings(events, model, scalings);
}
// Final scaling with trimmed alignment
```

**Tasks**:
1. Add iterative recalibration loop
2. Use alignment output to refine scalings
3. Match f5c's outlier trimming parameters

---

### Phase 5: Alignment QC and Recovery (Week 4-5)
**Goal**: Match f5c's alignment quality control behavior

**Current Behavior** (PyFIN):
```cpp
if (events_per_base < 0.5 || events_per_base > 20) {
    // Mark as failed
    mapping.status = 4; // Alignment QC failed
    return; // No output
}
```

**Target Behavior** (f5c):
```cpp
if (failed_qc) {
    // Try recovery with expanded bands
    realign_read_with_recovery(...);
    if (still_failed) {
        // Still output partial results with warning
        output_partial_alignment(...);
    }
}
```

**Tasks**:
1. Implement alignment QC recovery
2. Add partial output mode
3. Match f5c's QC thresholds exactly

---

## File-by-File Changes

### `fin/_eventalign/align.cpp`

| Function | Change Type | Description |
|----------|-------------|-------------|
| `align()` | Modify | Add Profile HMM fallback |
| `postalign()` | Keep | No changes needed |
| `realign_read()` | **NEW** | Complete 2-stage pipeline |
| `profile_hmm_align()` | **NEW** | Stage 2 HMM alignment |
| `merge_alignments()` | **NEW** | Combine ABEA + HMM results |

### `fin/_eventalign/common.h`

| Change | Description |
|--------|-------------|
| Add `AlignmentState::KMER` | Skip kmer state |
| Add `AlignmentState::BAD` | Bad event state |
| Add `AlignmentState::SHORTCUT` | Large skip state |
| Add HMM transition params | `p_skip_kmer`, etc. |

### `fin/_eventalign/events.cpp`

| Function | Change Type | Description |
|----------|-------------|-------------|
| `getevents()` | Keep | Already correct for RNA |

### `fin/_eventalign/model.cpp`

| Function | Change Type | Description |
|----------|-------------|-------------|
| `set_model_from_file()` | Keep | No changes |
| `set_model()` | Keep | No changes |

### `fin/_eventalign/event_api_wrapper.cpp`

| Function | Change Type | Description |
|----------|-------------|-------------|
| `run_eventalign()` | Modify | Call `realign_read()` |
| `profile_hmm_eventalign()` | **NEW** | Python binding for HMM |

---

## Priority Matrix

| Task | Priority | Complexity | Dependencies |
|------|----------|------------|--------------|
| Testing Framework | **P0** | Low | None |
| Profile HMM States | P1 | Medium | Testing |
| Profile HMM Align | P1 | High | States |
| realign_read() | P1 | High | Profile HMM |
| Iterative Recalib | P2 | Medium | realign_read |
| Alignment QC | P2 | Medium | realign_read |

---

## Success Criteria

### Minimum Viable Parity
- Position match rate >95%
- Event mean RMSE <0.5 pA
- Scaling parameters within 1%

### Full Parity
- Position match rate >99%
- Event mean RMSE <0.1 pA
- Scaling parameters within 0.1%
- Identical output format

---

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| f5c uses different constants | Medium | Extract exact values from f5c source |
| Numerical precision differences | Low | Use f5c's exact formulas |
| Edge cases in HMM | High | Comprehensive test coverage |
| Performance regression | Medium | Profile and optimize |

---

## Timeline

```
Week 1: Testing Framework (Phase 1)
  ├── Day 1-2: Create test_eventalign_vs_f5c.py
  ├── Day 3-4: Establish baseline metrics
  └── Day 5: Document current discrepancies

Week 2-3: Profile HMM (Phase 2)
  ├── Day 1-3: Add HMM states and structures
  ├── Day 4-7: Implement profile_hmm_align()
  └── Day 8-10: Unit tests and debugging

Week 3: Realign Read (Phase 3)
  ├── Day 1-3: Implement realign_read()
  └── Day 4-5: Integration testing

Week 4-5: Refinements (Phases 4-5)
  ├── Day 1-3: Iterative recalibration
  ├── Day 4-5: Alignment QC
  └── Day 6-7: Final testing and documentation
```

---

## Appendix: f5c Code References

### Key Files in f5c
- `src/align.cpp`: Main alignment logic
- `src/hmm.cpp`: Profile HMM implementation
- `src/event.c`: Event detection (from Scrappie)
- `src/model.c`: Pore model handling
- `src/meth.c`: Eventalign main logic

### Key Functions
- `realign_read()`: Main entry point
- `profile_hmm_align()`: Profile HMM aligner
- `align()`: ABEA implementation
- `adaptive_banded_simple_event_align()`: Core alignment

### Key Constants
```cpp
// From f5c
#define EVENT_DETECTION_THRESHOLD_MEAN 1.4
#define MIN_AVERAGE_EVENTS_PER_KMER 0.5
#define MAX_AVERAGE_EVENTS_PER_KMER 20.0
#define READS_FAILED_QC 4
```
