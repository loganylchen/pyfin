# Implementation of Full f5c 3-State HMM Alignment

## Overview

This document describes the implementation of f5c's full 3-state HMM alignment algorithm in `align.c`, replacing the simplified 2-state Viterbi algorithm.

## Problem

The initial implementation used a simplified alignment:
- Only 2 states per kmer (stay vs move)
- Simple emission model: `(event_mean - shift) / scale`
- No handling of bad events or kmer skips
- Fixed transition probabilities

F5C uses a sophisticated 3-state HMM:
- **3 states**: MATCH, BAD_EVENT, KMER_SKIP
- **10 transition types**: Match-to-match (self/next), bad-to-bad, bad-to-match (self/next), match-to-bad, match-to-skip, bad-to-skip, skip-to-skip, skip-to-match
- **Dynamic transition probabilities** based on events_per_base ratio
- **Sophisticated emission model**: `gp_mean = scale * model_mean + shift`

## Changes Made

### 1. Added HMM Structures (align.c)

```c
// Movement types for traceback
typedef enum {
    HMT_FROM_SAME_M = 0,    // From MATCH, same kmer
    HMT_FROM_PREV_M,        // From MATCH, previous kmer
    HMT_FROM_SAME_B,        // From BAD_EVENT, same kmer
    HMT_FROM_PREV_B,        // From BAD_EVENT, previous kmer
    HMT_FROM_PREV_K,        // From KMER_SKIP, previous kmer
    HMT_FROM_SOFT,          // From soft clip (flanking)
    HMT_NUM_MOVEMENT_TYPES
} HMMMovementType;

// Scores for each possible transition
typedef struct {
    float x[HMT_NUM_MOVEMENT_TYPES];
} HMMUpdateScores;

// Transition probabilities for each kmer block
typedef struct {
    float lp_mm_self;   // MATCH -> MATCH (same kmer)
    float lp_mm_next;   // MATCH -> MATCH (next kmer)
    float lp_mb;        // MATCH -> BAD_EVENT
    float lp_mk;        // MATCH -> KMER_SKIP
    float lp_bb;        // BAD_EVENT -> BAD_EVENT
    float lp_bk;        // BAD_EVENT -> KMER_SKIP
    float lp_bm_next;   // BAD_EVENT -> MATCH (next kmer)
    float lp_bm_self;   // BAD_EVENT -> MATCH (same kmer)
    float lp_kk;        // KMER_SKIP -> KMER_SKIP
    float lp_km;        // KMER_SKIP -> MATCH
} BlockTransitions;
```

### 2. Implemented F5C Emission Model

```c
static inline float log_probability_match(
    simple_scalings_t scaling,
    simple_model_t *model,
    event_t *events,
    int32_t event_idx,
    uint32_t rank)
{
    // F5C scaling model
    float gp_mean = scaling.scale * model[rank].level_mean + scaling.shift;
    float gp_stdv = model[rank].level_stdv * scaling.var;
    float gp_log_stdv = logf(gp_stdv);
    
    return log_normal_pdf(events[event_idx].mean, gp_mean, gp_stdv, gp_log_stdv);
}
```

### 3. Dynamic Transition Probabilities

```c
static void calculate_transitions(
    BlockTransitions *transitions,
    int n_kmers,
    float events_per_base)
{
    // F5C transition model
    float p_stay = 1.0f - (1.0f / events_per_base);
    float p_skip = 0.0025f;
    float p_bad = 0.001f;
    
    for (int i = 0; i < n_kmers; ++i)
    {
        transitions[i].lp_mm_self = logf(p_stay);
        transitions[i].lp_mm_next = logf(1.0f - p_stay - p_skip - p_bad);
        transitions[i].lp_mb = logf(p_bad);
        transitions[i].lp_mk = logf(p_skip);
        transitions[i].lp_bb = logf(p_stay);
        transitions[i].lp_bk = logf(p_skip);
        transitions[i].lp_bm_next = logf(1.0f - p_stay - p_skip);
        transitions[i].lp_bm_self = logf(0.01f);
        transitions[i].lp_kk = logf(p_stay);
        transitions[i].lp_km = logf(1.0f - p_stay);
    }
}
```

### 4. Three-State DP Table

Old: `dp[n_events][n_kmers]` - 2D table

New: `dp[n_events + 1][n_kmers * 3]` - 3 states per kmer
- State 0: MATCH (event matches kmer)
- State 1: BAD_EVENT (event should be ignored)
- State 2: KMER_SKIP (kmer with no events)

### 5. HMM Update Logic

For each event and each kmer, update all 3 states:

**STATE_MATCH**: Event aligns to this kmer
- Can come from: same MATCH, prev MATCH, same BAD, prev BAD, prev SKIP, soft clip
- Emission: log probability of event given kmer model

**STATE_BAD_EVENT**: Event is noise, ignore it
- Can come from: same MATCH, same BAD
- Emission: 0 (penalty already in transition)

**STATE_KMER_SKIP**: Kmer has no event
- Can come from: prev MATCH, prev BAD, prev SKIP
- Emission: 0 (no event consumed)
- Special: Doesn't consume event (same row)

### 6. Traceback Through 3 States

Track best previous state and kmer for each cell:
- `traceback_state[row][offset + state]` - which state we came from
- `traceback_kmer[row][offset + state]` - which kmer block we came from

Walk backwards from best ending state, recording:
- MATCH states → output as aligned pairs
- BAD_EVENT states → recorded but not output (skipped events)
- KMER_SKIP states → don't consume events (no decrement of row)

### 7. Updated Function Signatures

Changed from pre-allocated output buffer to dynamic allocation:

```c
// Old
int32_t align_with_flanking_cpu(
    simple_aligned_pair_t *out,  // Pre-allocated buffer
    ...);

// New
int32_t align_with_flanking_cpu(
    simple_aligned_pair_t **out_alignment,  // Pointer-to-pointer, function allocates
    ...);
```

## Key Differences from Simplified Algorithm

| Aspect | Simplified | F5C 3-State HMM |
|--------|-----------|-----------------|
| States | 1 per kmer | 3 per kmer (MATCH, BAD, SKIP) |
| Transitions | 2 types (stay/move) | 10 types |
| Bad events | Not handled | Explicitly modeled |
| Kmer skips | Not handled | Explicitly modeled |
| Emission | Simple scaling | Full variance model |
| Transitions | Fixed | Dynamic (events_per_base) |
| Memory | O(events × kmers) | O(events × kmers × 3) |

## Benefits

1. **Better alignment quality**: Properly handles noisy events
2. **Robust to bad signals**: BAD_EVENT state filters noise
3. **Handles missing events**: KMER_SKIP for unmeasured kmers
4. **Adaptive**: Transitions adjust to event density
5. **Matches f5c**: Same algorithm as original tool

## Testing Status

- ✓ Code implemented in `align.c`
- ✓ Updated `align_common.h` function signature
- ✓ Updated `eventalign.c` to use new signature
- ⚠ `align.cu` (GPU) needs matching update
- ⚠ Compilation not yet tested
- ⚠ Alignment quality not yet verified

## Next Steps

1. Update `align.cu` to match 3-state HMM
2. Test compilation of CPU version
3. Run on real nanopore data
4. Compare alignment with original f5c
5. Benchmark performance vs simplified version
