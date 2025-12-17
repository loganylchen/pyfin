# Analysis: Pyfin Eventalign vs F5C Differences

## Key Issue Identified

Your pyfin eventalign results differ from f5c because of **event ordering inconsistency**.

## The Problem

### Event Detection (get_events in Python)
- **Uses**: `getevents_raw_order()` in C
- **Returns**: Events in **RAW signal order** (no reversal)
- Event[0] = first event in raw signal = 3' end of RNA
- Event[n-1] = last event in raw signal = 5' end of RNA

### Eventalign Internal (profile_hmm_eventalign in C)
- **Uses**: `getevents_simple()` internally
- **Returns**: Events **REVERSED** to match 5'→3' sequence
- Event[0] = 5' end of RNA (matches sequence[0])
- Event[n-1] = 3' end of RNA

### The Mismatch
When eventalign converts event indices back to raw signal positions:
```c
// In eventalign.c line 263:
int raw_event_idx = (int)et.n - 1 - internal_event_idx;
```

This assumes internal events were **reversed**, but if external code uses `get_events()` directly, those events are **not reversed**, causing index mismatches.

## How F5C Works

F5C eventalign:
1. Detects events in raw order
2. **Reverses** them internally for alignment (5'→3' matching)
3. Outputs with start_idx/end_idx in **raw signal coordinates**

The key is that f5c's eventalign output (TSV file) has:
- `start_idx` and `end_idx`: Raw signal positions (unreversed)
- `event_index`: Index in the **reversed** event array (5'→3' order)

## Differences You're Seeing

### 1. Event Index Order
- **Pyfin**: May be using raw order event indices
- **F5C**: Uses reversed event indices in output

### 2. Signal Position Mapping
- **Pyfin**: `signal_start` calculation may not match f5c
- **F5C**: Uses actual raw signal start/end positions

### 3. Event-to-Kmer Alignment
- Both should align the same way IF events are properly reversed
- Differences suggest event reversal inconsistency

## Root Cause

The issue is in how event indices are interpreted:

**In pyfin eventalign (eventalign.c)**:
```c
// Line 154: Uses getevents_simple() which REVERSES events
event_table et = getevents_simple(nsample, rawptr);

// Line 263: Converts back assuming events were reversed
int raw_event_idx = (int)et.n - 1 - internal_event_idx;
```

**In pyfin visualization/output**:
- If you use `get_events()` separately (Python wrapper), those are NOT reversed
- This creates a mismatch when trying to correlate event indices

## Solution

You need to ensure consistency:

### Option 1: Always Use Reversed Events Internally
- Make sure all internal alignment uses reversed events
- Convert to raw signal coordinates only for output
- Document which coordinates are "event space" vs "signal space"

### Option 2: Standardize on Raw Order
- Remove event reversal entirely
- Reverse the **sequence** instead before alignment
- Keep all indices in raw signal order

### Option 3: Explicit Coordinate Systems
- Add flag to indicate "reversed" vs "raw" event order
- Always specify which coordinate system you're using
- Add conversion functions between systems

## Recommended Fix

The cleanest solution is to match f5c exactly:

1. **Internal alignment**: Use reversed events (5'→3' order)
2. **Output coordinates**: Always in raw signal space
3. **Event indices**: In reversed order (matching sequence positions)

This means:
- `event_index` in output = index in reversed event array
- `signal_start`, `signal_length` = raw signal coordinates (unreversed)
- Users provide 5'→3' sequences (no manual reversal needed)

Would you like me to implement this fix to make your eventalign match f5c's behavior exactly?
