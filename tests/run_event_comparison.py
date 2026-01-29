#!/usr/bin/env python3
"""
Comprehensive Event Detection Comparison: PyFIN vs f5c

Tests event detection on RNA004 data and produces a summary report.
"""
import sys
sys.path.insert(0, "/home/logan/Projects/pyfin")

from pathlib import Path
import gzip
from collections import defaultdict
import numpy as np

print("=" * 70)
print("EVENT DETECTION COMPARISON: PyFIN getevents vs f5c")
print("=" * 70)

TEST_DATA_DIR = Path("/home/logan/Projects/pyfin/tests/testdata")
RNA004_F5C_TSV_PATH = TEST_DATA_DIR / "RNA004.test.tsv.gz"
RNA004_POD5_PATH = TEST_DATA_DIR / "RNA004.test.pod5"

# Load f5c events (limit to 50 reads for testing)
MAX_READS = 50
print(f"\n[1] Loading f5c eventalign results (max {MAX_READS} reads)...")
events_by_read = defaultdict(list)
total_f5c_events = 0

with gzip.open(RNA004_F5C_TSV_PATH, "rt") as f:
    header = f.readline()  # Skip header
    
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 15:
            read_id = parts[3]
            
            # Limit to MAX_READS
            if read_id not in events_by_read and len(events_by_read) >= MAX_READS:
                continue
            
            events_by_read[read_id].append({
                "start": int(parts[13]),
                "end": int(parts[14]),
                "mean": float(parts[6]),
                "stdv": float(parts[7]),
                "event_idx": int(parts[5]),
            })
            total_f5c_events += 1

# Sort events by start position within each read
for read_id in events_by_read:
    events_by_read[read_id].sort(key=lambda e: e["start"])

print(f"  Loaded {len(events_by_read)} reads with {total_f5c_events:,} total events")

# Load POD5 signals
print("\n[2] Loading POD5 signals...")
import pod5

read_signals = {}
with pod5.Reader(str(RNA004_POD5_PATH)) as reader:
    for read in reader.reads():
        rid = str(read.read_id)
        if rid in events_by_read:
            read_signals[rid] = read.signal_pa.astype(np.float32)

print(f"  Loaded {len(read_signals)} matching signals")

# Run comparison
print("\n[3] Running PyFIN event detection...")
from fin._eventalign import getevents

# Collect metrics for each read
all_metrics = []

for read_id in events_by_read:
    if read_id not in read_signals:
        continue
    
    signal = read_signals[read_id]
    f5c_events = events_by_read[read_id]
    
    # Run PyFIN getevents
    result = getevents(signal)
    
    # Extract PyFIN events
    pyfin_n = result['n_events']
    pyfin_starts = result['starts']
    pyfin_ends = pyfin_starts + result['lengths']
    pyfin_means = result['means']
    pyfin_stdvs = result['stdvs']
    
    # Extract f5c events
    f5c_starts = np.array([e["start"] for e in f5c_events])
    f5c_ends = np.array([e["end"] for e in f5c_events])
    f5c_means = np.array([e["mean"] for e in f5c_events])
    f5c_stdvs = np.array([e["stdv"] for e in f5c_events])
    
    # Match events
    matched = 0
    start_diffs = []
    end_diffs = []
    mean_diffs = []
    stdv_diffs = []
    f5c_matched_means = []
    pyfin_matched_means = []
    
    for i in range(len(f5c_events)):
        f_start = f5c_starts[i]
        f_end = f5c_ends[i]
        f_mean = f5c_means[i]
        f_stdv = f5c_stdvs[i]
        
        # Find overlapping PyFIN event
        overlaps = (pyfin_starts < f_end) & (pyfin_ends > f_start)
        if np.any(overlaps):
            idx = np.where(overlaps)[0][0]
            start_diffs.append(int(pyfin_starts[idx]) - f_start)
            end_diffs.append(int(pyfin_ends[idx]) - f_end)
            mean_diffs.append(pyfin_means[idx] - f_mean)
            stdv_diffs.append(pyfin_stdvs[idx] - f_stdv)
            f5c_matched_means.append(f_mean)
            pyfin_matched_means.append(pyfin_means[idx])
            matched += 1
    
    # Calculate metrics
    metrics = {
        "read_id": read_id,
        "signal_len": len(signal),
        "f5c_events": len(f5c_events),
        "pyfin_events": pyfin_n,
        "matched": matched,
        "match_rate": matched / len(f5c_events) if f5c_events else 0,
        "ratio": pyfin_n / len(f5c_events) if f5c_events else 0,
    }
    
    if start_diffs:
        metrics["mean_start_diff"] = np.mean(start_diffs)
        metrics["std_start_diff"] = np.std(start_diffs)
        metrics["mean_end_diff"] = np.mean(end_diffs)
        metrics["std_end_diff"] = np.std(end_diffs)
        metrics["mean_mean_diff"] = np.mean(mean_diffs)
        metrics["std_mean_diff"] = np.std(mean_diffs)
        metrics["mean_stdv_diff"] = np.mean(stdv_diffs)
        metrics["std_stdv_diff"] = np.std(stdv_diffs)
        
        if len(f5c_matched_means) > 1:
            metrics["mean_correlation"] = np.corrcoef(f5c_matched_means, pyfin_matched_means)[0, 1]
        else:
            metrics["mean_correlation"] = np.nan
    
    all_metrics.append(metrics)

# Print per-read results
print("\n" + "=" * 70)
print("PER-READ RESULTS")
print("=" * 70)

for m in all_metrics:
    print(f"\n  Read: {m['read_id'][:12]}...")
    print(f"    Signal: {m['signal_len']:,} samples")
    print(f"    Events: f5c={m['f5c_events']:,}, PyFIN={m['pyfin_events']:,} (ratio={m['ratio']:.3f})")
    print(f"    Match rate: {m['match_rate']*100:.1f}% ({m['matched']:,}/{m['f5c_events']:,})")
    if "mean_start_diff" in m:
        print(f"    Boundary diff: start={m['mean_start_diff']:.1f}±{m['std_start_diff']:.1f}, end={m['mean_end_diff']:.1f}±{m['std_end_diff']:.1f}")
        print(f"    Mean diff: {m['mean_mean_diff']:.2f}±{m['std_mean_diff']:.2f} pA")
        if not np.isnan(m["mean_correlation"]):
            print(f"    Mean correlation: {m['mean_correlation']:.6f}")

# Print summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

total_f5c = sum(m["f5c_events"] for m in all_metrics)
total_pyfin = sum(m["pyfin_events"] for m in all_metrics)
total_matched = sum(m["matched"] for m in all_metrics)

print(f"\n  Total reads tested: {len(all_metrics)}")
print(f"  Total f5c events: {total_f5c:,}")
print(f"  Total PyFIN events: {total_pyfin:,}")
print(f"  Overall event ratio: {total_pyfin/total_f5c:.3f}")
print(f"  Overall match rate: {total_matched/total_f5c*100:.1f}%")

# Aggregate metrics
avg_start_diff = np.mean([m.get("mean_start_diff", np.nan) for m in all_metrics if "mean_start_diff" in m])
avg_end_diff = np.mean([m.get("mean_end_diff", np.nan) for m in all_metrics if "mean_end_diff" in m])
avg_mean_diff = np.mean([m.get("mean_mean_diff", np.nan) for m in all_metrics if "mean_mean_diff" in m])
avg_correlation = np.nanmean([m.get("mean_correlation", np.nan) for m in all_metrics if "mean_correlation" in m])

print(f"\n  Average boundary alignment:")
print(f"    Start: {avg_start_diff:.2f} samples")
print(f"    End: {avg_end_diff:.2f} samples")
print(f"\n  Average event mean difference: {avg_mean_diff:.2f} pA")
print(f"  Average mean correlation: {avg_correlation:.6f}")

# Interpretation
print("\n" + "-" * 70)
print("INTERPRETATION")
print("-" * 70)

if total_matched/total_f5c >= 0.99:
    print("\n  ✓ EXCELLENT: Event boundaries align perfectly (≥99% match)")
elif total_matched/total_f5c >= 0.90:
    print("\n  ◐ GOOD: Most event boundaries align (≥90% match)")
else:
    print("\n  ✗ NEEDS WORK: Event boundary alignment <90%")

if abs(avg_start_diff) <= 1 and abs(avg_end_diff) <= 1:
    print("  ✓ EXCELLENT: Boundary positions match within ±1 sample")
elif abs(avg_start_diff) <= 10 and abs(avg_end_diff) <= 10:
    print("  ◐ GOOD: Boundary positions match within ±10 samples")
else:
    print("  ✗ NEEDS WORK: Boundary positions differ by >10 samples")

if avg_correlation >= 0.999:
    print("  ✓ EXCELLENT: Event means are perfectly correlated (≥0.999)")
elif avg_correlation >= 0.99:
    print("  ◐ GOOD: Event means are highly correlated (≥0.99)")
else:
    print("  ✗ NEEDS WORK: Event mean correlation <0.99")

if abs(avg_mean_diff) > 1:
    print(f"\n  NOTE: Systematic mean offset of {avg_mean_diff:.2f} pA detected")
    print("        This may indicate different signal scaling/normalization")

print("\n" + "=" * 70)
print("TEST COMPLETE")
print("=" * 70)
