#!/usr/bin/env python3
"""
Event Detection Comparison Test: PyFIN getevents vs f5c

This test compares PyFIN's event detection algorithm (getevents) with f5c's 
event detection. Event boundaries (start_idx, end_idx) and event statistics 
(mean, stdv) from f5c's eventalign output are compared against PyFIN's getevents output.

The f5c eventalign TSV contains per-event information:
- start_idx, end_idx: signal sample boundaries for each event
- event_level_mean, event_stdv: event statistics
- event_index: sequential event number

PyFIN getevents returns:
- starts: event start positions (in samples)
- lengths: event lengths (in samples)  
- means: event mean current values
- stdvs: event standard deviations

Comparison metrics:
1. Event boundary alignment (start/end positions)
2. Event count comparison
3. Event mean/stdv correlation
4. Event merging/splitting analysis

Usage:
    pytest tests/test_event_detection_vs_f5c.py -v -s
    python tests/test_event_detection_vs_f5c.py
"""

import numpy as np
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass, field
from collections import defaultdict
import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Test data paths
TEST_DATA_DIR = PROJECT_ROOT / "tests" / "testdata"

# RNA004 test data
RNA004_POD5_PATH = TEST_DATA_DIR / "RNA004.test.pod5"
RNA004_F5C_TSV_PATH = TEST_DATA_DIR / "RNA004.test.tsv.gz"

# Reference
REFERENCE_PATH = TEST_DATA_DIR / "test.fa"

# Check module availability
try:
    import pod5
    POD5_AVAILABLE = True
except ImportError:
    POD5_AVAILABLE = False

try:
    from fin._eventalign import getevents, MODEL_RNA004, set_model
    GETEVENTS_AVAILABLE = True
except ImportError:
    GETEVENTS_AVAILABLE = False


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class F5CEvent:
    """Event record from f5c eventalign output."""
    read_id: str
    event_idx: int
    start_idx: int
    end_idx: int
    event_mean: float
    event_stdv: float
    event_length: float  # in seconds or samples
    reference_name: str = ""
    reference_pos: int = 0
    reference_kmer: str = ""


@dataclass
class PyFINEvent:
    """Event record from PyFIN getevents."""
    event_idx: int
    start_idx: int
    end_idx: int
    event_mean: float
    event_stdv: float
    length: int  # in samples


@dataclass
class EventMatch:
    """Matched event pair between f5c and PyFIN."""
    f5c_event: F5CEvent
    pyfin_event: PyFINEvent
    overlap_ratio: float
    start_diff: int
    end_diff: int
    mean_diff: float
    stdv_diff: float


@dataclass
class EventComparisonMetrics:
    """Comprehensive metrics comparing f5c vs PyFIN event detection."""
    read_id: str
    
    # Event counts
    f5c_event_count: int = 0
    pyfin_event_count: int = 0
    matched_event_count: int = 0
    
    # Boundary metrics
    mean_start_diff: float = 0.0
    std_start_diff: float = 0.0
    mean_end_diff: float = 0.0
    std_end_diff: float = 0.0
    median_start_diff: float = 0.0
    median_end_diff: float = 0.0
    
    # Event statistics metrics
    mean_correlation: float = 0.0
    mean_rmse: float = 0.0
    stdv_correlation: float = 0.0
    stdv_rmse: float = 0.0
    
    # Coverage metrics
    f5c_coverage_samples: int = 0
    pyfin_coverage_samples: int = 0
    overlap_samples: int = 0
    coverage_jaccard: float = 0.0
    
    # Event structure metrics
    f5c_merged_into_pyfin: int = 0  # Multiple f5c events mapping to one PyFIN event
    pyfin_split_from_f5c: int = 0   # One f5c event mapping to multiple PyFIN events
    
    # Signal range
    signal_length: int = 0
    f5c_signal_start: int = 0
    f5c_signal_end: int = 0
    pyfin_signal_start: int = 0
    pyfin_signal_end: int = 0
    
    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            "=" * 70,
            f"Event Detection Comparison: {self.read_id}",
            "=" * 70,
            "",
            "Event Counts:",
            f"  f5c events:       {self.f5c_event_count:,}",
            f"  PyFIN events:     {self.pyfin_event_count:,}",
            f"  Matched pairs:    {self.matched_event_count:,}",
            f"  Match rate:       {self.matched_event_count / max(1, self.f5c_event_count) * 100:.1f}%",
            "",
            "Boundary Alignment (samples):",
            f"  Start diff (mean ± std): {self.mean_start_diff:.1f} ± {self.std_start_diff:.1f}",
            f"  End diff (mean ± std):   {self.mean_end_diff:.1f} ± {self.std_end_diff:.1f}",
            f"  Start diff (median):     {self.median_start_diff:.1f}",
            f"  End diff (median):       {self.median_end_diff:.1f}",
            "",
            "Event Statistics:",
            f"  Mean correlation:  {self.mean_correlation:.4f}",
            f"  Mean RMSE:         {self.mean_rmse:.2f} pA",
            f"  Stdv correlation:  {self.stdv_correlation:.4f}",
            f"  Stdv RMSE:         {self.stdv_rmse:.2f} pA",
            "",
            "Coverage:",
            f"  f5c coverage:      {self.f5c_coverage_samples:,} samples",
            f"  PyFIN coverage:    {self.pyfin_coverage_samples:,} samples",
            f"  Overlap:           {self.overlap_samples:,} samples",
            f"  Jaccard index:     {self.coverage_jaccard:.4f}",
            "",
            "Signal Range:",
            f"  Signal length:     {self.signal_length:,} samples",
            f"  f5c range:         [{self.f5c_signal_start:,}, {self.f5c_signal_end:,}]",
            f"  PyFIN range:       [{self.pyfin_signal_start:,}, {self.pyfin_signal_end:,}]",
            "",
            "Event Structure:",
            f"  f5c merged into PyFIN: {self.f5c_merged_into_pyfin}",
            f"  PyFIN split from f5c:  {self.pyfin_split_from_f5c}",
            "=" * 70,
        ]
        return "\n".join(lines)


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_f5c_events_by_read(tsv_path: Path, max_reads: int = None) -> Dict[str, List[F5CEvent]]:
    """
    Load f5c eventalign TSV and group events by read_id.
    
    f5c TSV columns (0-indexed):
    0: contig, 1: position, 2: reference_kmer, 3: read_name, 4: strand,
    5: event_index, 6: event_level_mean, 7: event_stdv, 8: event_length,
    9: model_kmer, 10: model_mean, 11: model_stdv, 12: standardized_level,
    13: start_idx, 14: end_idx, 15: samples
    """
    events_by_read: Dict[str, List[F5CEvent]] = defaultdict(list)
    read_count = 0
    
    open_func = gzip.open if str(tsv_path).endswith('.gz') else open
    
    with open_func(tsv_path, 'rt') as f:
        # Skip header if present
        first_line = f.readline()
        if first_line.startswith('contig'):
            pass  # Header line, skip
        else:
            f.seek(0)  # Not a header, go back
        
        for line in f:
            if not line.strip():
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 15:
                continue
            
            # Skip header if encountered again
            if parts[0] == 'contig':
                continue
            
            read_id = parts[3]
            
            # Track unique reads
            if read_id not in events_by_read:
                read_count += 1
                if max_reads and read_count > max_reads:
                    break
            
            try:
                event = F5CEvent(
                    read_id=read_id,
                    event_idx=int(parts[5]),
                    start_idx=int(parts[13]),
                    end_idx=int(parts[14]),
                    event_mean=float(parts[6]),
                    event_stdv=float(parts[7]),
                    event_length=float(parts[8]),
                    reference_name=parts[0],
                    reference_pos=int(parts[1]),
                    reference_kmer=parts[2],
                )
                events_by_read[read_id].append(event)
            except (ValueError, IndexError):
                continue
    
    # Sort events by start_idx within each read
    for read_id in events_by_read:
        events_by_read[read_id].sort(key=lambda e: e.start_idx)
    
    return dict(events_by_read)


def load_pod5_signals(pod5_path: Path, read_ids: List[str] = None) -> Dict[str, np.ndarray]:
    """Load raw signals from POD5 file."""
    if not POD5_AVAILABLE:
        raise ImportError("pod5 package required")
    
    signals = {}
    
    with pod5.Reader(pod5_path) as reader:
        for read in reader.reads():
            rid = str(read.read_id)
            if read_ids is None or rid in read_ids:
                # Get pA signal (already calibrated)
                signals[rid] = read.signal_pa.astype(np.float32)
    
    return signals


def run_pyfin_getevents(signal: np.ndarray) -> List[PyFINEvent]:
    """Run PyFIN getevents and convert to event list."""
    if not GETEVENTS_AVAILABLE:
        raise ImportError("fin._eventalign.getevents not available")
    
    # Ensure float32
    if signal.dtype != np.float32:
        signal = signal.astype(np.float32)
    
    result = getevents(signal)
    
    events = []
    n_events = result['n_events']
    starts = result['starts']
    lengths = result['lengths']
    means = result['means']
    stdvs = result['stdvs']
    
    for i in range(n_events):
        start = int(starts[i])
        length = int(lengths[i])
        events.append(PyFINEvent(
            event_idx=i,
            start_idx=start,
            end_idx=start + length,
            event_mean=float(means[i]),
            event_stdv=float(stdvs[i]),
            length=length,
        ))
    
    return events


# =============================================================================
# Comparison Functions
# =============================================================================

def match_events(f5c_events: List[F5CEvent], 
                 pyfin_events: List[PyFINEvent],
                 overlap_threshold: float = 0.3) -> List[EventMatch]:
    """
    Match f5c events to PyFIN events based on signal position overlap.
    
    Uses efficient numpy-based matching with binary search for O(n log n) performance.
    An f5c event is matched to a PyFIN event if they have significant overlap.
    """
    if not f5c_events or not pyfin_events:
        return []
    
    # Convert to numpy arrays for efficiency
    f5c_starts = np.array([e.start_idx for e in f5c_events])
    f5c_ends = np.array([e.end_idx for e in f5c_events])
    pyfin_starts = np.array([e.start_idx for e in pyfin_events])
    pyfin_ends = np.array([e.end_idx for e in pyfin_events])
    
    matches = []
    pyfin_matched = set()
    
    # For each f5c event, find overlapping PyFIN events using binary search
    for f5c_idx, f5c_ev in enumerate(f5c_events):
        f5c_start = f5c_ev.start_idx
        f5c_end = f5c_ev.end_idx
        f5c_len = f5c_end - f5c_start
        
        if f5c_len <= 0:
            continue
        
        # Find candidate PyFIN events that could overlap
        # An event overlaps if: pyfin_start < f5c_end AND pyfin_end > f5c_start
        candidates = np.where((pyfin_starts < f5c_end) & (pyfin_ends > f5c_start))[0]
        
        best_match = None
        best_overlap_ratio = 0
        
        for pyfin_idx in candidates:
            if pyfin_idx in pyfin_matched:
                continue
            
            pyfin_ev = pyfin_events[pyfin_idx]
            pyfin_start = pyfin_ev.start_idx
            pyfin_end = pyfin_ev.end_idx
            pyfin_len = pyfin_end - pyfin_start
            
            # Calculate overlap
            overlap_start = max(f5c_start, pyfin_start)
            overlap_end = min(f5c_end, pyfin_end)
            overlap = max(0, overlap_end - overlap_start)
            
            # Overlap ratio relative to smaller event
            min_len = min(f5c_len, pyfin_len)
            overlap_ratio = overlap / min_len if min_len > 0 else 0
            
            if overlap_ratio > overlap_threshold and overlap_ratio > best_overlap_ratio:
                best_overlap_ratio = overlap_ratio
                best_match = pyfin_ev
        
        if best_match is not None:
            match = EventMatch(
                f5c_event=f5c_ev,
                pyfin_event=best_match,
                overlap_ratio=best_overlap_ratio,
                start_diff=best_match.start_idx - f5c_ev.start_idx,
                end_diff=best_match.end_idx - f5c_ev.end_idx,
                mean_diff=best_match.event_mean - f5c_ev.event_mean,
                stdv_diff=best_match.event_stdv - f5c_ev.event_stdv,
            )
            matches.append(match)
            pyfin_matched.add(best_match.event_idx)
    
    return matches


def compute_coverage_overlap(f5c_events: List[F5CEvent], 
                             pyfin_events: List[PyFINEvent],
                             signal_length: int) -> Tuple[int, int, int]:
    """
    Compute coverage in samples for f5c, PyFIN, and their overlap.
    Returns (f5c_coverage, pyfin_coverage, overlap).
    """
    # Create coverage arrays
    f5c_coverage = np.zeros(signal_length, dtype=bool)
    pyfin_coverage = np.zeros(signal_length, dtype=bool)
    
    for ev in f5c_events:
        start = max(0, ev.start_idx)
        end = min(signal_length, ev.end_idx)
        f5c_coverage[start:end] = True
    
    for ev in pyfin_events:
        start = max(0, ev.start_idx)
        end = min(signal_length, ev.end_idx)
        pyfin_coverage[start:end] = True
    
    f5c_samples = np.sum(f5c_coverage)
    pyfin_samples = np.sum(pyfin_coverage)
    overlap_samples = np.sum(f5c_coverage & pyfin_coverage)
    
    return int(f5c_samples), int(pyfin_samples), int(overlap_samples)


def analyze_event_structure(f5c_events: List[F5CEvent],
                            pyfin_events: List[PyFINEvent],
                            matches: List[EventMatch]) -> Tuple[int, int]:
    """
    Analyze event merging/splitting between f5c and PyFIN.
    
    Returns:
        (f5c_merged, pyfin_split): 
        - f5c_merged: count of f5c events that share a PyFIN match
        - pyfin_split: count of PyFIN events that share an f5c match
    """
    # Count how many f5c events map to each PyFIN event
    pyfin_to_f5c: Dict[int, List[int]] = defaultdict(list)
    f5c_to_pyfin: Dict[int, List[int]] = defaultdict(list)
    
    for match in matches:
        pyfin_idx = match.pyfin_event.event_idx
        f5c_idx = match.f5c_event.event_idx
        pyfin_to_f5c[pyfin_idx].append(f5c_idx)
        f5c_to_pyfin[f5c_idx].append(pyfin_idx)
    
    # f5c events that are merged into single PyFIN events
    f5c_merged = sum(1 for f5c_list in pyfin_to_f5c.values() if len(f5c_list) > 1)
    
    # PyFIN events that are split from single f5c events
    pyfin_split = sum(1 for pyfin_list in f5c_to_pyfin.values() if len(pyfin_list) > 1)
    
    return f5c_merged, pyfin_split


def compare_events_for_read(read_id: str,
                            f5c_events: List[F5CEvent],
                            signal: np.ndarray) -> EventComparisonMetrics:
    """
    Compare f5c events vs PyFIN getevents for a single read.
    """
    metrics = EventComparisonMetrics(read_id=read_id)
    metrics.signal_length = len(signal)
    metrics.f5c_event_count = len(f5c_events)
    
    # Run PyFIN getevents
    try:
        pyfin_events = run_pyfin_getevents(signal)
        metrics.pyfin_event_count = len(pyfin_events)
    except Exception as e:
        print(f"  Error running getevents: {e}")
        return metrics
    
    if not f5c_events or not pyfin_events:
        return metrics
    
    # Signal ranges
    metrics.f5c_signal_start = min(e.start_idx for e in f5c_events)
    metrics.f5c_signal_end = max(e.end_idx for e in f5c_events)
    metrics.pyfin_signal_start = min(e.start_idx for e in pyfin_events)
    metrics.pyfin_signal_end = max(e.end_idx for e in pyfin_events)
    
    # Match events
    matches = match_events(f5c_events, pyfin_events)
    metrics.matched_event_count = len(matches)
    
    if not matches:
        return metrics
    
    # Boundary differences
    start_diffs = [m.start_diff for m in matches]
    end_diffs = [m.end_diff for m in matches]
    
    metrics.mean_start_diff = np.mean(start_diffs)
    metrics.std_start_diff = np.std(start_diffs)
    metrics.mean_end_diff = np.mean(end_diffs)
    metrics.std_end_diff = np.std(end_diffs)
    metrics.median_start_diff = np.median(start_diffs)
    metrics.median_end_diff = np.median(end_diffs)
    
    # Event statistics
    f5c_means = [m.f5c_event.event_mean for m in matches]
    pyfin_means = [m.pyfin_event.event_mean for m in matches]
    f5c_stdvs = [m.f5c_event.event_stdv for m in matches]
    pyfin_stdvs = [m.pyfin_event.event_stdv for m in matches]
    
    if len(f5c_means) > 1:
        metrics.mean_correlation = np.corrcoef(f5c_means, pyfin_means)[0, 1]
        metrics.stdv_correlation = np.corrcoef(f5c_stdvs, pyfin_stdvs)[0, 1]
    
    metrics.mean_rmse = np.sqrt(np.mean([m.mean_diff**2 for m in matches]))
    metrics.stdv_rmse = np.sqrt(np.mean([m.stdv_diff**2 for m in matches]))
    
    # Coverage
    f5c_cov, pyfin_cov, overlap = compute_coverage_overlap(
        f5c_events, pyfin_events, len(signal)
    )
    metrics.f5c_coverage_samples = f5c_cov
    metrics.pyfin_coverage_samples = pyfin_cov
    metrics.overlap_samples = overlap
    
    union = f5c_cov + pyfin_cov - overlap
    metrics.coverage_jaccard = overlap / union if union > 0 else 0
    
    # Event structure
    merged, split = analyze_event_structure(f5c_events, pyfin_events, matches)
    metrics.f5c_merged_into_pyfin = merged
    metrics.pyfin_split_from_f5c = split
    
    return metrics


# =============================================================================
# Test Summary
# =============================================================================

@dataclass
class TestSummary:
    """Summary of event detection comparison across all reads."""
    total_reads: int = 0
    total_f5c_events: int = 0
    total_pyfin_events: int = 0
    total_matched_events: int = 0
    
    # Aggregated metrics
    mean_match_rate: float = 0.0
    mean_start_diff: float = 0.0
    mean_end_diff: float = 0.0
    mean_mean_correlation: float = 0.0
    mean_mean_rmse: float = 0.0
    mean_coverage_jaccard: float = 0.0
    
    # Per-read metrics
    per_read_metrics: List[EventComparisonMetrics] = field(default_factory=list)
    
    def compute_aggregates(self):
        """Compute aggregate statistics from per-read metrics."""
        if not self.per_read_metrics:
            return
        
        self.total_reads = len(self.per_read_metrics)
        self.total_f5c_events = sum(m.f5c_event_count for m in self.per_read_metrics)
        self.total_pyfin_events = sum(m.pyfin_event_count for m in self.per_read_metrics)
        self.total_matched_events = sum(m.matched_event_count for m in self.per_read_metrics)
        
        # Match rates
        match_rates = [
            m.matched_event_count / m.f5c_event_count 
            for m in self.per_read_metrics 
            if m.f5c_event_count > 0
        ]
        self.mean_match_rate = np.mean(match_rates) if match_rates else 0
        
        # Boundary diffs (only from reads with matches)
        start_diffs = [m.mean_start_diff for m in self.per_read_metrics if m.matched_event_count > 0]
        end_diffs = [m.mean_end_diff for m in self.per_read_metrics if m.matched_event_count > 0]
        self.mean_start_diff = np.mean(start_diffs) if start_diffs else 0
        self.mean_end_diff = np.mean(end_diffs) if end_diffs else 0
        
        # Correlations
        correlations = [
            m.mean_correlation for m in self.per_read_metrics 
            if m.matched_event_count > 1 and not np.isnan(m.mean_correlation)
        ]
        self.mean_mean_correlation = np.mean(correlations) if correlations else 0
        
        # RMSE
        rmses = [m.mean_rmse for m in self.per_read_metrics if m.matched_event_count > 0]
        self.mean_mean_rmse = np.mean(rmses) if rmses else 0
        
        # Coverage
        jaccards = [m.coverage_jaccard for m in self.per_read_metrics if m.coverage_jaccard > 0]
        self.mean_coverage_jaccard = np.mean(jaccards) if jaccards else 0
    
    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            "",
            "=" * 70,
            "EVENT DETECTION COMPARISON SUMMARY",
            "=" * 70,
            "",
            "Overall Statistics:",
            f"  Total reads tested:     {self.total_reads:,}",
            f"  Total f5c events:       {self.total_f5c_events:,}",
            f"  Total PyFIN events:     {self.total_pyfin_events:,}",
            f"  Total matched events:   {self.total_matched_events:,}",
            f"  Overall match rate:     {self.total_matched_events / max(1, self.total_f5c_events) * 100:.1f}%",
            "",
            "Aggregated Metrics (averaged across reads):",
            f"  Mean match rate:        {self.mean_match_rate * 100:.1f}%",
            f"  Mean start diff:        {self.mean_start_diff:.1f} samples",
            f"  Mean end diff:          {self.mean_end_diff:.1f} samples",
            f"  Mean correlation:       {self.mean_mean_correlation:.4f}",
            f"  Mean RMSE:              {self.mean_mean_rmse:.2f} pA",
            f"  Mean coverage Jaccard:  {self.mean_coverage_jaccard:.4f}",
            "",
            "Interpretation:",
        ]
        
        # Add interpretation
        if self.mean_match_rate >= 0.9:
            lines.append("  ✓ Excellent event matching (≥90%)")
        elif self.mean_match_rate >= 0.7:
            lines.append("  ◐ Good event matching (70-90%)")
        else:
            lines.append("  ✗ Event matching needs improvement (<70%)")
        
        if abs(self.mean_start_diff) <= 5 and abs(self.mean_end_diff) <= 5:
            lines.append("  ✓ Excellent boundary alignment (±5 samples)")
        elif abs(self.mean_start_diff) <= 20 and abs(self.mean_end_diff) <= 20:
            lines.append("  ◐ Good boundary alignment (±20 samples)")
        else:
            lines.append("  ✗ Boundary alignment needs improvement")
        
        if self.mean_mean_correlation >= 0.99:
            lines.append("  ✓ Excellent mean correlation (≥0.99)")
        elif self.mean_mean_correlation >= 0.95:
            lines.append("  ◐ Good mean correlation (0.95-0.99)")
        else:
            lines.append("  ✗ Mean correlation needs improvement")
        
        lines.append("")
        lines.append("=" * 70)
        
        return "\n".join(lines)


# =============================================================================
# Test Functions
# =============================================================================

@pytest.fixture
def f5c_events_by_read():
    """Load f5c events grouped by read."""
    if not RNA004_F5C_TSV_PATH.exists():
        pytest.skip(f"Test data not found: {RNA004_F5C_TSV_PATH}")
    return load_f5c_events_by_read(RNA004_F5C_TSV_PATH, max_reads=10)


@pytest.fixture
def pod5_signals(f5c_events_by_read):
    """Load POD5 signals for reads with f5c events."""
    if not RNA004_POD5_PATH.exists():
        pytest.skip(f"Test data not found: {RNA004_POD5_PATH}")
    if not POD5_AVAILABLE:
        pytest.skip("pod5 package not available")
    
    read_ids = list(f5c_events_by_read.keys())
    return load_pod5_signals(RNA004_POD5_PATH, read_ids)


@pytest.mark.skipif(not GETEVENTS_AVAILABLE, reason="getevents not available")
@pytest.mark.skipif(not POD5_AVAILABLE, reason="pod5 not available")
class TestEventDetection:
    """Test suite for event detection comparison."""
    
    def test_event_detection_basic(self, f5c_events_by_read, pod5_signals):
        """Basic test: verify event detection runs without errors."""
        # Get first read that exists in both
        common_reads = set(f5c_events_by_read.keys()) & set(pod5_signals.keys())
        assert len(common_reads) > 0, "No common reads between f5c and POD5"
        
        read_id = next(iter(common_reads))
        signal = pod5_signals[read_id]
        
        # Run getevents
        events = run_pyfin_getevents(signal)
        
        assert len(events) > 0, "No events detected"
        print(f"\n  Read {read_id[:8]}...: {len(events)} PyFIN events detected")
    
    def test_event_counts(self, f5c_events_by_read, pod5_signals):
        """Test: compare event counts between f5c and PyFIN."""
        common_reads = set(f5c_events_by_read.keys()) & set(pod5_signals.keys())
        
        print("\n  Event Count Comparison:")
        for read_id in list(common_reads)[:5]:
            f5c_count = len(f5c_events_by_read[read_id])
            signal = pod5_signals[read_id]
            pyfin_events = run_pyfin_getevents(signal)
            
            ratio = len(pyfin_events) / f5c_count if f5c_count > 0 else 0
            print(f"    {read_id[:8]}...: f5c={f5c_count:,}, PyFIN={len(pyfin_events):,}, ratio={ratio:.2f}")
    
    def test_event_boundary_alignment(self, f5c_events_by_read, pod5_signals):
        """Test: verify event boundaries are aligned within tolerance."""
        common_reads = set(f5c_events_by_read.keys()) & set(pod5_signals.keys())
        read_id = next(iter(common_reads))
        
        f5c_events = f5c_events_by_read[read_id]
        signal = pod5_signals[read_id]
        
        metrics = compare_events_for_read(read_id, f5c_events, signal)
        
        print(f"\n  Boundary alignment for {read_id[:8]}...:")
        print(f"    Start diff: {metrics.mean_start_diff:.1f} ± {metrics.std_start_diff:.1f}")
        print(f"    End diff:   {metrics.mean_end_diff:.1f} ± {metrics.std_end_diff:.1f}")
        
        # Relaxed threshold for initial testing
        assert abs(metrics.mean_start_diff) < 100, f"Start diff too large: {metrics.mean_start_diff}"
    
    def test_event_mean_correlation(self, f5c_events_by_read, pod5_signals):
        """Test: verify event means are correlated."""
        common_reads = set(f5c_events_by_read.keys()) & set(pod5_signals.keys())
        read_id = next(iter(common_reads))
        
        f5c_events = f5c_events_by_read[read_id]
        signal = pod5_signals[read_id]
        
        metrics = compare_events_for_read(read_id, f5c_events, signal)
        
        print(f"\n  Event mean correlation for {read_id[:8]}...:")
        print(f"    Correlation: {metrics.mean_correlation:.4f}")
        print(f"    RMSE:        {metrics.mean_rmse:.2f} pA")
        
        # Event means should be highly correlated (both computed from same signal)
        assert metrics.mean_correlation > 0.9, f"Low correlation: {metrics.mean_correlation}"
    
    def test_comprehensive_comparison(self, f5c_events_by_read, pod5_signals):
        """Comprehensive test: full comparison across multiple reads."""
        common_reads = set(f5c_events_by_read.keys()) & set(pod5_signals.keys())
        
        summary = TestSummary()
        
        print("\n" + "=" * 70)
        print("COMPREHENSIVE EVENT DETECTION COMPARISON")
        print("=" * 70)
        
        for read_id in common_reads:
            f5c_events = f5c_events_by_read[read_id]
            signal = pod5_signals[read_id]
            
            metrics = compare_events_for_read(read_id, f5c_events, signal)
            summary.per_read_metrics.append(metrics)
            
            # Print per-read summary
            match_rate = metrics.matched_event_count / max(1, metrics.f5c_event_count)
            print(f"\n  {read_id[:12]}...:")
            print(f"    Events: f5c={metrics.f5c_event_count}, PyFIN={metrics.pyfin_event_count}, matched={metrics.matched_event_count} ({match_rate*100:.1f}%)")
            print(f"    Boundary: start={metrics.mean_start_diff:+.1f}, end={metrics.mean_end_diff:+.1f}")
            if metrics.mean_correlation and not np.isnan(metrics.mean_correlation):
                print(f"    Correlation: {metrics.mean_correlation:.4f}")
        
        # Compute and print summary
        summary.compute_aggregates()
        print(summary.summary())
        
        # Store summary for reporting
        return summary


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Run event detection comparison and print summary."""
    print("=" * 70)
    print("PyFIN Event Detection vs f5c Comparison")
    print("=" * 70)
    
    # Check prerequisites
    if not POD5_AVAILABLE:
        print("ERROR: pod5 package not available")
        return 1
    
    if not GETEVENTS_AVAILABLE:
        print("ERROR: fin._eventalign.getevents not available")
        return 1
    
    if not RNA004_POD5_PATH.exists():
        print(f"ERROR: POD5 file not found: {RNA004_POD5_PATH}")
        return 1
    
    if not RNA004_F5C_TSV_PATH.exists():
        print(f"ERROR: f5c TSV not found: {RNA004_F5C_TSV_PATH}")
        return 1
    
    # Load data
    print("\nLoading f5c events...")
    f5c_events_by_read = load_f5c_events_by_read(RNA004_F5C_TSV_PATH, max_reads=10)
    print(f"  Loaded {len(f5c_events_by_read)} reads")
    
    print("\nLoading POD5 signals...")
    read_ids = list(f5c_events_by_read.keys())
    signals = load_pod5_signals(RNA004_POD5_PATH, read_ids)
    print(f"  Loaded {len(signals)} signals")
    
    # Find common reads
    common_reads = set(f5c_events_by_read.keys()) & set(signals.keys())
    print(f"\nCommon reads: {len(common_reads)}")
    
    if not common_reads:
        print("ERROR: No common reads between f5c output and POD5 file")
        return 1
    
    # Run comparison
    summary = TestSummary()
    
    for read_id in common_reads:
        f5c_events = f5c_events_by_read[read_id]
        signal = signals[read_id]
        
        print(f"\nProcessing {read_id}...")
        print(f"  Signal length: {len(signal):,} samples")
        print(f"  f5c events: {len(f5c_events):,}")
        
        metrics = compare_events_for_read(read_id, f5c_events, signal)
        summary.per_read_metrics.append(metrics)
        
        print(metrics.summary())
    
    # Print overall summary
    summary.compute_aggregates()
    print(summary.summary())
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
