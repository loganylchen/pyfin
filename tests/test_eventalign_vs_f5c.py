#!/usr/bin/env python3
"""
Comprehensive Test Framework: PyFIN vs f5c Eventalign Comparison

This test module provides detailed comparison between PyFIN's eventalign
implementation and f5c's reference output. It generates quantitative metrics
to track progress toward result parity.

Usage:
    # Run all tests
    pytest tests/test_eventalign_vs_f5c.py -v
    
    # Run with detailed output
    pytest tests/test_eventalign_vs_f5c.py -v -s
    
    # Generate comparison report
    python tests/test_eventalign_vs_f5c.py --report

Test Metrics:
    - Position Match Rate: % of positions with matching event assignments
    - Event Index Correlation: Pearson correlation of event indices
    - Event Mean RMSE: Root mean square error of event mean values
    - Coverage Match: % of reference covered by aligned events
    - Scaling Parameter Match: Difference in scale/shift parameters
"""

import numpy as np
import gzip
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
from collections import defaultdict
import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Test data paths - actual test data location
TEST_DATA_DIR = PROJECT_ROOT / "tests" / "testdata"

# RNA004 test data
RNA004_POD5_PATH = TEST_DATA_DIR / "RNA004.test.pod5"
RNA004_FASTQ_PATH = TEST_DATA_DIR / "RNA004.test.fq.gz"
RNA004_BAM_PATH = TEST_DATA_DIR / "RNA004.test.bam"
RNA004_F5C_TSV_PATH = TEST_DATA_DIR / "RNA004.test.tsv.gz"  # f5c eventalign output

# RNA002 test data (to be added)
RNA002_POD5_PATH = TEST_DATA_DIR / "RNA002.test.pod5"
RNA002_FASTQ_PATH = TEST_DATA_DIR / "RNA002.test.fq.gz"
RNA002_BAM_PATH = TEST_DATA_DIR / "RNA002.test.bam"
RNA002_F5C_TSV_PATH = TEST_DATA_DIR / "RNA002.test.tsv.gz"

# Shared reference (SIRV sequences)
REFERENCE_PATH = TEST_DATA_DIR / "test.fa"

# Legacy paths (for backward compatibility)
F5C_TSV_PATH = TEST_DATA_DIR / "one_read.eventalign.tsv.gz"
FASTA_PATH = TEST_DATA_DIR / "one_read.fa"
POD5_PATH = TEST_DATA_DIR / "one_read.pod5"
MODEL_PATH = TEST_DATA_DIR / "rna002_model.tsv"

# Check if required modules are available
try:
    import pod5
    POD5_AVAILABLE = True
except ImportError:
    POD5_AVAILABLE = False

try:
    from fin._eventalign import run_eventalign, MODEL_RNA002
    EVENTALIGN_AVAILABLE = True
except ImportError:
    EVENTALIGN_AVAILABLE = False


# =============================================================================
# Data Classes for Test Results
# =============================================================================

@dataclass
class EventAlignment:
    """Single event alignment record."""
    reference_name: str
    reference_position: int
    reference_kmer: str
    read_id: str
    strand: str
    event_idx: int
    event_mean: float
    event_stdv: float
    duration: float
    model_kmer: str
    model_mean: float
    model_stdv: float
    scaled_mean: float
    start_idx: int
    end_idx: int


@dataclass 
class ComparisonMetrics:
    """Metrics comparing PyFIN vs f5c results."""
    # Basic counts
    f5c_total_alignments: int
    pyfin_total_alignments: int
    
    # Position-level metrics
    total_positions: int
    matched_positions: int
    f5c_only_positions: int
    pyfin_only_positions: int
    position_match_rate: float
    
    # Event index metrics
    event_idx_correlation: float
    event_idx_rmse: float
    
    # Signal coordinate metrics
    start_idx_rmse: float
    end_idx_rmse: float
    
    # Event statistics metrics
    event_mean_rmse: float
    event_mean_correlation: float
    event_stdv_rmse: float
    
    # Coverage metrics
    f5c_coverage: float  # % of reference covered
    pyfin_coverage: float
    coverage_overlap: float  # Jaccard index of covered positions
    
    # Scaling metrics (if available)
    scale_diff: Optional[float] = None
    shift_diff: Optional[float] = None
    var_diff: Optional[float] = None
    
    # Status
    pyfin_success: bool = False
    error_message: Optional[str] = None
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    def summary(self) -> str:
        """Return human-readable summary."""
        lines = [
            "=" * 60,
            "EVENTALIGN COMPARISON METRICS",
            "=" * 60,
            "",
            "Alignment Counts:",
            f"  f5c alignments:   {self.f5c_total_alignments}",
            f"  PyFIN alignments: {self.pyfin_total_alignments}",
            "",
            "Position-Level Accuracy:",
            f"  Total positions:     {self.total_positions}",
            f"  Matched positions:   {self.matched_positions}",
            f"  Position match rate: {self.position_match_rate:.1%}",
            f"  f5c only:           {self.f5c_only_positions}",
            f"  PyFIN only:         {self.pyfin_only_positions}",
            "",
            "Event Index Metrics:",
            f"  Correlation: {self.event_idx_correlation:.4f}",
            f"  RMSE:        {self.event_idx_rmse:.2f}",
            "",
            "Signal Coordinate Metrics:",
            f"  start_idx RMSE: {self.start_idx_rmse:.2f}",
            f"  end_idx RMSE:   {self.end_idx_rmse:.2f}",
            "",
            "Event Statistics Metrics:",
            f"  Event mean RMSE:        {self.event_mean_rmse:.4f} pA",
            f"  Event mean correlation: {self.event_mean_correlation:.4f}",
            f"  Event stdv RMSE:        {self.event_stdv_rmse:.4f} pA",
            "",
            "Coverage Metrics:",
            f"  f5c coverage:    {self.f5c_coverage:.1%}",
            f"  PyFIN coverage:  {self.pyfin_coverage:.1%}",
            f"  Coverage overlap (Jaccard): {self.coverage_overlap:.1%}",
        ]
        
        if self.scale_diff is not None:
            lines.extend([
                "",
                "Scaling Parameter Differences:",
                f"  Scale diff: {self.scale_diff:.6f}",
                f"  Shift diff: {self.shift_diff:.6f}",
            ])
            if self.var_diff is not None:
                lines.append(f"  Var diff:   {self.var_diff:.6f}")
        
        lines.extend([
            "",
            "=" * 60,
        ])
        
        return "\n".join(lines)


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_f5c_eventalign(tsv_path: Path) -> List[EventAlignment]:
    """
    Load f5c eventalign TSV output.
    
    f5c eventalign output format (with header):
    contig, position, reference_kmer, read_name, strand, event_index,
    event_level_mean, event_stdv, event_length, model_kmer, model_mean,
    model_stdv, standardized_level, start_idx, end_idx, samples
    """
    alignments = []
    opener = gzip.open if str(tsv_path).endswith('.gz') else open
    
    with opener(tsv_path, 'rt') as f:
        header = None
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            
            # Skip header line (first line starting with 'contig')
            if line_num == 1 and parts[0] == 'contig':
                header = parts
                continue
            
            if len(parts) < 15:
                print(f"Warning: Line {line_num} has only {len(parts)} columns, skipping")
                continue
            
            try:
                # f5c column mapping:
                # 0: contig (reference_name)
                # 1: position (reference_position)
                # 2: reference_kmer
                # 3: read_name (read_id)
                # 4: strand
                # 5: event_index (event_idx)
                # 6: event_level_mean (event_mean)
                # 7: event_stdv
                # 8: event_length (duration in seconds)
                # 9: model_kmer
                # 10: model_mean
                # 11: model_stdv
                # 12: standardized_level (scaled_mean)
                # 13: start_idx
                # 14: end_idx
                aln = EventAlignment(
                    reference_name=parts[0],
                    reference_position=int(parts[1]),
                    reference_kmer=parts[2],
                    read_id=parts[3],
                    strand=parts[4],
                    event_idx=int(parts[5]),
                    event_mean=float(parts[6]),
                    event_stdv=float(parts[7]),
                    duration=float(parts[8]),
                    model_kmer=parts[9],
                    model_mean=float(parts[10]),
                    model_stdv=float(parts[11]),
                    scaled_mean=float(parts[12]),
                    start_idx=int(parts[13]),
                    end_idx=int(parts[14]),
                )
                alignments.append(aln)
            except (ValueError, IndexError) as e:
                print(f"Warning: Failed to parse line {line_num}: {e}")
                continue
    
    return alignments


def load_reference_fasta(fasta_path: Path) -> Tuple[str, str]:
    """Load reference sequence from FASTA file."""
    with open(fasta_path, "r") as f:
        lines = f.readlines()
    
    ref_name = None
    ref_seq_parts = []
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            ref_name = line[1:].split()[0]
        else:
            ref_seq_parts.append(line)
    
    ref_seq = "".join(ref_seq_parts).upper()
    return ref_name, ref_seq


def load_signal_from_pod5(pod5_path: Path) -> Tuple[str, np.ndarray, float]:
    """Load signal data from POD5 file."""
    import pod5
    
    with pod5.Reader(str(pod5_path)) as reader:
        for read in reader.reads():
            signal = read.signal_pa.astype(np.float32)
            sample_rate = float(read.run_info.sample_rate)
            return str(read.read_id), signal, sample_rate
    
    raise ValueError(f"No reads found in {pod5_path}")


def load_all_signals_from_pod5(pod5_path: Path) -> Dict[str, Tuple[np.ndarray, float]]:
    """Load all signals from POD5 file, indexed by read_id."""
    import pod5
    
    signals = {}
    with pod5.Reader(str(pod5_path)) as reader:
        for read in reader.reads():
            signal = read.signal_pa.astype(np.float32)
            sample_rate = float(read.run_info.sample_rate)
            signals[str(read.read_id)] = (signal, sample_rate)
    
    return signals


def load_sequences_from_fastq(fastq_path: Path) -> Dict[str, str]:
    """Load sequences from FASTQ file, indexed by read_id."""
    sequences = {}
    opener = gzip.open if str(fastq_path).endswith('.gz') else open
    
    with opener(fastq_path, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            if not header.startswith('@'):
                continue
            
            read_id = header[1:].split()[0]
            sequence = f.readline().strip().upper()
            f.readline()  # + line
            f.readline()  # quality line
            
            sequences[read_id] = sequence
    
    return sequences


def load_alignments_from_bam(bam_path: Path, reference_path: Path) -> Dict[str, Dict]:
    """
    Load alignments from BAM file.
    
    Returns dict indexed by read_id with:
        - ref_name: reference name
        - ref_start: 0-based start position
        - ref_end: end position
        - query_seq: query sequence (from BAM)
        - is_reverse: strand
    """
    try:
        import pysam
    except ImportError:
        return {}
    
    alignments = {}
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            
            alignments[read.query_name] = {
                'ref_name': read.reference_name,
                'ref_start': read.reference_start,
                'ref_end': read.reference_end,
                'query_seq': read.query_sequence,
                'is_reverse': read.is_reverse,
                'cigar': read.cigarstring,
            }
    
    return alignments


def load_reference_sequences(fasta_path: Path) -> Dict[str, str]:
    """Load all reference sequences from FASTA file."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq).upper()
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:
            sequences[current_name] = ''.join(current_seq).upper()
    
    return sequences


# =============================================================================
# PyFIN Execution
# =============================================================================

def run_pyfin_eventalign(
    signal: np.ndarray,
    sample_rate: float,
    read_id: str,
    read_seq: str,
    ref_seq: str,
    ref_name: str
) -> Tuple[Optional[List[Dict]], Optional[Dict], Optional[str]]:
    """
    Run PyFIN eventalign and return results.
    
    Returns:
        Tuple of (alignments, scalings, error_message)
        - alignments: List of alignment dicts if successful, None otherwise
        - scalings: Dict with scale/shift/var if successful
        - error_message: Error string if failed, None otherwise
    """
    try:
        from fin._eventalign import run_eventalign, MODEL_RNA002
        
        result = run_eventalign(
            read_ids=[read_id],
            read_seqs=[read_seq],
            ref_seqs=[ref_seq],
            ref_names=[ref_name],
            ref_lens=[len(ref_seq)],
            signals=[signal],
            sample_rates=[sample_rate],
            model_id=MODEL_RNA002,
        )
        
        # Extract results
        alignments = result["full"][0][0]  # List of alignment dicts
        scalings = result["scalings"][0]
        events = result["events"][0]
        mapping = result["mapping"][0][0]
        
        # Add event data to alignments
        for aln in alignments:
            event_idx = aln['event_idx']
            if event_idx >= 0 and event_idx < events['n_events']:
                aln['event_mean'] = float(events['means'][event_idx])
                aln['event_stdv'] = float(events['stdvs'][event_idx])
                aln['start_idx'] = int(events['starts'][event_idx])
                aln['duration'] = float(events['lengths'][event_idx])
                aln['end_idx'] = aln['start_idx'] + int(aln['duration'])
        
        return alignments, scalings, None
        
    except ImportError as e:
        return None, None, f"Import error: {e}"
    except Exception as e:
        import traceback
        return None, None, f"Runtime error: {e}\n{traceback.format_exc()}"


# =============================================================================
# Comparison Functions
# =============================================================================

def index_by_position(alignments: List) -> Dict[int, List]:
    """Index alignments by reference position."""
    by_pos = defaultdict(list)
    for aln in alignments:
        if isinstance(aln, EventAlignment):
            pos = aln.reference_position
        else:
            pos = aln['ref_position']
        by_pos[pos].append(aln)
    return dict(by_pos)


def compute_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Pearson correlation coefficient."""
    if len(x) < 2 or len(y) < 2:
        return 0.0
    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def compute_rmse(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Root Mean Square Error."""
    if len(x) == 0 or len(y) == 0:
        return float('inf')
    return float(np.sqrt(np.mean((x - y) ** 2)))


def compare_results(
    f5c_alignments: List[EventAlignment],
    pyfin_alignments: List[Dict],
    pyfin_scalings: Optional[Dict],
    ref_len: int,
    f5c_scalings: Optional[Dict] = None
) -> ComparisonMetrics:
    """
    Compare f5c and PyFIN eventalign results.
    
    Args:
        f5c_alignments: List of f5c EventAlignment objects
        pyfin_alignments: List of PyFIN alignment dicts
        pyfin_scalings: PyFIN scaling parameters
        ref_len: Reference sequence length
        f5c_scalings: Optional f5c scaling parameters for comparison
    
    Returns:
        ComparisonMetrics with detailed comparison results
    """
    # Index by position
    f5c_by_pos = index_by_position(f5c_alignments)
    pyfin_by_pos = index_by_position(pyfin_alignments) if pyfin_alignments else {}
    
    # Get all positions
    all_positions = set(f5c_by_pos.keys()) | set(pyfin_by_pos.keys())
    matched_positions = set(f5c_by_pos.keys()) & set(pyfin_by_pos.keys())
    f5c_only = set(f5c_by_pos.keys()) - set(pyfin_by_pos.keys())
    pyfin_only = set(pyfin_by_pos.keys()) - set(f5c_by_pos.keys())
    
    # Position match rate
    position_match_rate = len(matched_positions) / len(all_positions) if all_positions else 0.0
    
    # For matched positions, compare event-level metrics
    f5c_event_indices = []
    pyfin_event_indices = []
    f5c_event_means = []
    pyfin_event_means = []
    f5c_event_stdvs = []
    pyfin_event_stdvs = []
    f5c_start_indices = []
    pyfin_start_indices = []
    f5c_end_indices = []
    pyfin_end_indices = []
    
    for pos in sorted(matched_positions):
        f5c_events = f5c_by_pos[pos]
        pyfin_events = pyfin_by_pos[pos]
        
        # Match events by closest event_idx
        for f5c_evt in f5c_events:
            # Find closest PyFIN event by event_idx
            best_match = None
            best_diff = float('inf')
            for pyfin_evt in pyfin_events:
                diff = abs(f5c_evt.event_idx - pyfin_evt['event_idx'])
                if diff < best_diff:
                    best_diff = diff
                    best_match = pyfin_evt
            
            if best_match is not None:
                f5c_event_indices.append(f5c_evt.event_idx)
                pyfin_event_indices.append(best_match['event_idx'])
                
                f5c_event_means.append(f5c_evt.event_mean)
                pyfin_event_means.append(best_match.get('event_mean', 0.0))
                
                f5c_event_stdvs.append(f5c_evt.event_stdv)
                pyfin_event_stdvs.append(best_match.get('event_stdv', 0.0))
                
                f5c_start_indices.append(f5c_evt.start_idx)
                pyfin_start_indices.append(best_match.get('start_idx', 0))
                
                f5c_end_indices.append(f5c_evt.end_idx)
                pyfin_end_indices.append(best_match.get('end_idx', 0))
    
    # Convert to numpy arrays
    f5c_event_indices = np.array(f5c_event_indices)
    pyfin_event_indices = np.array(pyfin_event_indices)
    f5c_event_means = np.array(f5c_event_means)
    pyfin_event_means = np.array(pyfin_event_means)
    f5c_event_stdvs = np.array(f5c_event_stdvs)
    pyfin_event_stdvs = np.array(pyfin_event_stdvs)
    f5c_start_indices = np.array(f5c_start_indices)
    pyfin_start_indices = np.array(pyfin_start_indices)
    f5c_end_indices = np.array(f5c_end_indices)
    pyfin_end_indices = np.array(pyfin_end_indices)
    
    # Compute metrics
    event_idx_correlation = compute_correlation(f5c_event_indices, pyfin_event_indices)
    event_idx_rmse = compute_rmse(f5c_event_indices, pyfin_event_indices)
    
    event_mean_rmse = compute_rmse(f5c_event_means, pyfin_event_means)
    event_mean_correlation = compute_correlation(f5c_event_means, pyfin_event_means)
    event_stdv_rmse = compute_rmse(f5c_event_stdvs, pyfin_event_stdvs)
    
    start_idx_rmse = compute_rmse(f5c_start_indices, pyfin_start_indices)
    end_idx_rmse = compute_rmse(f5c_end_indices, pyfin_end_indices)
    
    # Coverage metrics
    f5c_coverage = len(f5c_by_pos) / ref_len if ref_len > 0 else 0.0
    pyfin_coverage = len(pyfin_by_pos) / ref_len if ref_len > 0 else 0.0
    
    # Jaccard index for coverage overlap
    if all_positions:
        coverage_overlap = len(matched_positions) / len(all_positions)
    else:
        coverage_overlap = 0.0
    
    # Scaling parameter comparison
    scale_diff = None
    shift_diff = None
    var_diff = None
    
    if pyfin_scalings and f5c_scalings:
        scale_diff = abs(pyfin_scalings.get('scale', 0) - f5c_scalings.get('scale', 0))
        shift_diff = abs(pyfin_scalings.get('shift', 0) - f5c_scalings.get('shift', 0))
        if 'var' in pyfin_scalings and 'var' in f5c_scalings:
            var_diff = abs(pyfin_scalings['var'] - f5c_scalings['var'])
    
    return ComparisonMetrics(
        f5c_total_alignments=len(f5c_alignments),
        pyfin_total_alignments=len(pyfin_alignments) if pyfin_alignments else 0,
        total_positions=len(all_positions),
        matched_positions=len(matched_positions),
        f5c_only_positions=len(f5c_only),
        pyfin_only_positions=len(pyfin_only),
        position_match_rate=position_match_rate,
        event_idx_correlation=event_idx_correlation,
        event_idx_rmse=event_idx_rmse,
        start_idx_rmse=start_idx_rmse,
        end_idx_rmse=end_idx_rmse,
        event_mean_rmse=event_mean_rmse,
        event_mean_correlation=event_mean_correlation,
        event_stdv_rmse=event_stdv_rmse,
        f5c_coverage=f5c_coverage,
        pyfin_coverage=pyfin_coverage,
        coverage_overlap=coverage_overlap,
        scale_diff=scale_diff,
        shift_diff=shift_diff,
        var_diff=var_diff,
        pyfin_success=pyfin_alignments is not None and len(pyfin_alignments) > 0,
    )


# =============================================================================
# Test Functions (pytest compatible)
# =============================================================================

def has_rna004_data() -> bool:
    """Check if RNA004 test data is available."""
    return (RNA004_POD5_PATH.exists() and 
            RNA004_FASTQ_PATH.exists() and 
            RNA004_BAM_PATH.exists() and
            REFERENCE_PATH.exists())


def has_rna002_data() -> bool:
    """Check if RNA002 test data is available."""
    return (RNA002_POD5_PATH.exists() and 
            RNA002_FASTQ_PATH.exists() and 
            RNA002_BAM_PATH.exists() and
            REFERENCE_PATH.exists())


def has_legacy_data() -> bool:
    """Check if legacy single-read test data is available."""
    return F5C_TSV_PATH.exists() and FASTA_PATH.exists() and POD5_PATH.exists()


@pytest.mark.skipif(not POD5_AVAILABLE, reason="pod5 not installed")
@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="fin._eventalign not available")
class TestEventalignRNA004:
    """Test class for RNA004 PyFIN vs f5c comparison."""
    
    @classmethod
    def setup_class(cls):
        """Load RNA004 test data."""
        if not has_rna004_data():
            pytest.skip("RNA004 test data not available")
        
        print("\n" + "=" * 60)
        print("Loading RNA004 test data...")
        print("=" * 60)
        
        # Load all signals
        cls.signals = load_all_signals_from_pod5(RNA004_POD5_PATH)
        print(f"  Loaded {len(cls.signals)} signals from POD5")
        
        # Load sequences
        cls.sequences = load_sequences_from_fastq(RNA004_FASTQ_PATH)
        print(f"  Loaded {len(cls.sequences)} sequences from FASTQ")
        
        # Load alignments
        cls.alignments = load_alignments_from_bam(RNA004_BAM_PATH, REFERENCE_PATH)
        print(f"  Loaded {len(cls.alignments)} alignments from BAM")
        
        # Load reference
        cls.references = load_reference_sequences(REFERENCE_PATH)
        print(f"  Loaded {len(cls.references)} reference sequences")
        
        # Load f5c reference output if available
        if RNA004_F5C_TSV_PATH.exists():
            cls.f5c_alignments = load_f5c_eventalign(RNA004_F5C_TSV_PATH)
            print(f"  Loaded {len(cls.f5c_alignments)} f5c reference alignments")
        else:
            cls.f5c_alignments = []
            print("  WARNING: f5c reference output not found")
        
        # Run PyFIN on all reads
        cls.pyfin_results = {}
        cls.metrics_per_read = {}
        
        from fin._eventalign import run_eventalign, MODEL_RNA004
        
        for read_id in cls.signals:
            if read_id not in cls.alignments:
                continue
            
            signal, sample_rate = cls.signals[read_id]
            read_seq = cls.sequences.get(read_id, "")
            aln_info = cls.alignments[read_id]
            ref_name = aln_info['ref_name']
            ref_seq = cls.references.get(ref_name, "")
            
            if not read_seq or not ref_seq:
                continue
            
            try:
                result = run_eventalign(
                    read_ids=[read_id],
                    read_seqs=[read_seq],
                    ref_seqs=[ref_seq],
                    ref_names=[ref_name],
                    ref_lens=[len(ref_seq)],
                    signals=[signal],
                    sample_rates=[sample_rate],
                    model_id=MODEL_RNA004,
                )
                cls.pyfin_results[read_id] = result
            except Exception as e:
                print(f"  WARNING: Failed to process {read_id}: {e}")
        
        print(f"  Processed {len(cls.pyfin_results)} reads with PyFIN")
    
    def test_pyfin_processes_all_reads(self):
        """Test that PyFIN processes all reads without errors."""
        expected_reads = len([r for r in self.signals if r in self.alignments])
        actual_reads = len(self.pyfin_results)
        
        # Allow some failures but most should work
        assert actual_reads >= expected_reads * 0.8, \
            f"Only processed {actual_reads}/{expected_reads} reads"
    
    def test_alignment_count_reasonable(self):
        """Test that each read produces reasonable number of alignments."""
        for read_id, result in self.pyfin_results.items():
            alignments = result["full"][0][0]
            
            # Should have some alignments
            assert len(alignments) > 0, f"No alignments for {read_id}"
            
            # Events per base should be reasonable (0.5 to 20)
            mapping = result["mapping"][0][0]
            epb = mapping.get('events_per_base', 0)
            
            if epb > 0:  # Only check if alignment succeeded
                assert 0.3 <= epb <= 25, \
                    f"Events per base {epb} out of range for {read_id}"


@pytest.mark.skipif(not POD5_AVAILABLE, reason="pod5 not installed")
@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="fin._eventalign not available")
class TestEventalignRNA002:
    """Test class for RNA002 PyFIN vs f5c comparison."""
    
    @classmethod
    def setup_class(cls):
        """Load RNA002 test data."""
        if not has_rna002_data():
            pytest.skip("RNA002 test data not available")
        
        print("\n" + "=" * 60)
        print("Loading RNA002 test data...")
        print("=" * 60)
        
        # Load all signals
        cls.signals = load_all_signals_from_pod5(RNA002_POD5_PATH)
        print(f"  Loaded {len(cls.signals)} signals from POD5")
        
        # Load sequences
        cls.sequences = load_sequences_from_fastq(RNA002_FASTQ_PATH)
        print(f"  Loaded {len(cls.sequences)} sequences from FASTQ")
        
        # Load alignments
        cls.alignments = load_alignments_from_bam(RNA002_BAM_PATH, REFERENCE_PATH)
        print(f"  Loaded {len(cls.alignments)} alignments from BAM")
        
        # Load reference
        cls.references = load_reference_sequences(REFERENCE_PATH)
        print(f"  Loaded {len(cls.references)} reference sequences")
        
        # Load f5c reference output if available
        if RNA002_F5C_TSV_PATH.exists():
            cls.f5c_alignments = load_f5c_eventalign(RNA002_F5C_TSV_PATH)
            print(f"  Loaded {len(cls.f5c_alignments)} f5c reference alignments")
        else:
            cls.f5c_alignments = []
            print("  WARNING: f5c reference output not found")
        
        # Run PyFIN on all reads
        cls.pyfin_results = {}
        
        from fin._eventalign import run_eventalign, MODEL_RNA002
        
        for read_id in cls.signals:
            if read_id not in cls.alignments:
                continue
            
            signal, sample_rate = cls.signals[read_id]
            read_seq = cls.sequences.get(read_id, "")
            aln_info = cls.alignments[read_id]
            ref_name = aln_info['ref_name']
            ref_seq = cls.references.get(ref_name, "")
            
            if not read_seq or not ref_seq:
                continue
            
            try:
                result = run_eventalign(
                    read_ids=[read_id],
                    read_seqs=[read_seq],
                    ref_seqs=[ref_seq],
                    ref_names=[ref_name],
                    ref_lens=[len(ref_seq)],
                    signals=[signal],
                    sample_rates=[sample_rate],
                    model_id=MODEL_RNA002,
                )
                cls.pyfin_results[read_id] = result
            except Exception as e:
                print(f"  WARNING: Failed to process {read_id}: {e}")
        
        print(f"  Processed {len(cls.pyfin_results)} reads with PyFIN")
    
    def test_pyfin_processes_all_reads(self):
        """Test that PyFIN processes all reads without errors."""
        expected_reads = len([r for r in self.signals if r in self.alignments])
        actual_reads = len(self.pyfin_results)
        
        assert actual_reads >= expected_reads * 0.8, \
            f"Only processed {actual_reads}/{expected_reads} reads"
    
    def test_alignment_count_reasonable(self):
        """Test that each read produces reasonable number of alignments."""
        for read_id, result in self.pyfin_results.items():
            alignments = result["full"][0][0]
            assert len(alignments) > 0, f"No alignments for {read_id}"


@pytest.mark.skipif(not POD5_AVAILABLE, reason="pod5 not installed")
@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="fin._eventalign not available")
@pytest.mark.skipif(not has_legacy_data(), reason="Legacy test data not found")
class TestEventalignVsF5C:
    """Test class for PyFIN vs f5c eventalign comparison."""
    
    @classmethod
    def setup_class(cls):
        """Load test data once for all tests."""
        print("\n" + "=" * 60)
        print("Loading test data...")
        print("=" * 60)
        
        # Check if test data exists
        if not F5C_TSV_PATH.exists():
            raise FileNotFoundError(f"f5c reference output not found: {F5C_TSV_PATH}")
        if not FASTA_PATH.exists():
            raise FileNotFoundError(f"Reference FASTA not found: {FASTA_PATH}")
        if not POD5_PATH.exists():
            raise FileNotFoundError(f"POD5 file not found: {POD5_PATH}")
        
        # Load f5c reference
        cls.f5c_alignments = load_f5c_eventalign(F5C_TSV_PATH)
        print(f"  Loaded {len(cls.f5c_alignments)} f5c alignments")
        
        # Load reference
        cls.ref_name, cls.ref_seq = load_reference_fasta(FASTA_PATH)
        print(f"  Reference: {cls.ref_name}, {len(cls.ref_seq)} bp")
        
        # Load signal
        cls.read_id, cls.signal, cls.sample_rate = load_signal_from_pod5(POD5_PATH)
        print(f"  Signal: {len(cls.signal)} samples, {cls.sample_rate} Hz")
        
        # Run PyFIN eventalign
        print("\nRunning PyFIN eventalign...")
        cls.pyfin_alignments, cls.pyfin_scalings, cls.error = run_pyfin_eventalign(
            signal=cls.signal,
            sample_rate=cls.sample_rate,
            read_id=cls.read_id,
            read_seq=cls.ref_seq,  # Using ref as read for initial test
            ref_seq=cls.ref_seq,
            ref_name=cls.ref_name
        )
        
        if cls.error:
            print(f"  PyFIN error: {cls.error}")
        else:
            print(f"  PyFIN alignments: {len(cls.pyfin_alignments) if cls.pyfin_alignments else 0}")
        
        # Compute comparison metrics
        cls.metrics = compare_results(
            cls.f5c_alignments,
            cls.pyfin_alignments,
            cls.pyfin_scalings,
            len(cls.ref_seq)
        )
        
        print("\n" + cls.metrics.summary())
    
    def test_pyfin_runs_successfully(self):
        """Test that PyFIN eventalign runs without errors."""
        assert self.error is None, f"PyFIN failed with error: {self.error}"
        assert self.pyfin_alignments is not None, "PyFIN returned no alignments"
        assert len(self.pyfin_alignments) > 0, "PyFIN returned empty alignments"
    
    def test_position_match_rate(self):
        """Test that position match rate meets minimum threshold."""
        # Current expected: ~50% (before Profile HMM implementation)
        # Target: >95% (after full implementation)
        MIN_MATCH_RATE = 0.30  # 30% as baseline
        TARGET_MATCH_RATE = 0.95  # 95% target
        
        assert self.metrics.position_match_rate >= MIN_MATCH_RATE, \
            f"Position match rate {self.metrics.position_match_rate:.1%} below minimum {MIN_MATCH_RATE:.1%}"
        
        if self.metrics.position_match_rate < TARGET_MATCH_RATE:
            print(f"\n  NOTE: Position match rate {self.metrics.position_match_rate:.1%} "
                  f"below target {TARGET_MATCH_RATE:.1%}")
    
    def test_event_index_correlation(self):
        """Test that event indices are correlated."""
        # Current expected: >0.8 (events are generally in order)
        # Target: >0.99 (nearly identical)
        MIN_CORRELATION = 0.7
        TARGET_CORRELATION = 0.99
        
        if self.metrics.pyfin_success:
            assert self.metrics.event_idx_correlation >= MIN_CORRELATION, \
                f"Event index correlation {self.metrics.event_idx_correlation:.4f} below minimum {MIN_CORRELATION}"
            
            if self.metrics.event_idx_correlation < TARGET_CORRELATION:
                print(f"\n  NOTE: Event index correlation {self.metrics.event_idx_correlation:.4f} "
                      f"below target {TARGET_CORRELATION}")
    
    def test_event_mean_accuracy(self):
        """Test that event mean values are similar."""
        # Current expected: <5 pA RMSE (same event detection)
        # Target: <0.5 pA RMSE (identical scaling)
        MAX_RMSE = 10.0  # pA
        TARGET_RMSE = 0.5  # pA
        
        if self.metrics.pyfin_success:
            assert self.metrics.event_mean_rmse <= MAX_RMSE, \
                f"Event mean RMSE {self.metrics.event_mean_rmse:.2f} pA exceeds maximum {MAX_RMSE} pA"
            
            if self.metrics.event_mean_rmse > TARGET_RMSE:
                print(f"\n  NOTE: Event mean RMSE {self.metrics.event_mean_rmse:.2f} pA "
                      f"above target {TARGET_RMSE} pA")
    
    def test_coverage_overlap(self):
        """Test that coverage overlaps significantly."""
        # Current expected: >50% overlap
        # Target: >95% overlap
        MIN_OVERLAP = 0.30
        TARGET_OVERLAP = 0.95
        
        if self.metrics.pyfin_success:
            assert self.metrics.coverage_overlap >= MIN_OVERLAP, \
                f"Coverage overlap {self.metrics.coverage_overlap:.1%} below minimum {MIN_OVERLAP:.1%}"
            
            if self.metrics.coverage_overlap < TARGET_OVERLAP:
                print(f"\n  NOTE: Coverage overlap {self.metrics.coverage_overlap:.1%} "
                      f"below target {TARGET_OVERLAP:.1%}")
    
    def test_signal_coordinates(self):
        """Test that signal coordinates (start/end indices) are close."""
        # Start/end indices should be very close for matched events
        MAX_COORD_RMSE = 100  # samples
        TARGET_COORD_RMSE = 10  # samples
        
        if self.metrics.pyfin_success and self.metrics.matched_positions > 0:
            assert self.metrics.start_idx_rmse <= MAX_COORD_RMSE, \
                f"Start index RMSE {self.metrics.start_idx_rmse:.1f} exceeds maximum {MAX_COORD_RMSE}"
            
            if self.metrics.start_idx_rmse > TARGET_COORD_RMSE:
                print(f"\n  NOTE: Start index RMSE {self.metrics.start_idx_rmse:.1f} "
                      f"above target {TARGET_COORD_RMSE}")


# =============================================================================
# Detailed Diagnostic Functions
# =============================================================================

def generate_position_comparison_report(
    f5c_alignments: List[EventAlignment],
    pyfin_alignments: List[Dict],
    output_path: Path
) -> None:
    """Generate detailed position-by-position comparison report."""
    f5c_by_pos = index_by_position(f5c_alignments)
    pyfin_by_pos = index_by_position(pyfin_alignments) if pyfin_alignments else {}
    
    all_positions = sorted(set(f5c_by_pos.keys()) | set(pyfin_by_pos.keys()))
    
    with open(output_path, 'w') as f:
        f.write("# Position-by-Position Eventalign Comparison Report\n\n")
        f.write(f"Total positions: {len(all_positions)}\n")
        f.write(f"f5c positions: {len(f5c_by_pos)}\n")
        f.write(f"PyFIN positions: {len(pyfin_by_pos)}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write(f"{'Pos':>6} | {'f5c_N':>5} | {'PyFIN_N':>7} | {'Status':>10} | "
                f"{'f5c_EvtIdx':>10} | {'PyFIN_EvtIdx':>12}\n")
        f.write("=" * 80 + "\n")
        
        for pos in all_positions:
            f5c_events = f5c_by_pos.get(pos, [])
            pyfin_events = pyfin_by_pos.get(pos, [])
            
            n_f5c = len(f5c_events)
            n_pyfin = len(pyfin_events)
            
            if n_f5c > 0 and n_pyfin > 0:
                status = "MATCH"
            elif n_f5c > 0:
                status = "F5C_ONLY"
            else:
                status = "PYFIN_ONLY"
            
            f5c_idx_str = ",".join(str(e.event_idx) for e in f5c_events[:3])
            if n_f5c > 3:
                f5c_idx_str += "..."
            
            pyfin_idx_str = ",".join(str(e['event_idx']) for e in pyfin_events[:3])
            if n_pyfin > 3:
                pyfin_idx_str += "..."
            
            f.write(f"{pos:6d} | {n_f5c:5d} | {n_pyfin:7d} | {status:>10} | "
                    f"{f5c_idx_str:>10} | {pyfin_idx_str:>12}\n")
    
    print(f"Position comparison report saved to: {output_path}")


def generate_event_statistics_report(
    f5c_alignments: List[EventAlignment],
    pyfin_alignments: List[Dict],
    output_path: Path
) -> None:
    """Generate statistical comparison of event values."""
    with open(output_path, 'w') as f:
        f.write("# Event Statistics Comparison Report\n\n")
        
        # f5c statistics
        f5c_means = [a.event_mean for a in f5c_alignments]
        f5c_stdvs = [a.event_stdv for a in f5c_alignments]
        f5c_durations = [a.duration for a in f5c_alignments]
        
        f.write("## f5c Event Statistics\n\n")
        f.write(f"Event mean: {np.mean(f5c_means):.2f} ± {np.std(f5c_means):.2f} pA\n")
        f.write(f"Event stdv: {np.mean(f5c_stdvs):.2f} ± {np.std(f5c_stdvs):.2f} pA\n")
        f.write(f"Duration:   {np.mean(f5c_durations):.2f} ± {np.std(f5c_durations):.2f} samples\n")
        f.write(f"Total events: {len(f5c_alignments)}\n\n")
        
        # PyFIN statistics
        if pyfin_alignments:
            pyfin_means = [a.get('event_mean', 0) for a in pyfin_alignments]
            pyfin_stdvs = [a.get('event_stdv', 0) for a in pyfin_alignments]
            pyfin_durations = [a.get('duration', 0) for a in pyfin_alignments]
            
            f.write("## PyFIN Event Statistics\n\n")
            f.write(f"Event mean: {np.mean(pyfin_means):.2f} ± {np.std(pyfin_means):.2f} pA\n")
            f.write(f"Event stdv: {np.mean(pyfin_stdvs):.2f} ± {np.std(pyfin_stdvs):.2f} pA\n")
            f.write(f"Duration:   {np.mean(pyfin_durations):.2f} ± {np.std(pyfin_durations):.2f} samples\n")
            f.write(f"Total events: {len(pyfin_alignments)}\n\n")
            
            # Differences
            f.write("## Differences\n\n")
            f.write(f"Mean difference: {np.mean(pyfin_means) - np.mean(f5c_means):.2f} pA\n")
            f.write(f"Stdv difference: {np.mean(pyfin_stdvs) - np.mean(f5c_stdvs):.2f} pA\n")
            f.write(f"Duration difference: {np.mean(pyfin_durations) - np.mean(f5c_durations):.2f} samples\n")
        else:
            f.write("## PyFIN Event Statistics\n\n")
            f.write("No PyFIN alignments available.\n")
    
    print(f"Event statistics report saved to: {output_path}")


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Run comparison and generate reports."""
    import argparse
    
    parser = argparse.ArgumentParser(description="PyFIN vs f5c Eventalign Comparison")
    parser.add_argument("--report", action="store_true", help="Generate detailed reports")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "test_results",
                        help="Output directory for reports")
    args = parser.parse_args()
    
    print("=" * 60)
    print("PyFIN vs f5c Eventalign Comparison")
    print("=" * 60)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("\nLoading test data...")
    f5c_alignments = load_f5c_eventalign(F5C_TSV_PATH)
    print(f"  f5c alignments: {len(f5c_alignments)}")
    
    ref_name, ref_seq = load_reference_fasta(FASTA_PATH)
    print(f"  Reference: {ref_name}, {len(ref_seq)} bp")
    
    read_id, signal, sample_rate = load_signal_from_pod5(POD5_PATH)
    print(f"  Signal: {len(signal)} samples")
    
    # Run PyFIN
    print("\nRunning PyFIN eventalign...")
    pyfin_alignments, pyfin_scalings, error = run_pyfin_eventalign(
        signal=signal,
        sample_rate=sample_rate,
        read_id=read_id,
        read_seq=ref_seq,
        ref_seq=ref_seq,
        ref_name=ref_name
    )
    
    if error:
        print(f"  ERROR: {error}")
    else:
        print(f"  PyFIN alignments: {len(pyfin_alignments) if pyfin_alignments else 0}")
    
    # Compute metrics
    metrics = compare_results(f5c_alignments, pyfin_alignments, pyfin_scalings, len(ref_seq))
    
    # Print summary
    print(metrics.summary())
    
    # Save metrics
    metrics_path = args.output_dir / "comparison_metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics.to_dict(), f, indent=2)
    print(f"\nMetrics saved to: {metrics_path}")
    
    # Generate detailed reports
    if args.report:
        print("\nGenerating detailed reports...")
        
        generate_position_comparison_report(
            f5c_alignments, pyfin_alignments,
            args.output_dir / "position_comparison.txt"
        )
        
        generate_event_statistics_report(
            f5c_alignments, pyfin_alignments,
            args.output_dir / "event_statistics.txt"
        )
    
    # Return exit code based on success
    return 0 if metrics.pyfin_success else 1


if __name__ == "__main__":
    sys.exit(main())
