#!/usr/bin/env python3
"""
Eventalign Comparison Test: PyFIN run_eventalign vs f5c

This test compares the full eventalign pipeline output between PyFIN and f5c.
It loads the necessary inputs from POD5, BAM, and FASTA files, runs PyFIN's
eventalign, and compares with pre-computed f5c results.

Required inputs for PyFIN run_eventalign:
- read_ids: list of read identifier strings
- read_seqs: list of basecalled read sequences
- ref_seqs: list of reference sequences
- ref_names: list of reference names
- ref_lens: list of reference lengths
- signals: list of raw signal arrays (float32)
- sample_rates: list of sample rates (Hz)
- model_id: MODEL_RNA002 or MODEL_RNA004

Usage:
    python tests/run_eventalign_comparison.py
"""
import sys
sys.path.insert(0, "/home/logan/Projects/pyfin")

import gzip
from pathlib import Path
from collections import defaultdict
import numpy as np

print("=" * 70)
print("EVENTALIGN COMPARISON: PyFIN run_eventalign vs f5c")
print("=" * 70)

# Test data paths
TEST_DATA_DIR = Path("/home/logan/Projects/pyfin/tests/testdata")
RNA004_POD5_PATH = TEST_DATA_DIR / "RNA004.test.pod5"
RNA004_BAM_PATH = TEST_DATA_DIR / "RNA004.test.bam"
RNA004_F5C_TSV_PATH = TEST_DATA_DIR / "RNA004.test.tsv.gz"
REFERENCE_PATH = TEST_DATA_DIR / "test.fa"

# ============================================================================
# Step 1: Load f5c eventalign results (ground truth)
# ============================================================================
print("\n[1] Loading f5c eventalign results...")

f5c_events = defaultdict(list)  # {read_id: [events]}
with gzip.open(RNA004_F5C_TSV_PATH, "rt") as f:
    header = f.readline()  # Skip header
    
    count = 0
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 16:
            continue
        
        read_id = parts[3]
        
        # Limit reads for testing
        if read_id not in f5c_events and len(f5c_events) >= 5:
            continue
        
        event = {
            "contig": parts[0],
            "position": int(parts[1]),
            "reference_kmer": parts[2],
            "read_name": parts[3],
            "strand": parts[4],
            "event_index": int(parts[5]),
            "event_level_mean": float(parts[6]),
            "event_stdv": float(parts[7]),
            "event_length": float(parts[8]),
            "model_kmer": parts[9],
            "model_mean": float(parts[10]),
            "model_stdv": float(parts[11]),
            "standardized_level": float(parts[12]) if parts[12] != "inf" else 0.0,
            "start_idx": int(parts[13]),
            "end_idx": int(parts[14]),
        }
        f5c_events[read_id].append(event)
        count += 1

print(f"  Loaded {len(f5c_events)} reads with {count:,} events from f5c")

# Get the reference names used by each read
read_to_ref = {}
for read_id, events in f5c_events.items():
    refs = set(e["contig"] for e in events)
    read_to_ref[read_id] = list(refs)[0]  # Assume single reference per read
    print(f"    {read_id[:12]}...: {len(events):,} events -> {read_to_ref[read_id]}")

# ============================================================================
# Step 2: Load reference sequences
# ============================================================================
print("\n[2] Loading reference sequences...")

def load_fasta(fasta_path):
    """Load FASTA file into dict."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:
            sequences[current_name] = "".join(current_seq)
    
    return sequences

reference_seqs = load_fasta(REFERENCE_PATH)
print(f"  Loaded {len(reference_seqs)} reference sequences")

# ============================================================================
# Step 3: Load read sequences from BAM
# ============================================================================
print("\n[3] Loading read sequences from BAM...")

import pysam

read_seqs = {}
bam = pysam.AlignmentFile(str(RNA004_BAM_PATH), "rb")
for read in bam:
    if read.query_name in f5c_events:
        read_seqs[read.query_name] = read.query_sequence
bam.close()

print(f"  Loaded {len(read_seqs)} read sequences")

# ============================================================================
# Step 4: Load raw signals from POD5
# ============================================================================
print("\n[4] Loading raw signals from POD5...")

import pod5

signals = {}
sample_rates = {}

with pod5.Reader(str(RNA004_POD5_PATH)) as reader:
    for read in reader.reads():
        rid = str(read.read_id)
        if rid in f5c_events:
            signals[rid] = read.signal_pa.astype(np.float32)
            sample_rates[rid] = float(read.run_info.sample_rate)

print(f"  Loaded {len(signals)} signals")
for rid in list(signals.keys())[:3]:
    print(f"    {rid[:12]}...: {len(signals[rid]):,} samples, {sample_rates[rid]} Hz")

# ============================================================================
# Step 5: Prepare inputs for PyFIN run_eventalign
# ============================================================================
print("\n[5] Preparing inputs for PyFIN eventalign...")

# Find common reads
common_reads = set(f5c_events.keys()) & set(read_seqs.keys()) & set(signals.keys())
print(f"  Common reads: {len(common_reads)}")

if not common_reads:
    print("ERROR: No common reads found!")
    sys.exit(1)

# Prepare input lists
read_ids_list = []
read_seqs_list = []
ref_seqs_list = []
ref_names_list = []
ref_lens_list = []
signals_list = []
sample_rates_list = []

for read_id in common_reads:
    ref_name = read_to_ref[read_id]
    
    if ref_name not in reference_seqs:
        print(f"  WARNING: Reference {ref_name} not found, skipping {read_id}")
        continue
    
    ref_seq = reference_seqs[ref_name]
    
    read_ids_list.append(read_id)
    read_seqs_list.append(read_seqs[read_id])
    ref_seqs_list.append(ref_seq)
    ref_names_list.append(ref_name)
    ref_lens_list.append(len(ref_seq))
    signals_list.append(signals[read_id])
    sample_rates_list.append(sample_rates[read_id])

print(f"  Prepared {len(read_ids_list)} reads for eventalign")

# ============================================================================
# Step 6: Run PyFIN eventalign
# ============================================================================
print("\n[6] Running PyFIN eventalign...")

from fin._eventalign import run_eventalign, MODEL_RNA004

try:
    result = run_eventalign(
        read_ids=read_ids_list,
        read_seqs=read_seqs_list,
        ref_seqs=ref_seqs_list,
        ref_names=ref_names_list,
        ref_lens=ref_lens_list,
        signals=signals_list,
        sample_rates=sample_rates_list,
        model_id=MODEL_RNA004,  # Use RNA004 model for RNA004 data
    )
    print("  ✓ PyFIN eventalign completed successfully")
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# Step 7: Compare results
# ============================================================================
print("\n[7] Comparing PyFIN vs f5c results...")
print("=" * 70)

# Extract result components
pyfin_full = result["full"]  # [read_idx][ref_idx] = list of alignment dicts
pyfin_scalings = result["scalings"]  # [read_idx] = {scale, shift, var}
pyfin_events = result["events"]  # [read_idx] = {starts, lengths, means, stdvs}
summary = result["summary"]

print(f"\nPyFIN Summary:")
print(f"  Reads processed: {summary['num_reads']}")
print(f"  References: {summary['num_refs']}")

for read_idx, read_id in enumerate(read_ids_list):
    print(f"\n{'='*60}")
    print(f"Read: {read_id}")
    print(f"{'='*60}")
    
    # f5c events for this read
    f5c_read_events = f5c_events[read_id]
    
    # PyFIN events
    pyfin_ev = pyfin_events[read_idx]
    pyfin_align = pyfin_full[read_idx][0]  # First (only) reference
    pyfin_scale = pyfin_scalings[read_idx]
    
    print(f"\n  Event Detection:")
    print(f"    f5c events (in alignment): {len(f5c_read_events):,}")
    print(f"    PyFIN events (detected):   {pyfin_ev['n_events']:,}")
    
    print(f"\n  Scaling Parameters:")
    print(f"    PyFIN: scale={pyfin_scale['scale']:.4f}, shift={pyfin_scale['shift']:.4f}, var={pyfin_scale['var']:.4f}")
    
    print(f"\n  Alignment Results:")
    print(f"    PyFIN alignments: {len(pyfin_align):,}")
    
    if pyfin_align:
        # Compare alignment positions
        f5c_positions = set((e["position"], e["event_index"]) for e in f5c_read_events)
        pyfin_positions = set((a["ref_position"], a["event_idx"]) for a in pyfin_align)
        
        common = f5c_positions & pyfin_positions
        f5c_only = f5c_positions - pyfin_positions
        pyfin_only = pyfin_positions - f5c_positions
        
        print(f"    Common (pos, event_idx): {len(common):,}")
        print(f"    f5c only:                {len(f5c_only):,}")
        print(f"    PyFIN only:              {len(pyfin_only):,}")
        
        if f5c_positions:
            match_rate = len(common) / len(f5c_positions) * 100
            print(f"    Match rate:              {match_rate:.1f}%")
        
        # Compare k-mers for matching positions
        if common:
            f5c_kmer_dict = {(e["position"], e["event_index"]): e["reference_kmer"] for e in f5c_read_events}
            pyfin_kmer_dict = {(a["ref_position"], a["event_idx"]): a["ref_kmer"] for a in pyfin_align}
            
            kmer_matches = 0
            for pos_ev in list(common)[:100]:  # Check first 100
                if f5c_kmer_dict[pos_ev] == pyfin_kmer_dict[pos_ev]:
                    kmer_matches += 1
            
            print(f"    K-mer match (sample):    {kmer_matches}/100")
    
    # Show first few alignments from each
    print(f"\n  Sample Alignments (first 5):")
    print(f"    f5c:   {[(e['position'], e['reference_kmer'], e['event_index']) for e in f5c_read_events[:5]]}")
    if pyfin_align:
        print(f"    PyFIN: {[(a['ref_position'], a['ref_kmer'], a['event_idx']) for a in pyfin_align[:5]]}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

total_f5c_align = sum(len(f5c_events[rid]) for rid in read_ids_list)
total_pyfin_align = sum(len(pyfin_full[i][0]) for i in range(len(read_ids_list)))

print(f"\n  Total f5c alignments:   {total_f5c_align:,}")
print(f"  Total PyFIN alignments: {total_pyfin_align:,}")
print(f"  Ratio:                  {total_pyfin_align/total_f5c_align:.3f}")

print("\n" + "=" * 70)
print("TEST COMPLETE")
print("=" * 70)
