#!/usr/bin/env python
"""
Debug script to isolate run_eventalign segfault.
"""
import sys
import gzip
import numpy as np
from pathlib import Path

# Paths - use the existing test data
TEST_DATA_DIR = Path("/home/logan/Projects/pyfin/tests/testdata")
RNA004_POD5_PATH = TEST_DATA_DIR / "RNA004.test.pod5"
RNA004_BAM_PATH = TEST_DATA_DIR / "RNA004.test.bam"
REFERENCE_PATH = TEST_DATA_DIR / "test.fa"
F5C_TSV = TEST_DATA_DIR / "RNA004.test.tsv.gz"

print("=" * 60)
print("DEBUG: run_eventalign segfault investigation")
print("=" * 60)

# Step 1: Load ONE read from f5c to identify a test case
print("\n[1] Loading ONE read from f5c...")
f5c_read_id = None
f5c_ref_name = None

with gzip.open(F5C_TSV, "rt") as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            f5c_read_id = parts[3]
            f5c_ref_name = parts[0]
            break

print(f"  Test read: {f5c_read_id}")
print(f"  Test reference: {f5c_ref_name}")

# Step 2: Load reference sequence
print("\n[2] Loading reference sequence...")
def load_fasta(path):
    seqs = {}
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs

refs = load_fasta(REFERENCE_PATH)
ref_seq = refs.get(f5c_ref_name, "")
print(f"  Reference {f5c_ref_name}: {len(ref_seq)} bases")

if not ref_seq:
    print(f"  ERROR: Reference {f5c_ref_name} not found in FASTA!")
    sys.exit(1)

# Step 3: Load read sequence from BAM
print("\n[3] Loading read sequence from BAM...")
import pysam

read_seq = None
bam = pysam.AlignmentFile(str(RNA004_BAM_PATH), "rb")
for read in bam:
    if read.query_name == f5c_read_id:
        read_seq = read.query_sequence
        break
bam.close()

print(f"  Read sequence: {len(read_seq) if read_seq else 0} bases")

if not read_seq:
    print(f"  ERROR: Read {f5c_read_id} not found in BAM!")
    sys.exit(1)

# Step 4: Load signal from POD5
print("\n[4] Loading signal from POD5...")
import pod5

signal = None
sample_rate = None

with pod5.Reader(str(RNA004_POD5_PATH)) as reader:
    for read in reader.reads():
        if str(read.read_id) == f5c_read_id:
            signal = read.signal_pa.astype(np.float32)
            sample_rate = float(read.run_info.sample_rate)
            break

print(f"  Signal: {len(signal) if signal is not None else 0} samples")
print(f"  Sample rate: {sample_rate} Hz")

if signal is None:
    print(f"  ERROR: Read {f5c_read_id} not found in POD5!")
    sys.exit(1)

# Step 5: Test run_eventalign with explicit type checks
print("\n[5] Preparing inputs with type checks...")

# Build the lists
read_ids = [f5c_read_id]
read_seqs = [read_seq]
ref_seqs = [ref_seq]
ref_names = [f5c_ref_name]
ref_lens = [len(ref_seq)]
signals = [signal]
sample_rates = [sample_rate]

print(f"  read_ids type: {type(read_ids)}, element type: {type(read_ids[0])}")
print(f"  read_seqs type: {type(read_seqs)}, element type: {type(read_seqs[0])}")
print(f"  ref_seqs type: {type(ref_seqs)}, element type: {type(ref_seqs[0])}")
print(f"  ref_names type: {type(ref_names)}, element type: {type(ref_names[0])}")
print(f"  ref_lens type: {type(ref_lens)}, element type: {type(ref_lens[0])}")
print(f"  signals type: {type(signals)}, element type: {type(signals[0])}, dtype: {signals[0].dtype}")
print(f"  sample_rates type: {type(sample_rates)}, element type: {type(sample_rates[0])}")

# Check for any None values
print("\n  Checking for None/empty values:")
print(f"    read_id: {read_ids[0] is not None}, len={len(read_ids[0])}")
print(f"    read_seq: {read_seqs[0] is not None}, len={len(read_seqs[0])}")
print(f"    ref_seq: {ref_seqs[0] is not None}, len={len(ref_seqs[0])}")
print(f"    ref_name: {ref_names[0] is not None}, len={len(ref_names[0])}")
print(f"    ref_len: {ref_lens[0]} (expected: {len(ref_seqs[0])})")
print(f"    signal: {signals[0] is not None}, len={len(signals[0])}")
print(f"    sample_rate: {sample_rates[0]}")

# Step 6: Try run_eventalign
print("\n[6] Testing run_eventalign...")
sys.stdout.flush()

from fin._eventalign import run_eventalign, MODEL_RNA004

print("  Importing successful, calling run_eventalign...")
sys.stdout.flush()

try:
    result = run_eventalign(
        read_ids=read_ids,
        read_seqs=read_seqs,
        ref_seqs=ref_seqs,
        ref_names=ref_names,
        ref_lens=ref_lens,
        signals=signals,
        sample_rates=sample_rates,
        model_id=MODEL_RNA004,
    )
    
    print("  ✓ run_eventalign completed!")
    print(f"\n  Result keys: {list(result.keys())}")
    print(f"  Summary: {result.get('summary', {})}")
    
    # Check events
    if 'events' in result and len(result['events']) > 0:
        ev = result['events'][0]
        print(f"\n  Events detected: {ev.get('n_events', len(ev.get('means', [])))} events")
    
    # Check alignment
    if 'full' in result and len(result['full']) > 0:
        align = result['full'][0][0]
        print(f"  Alignment entries: {len(align)}")
        if len(align) > 0:
            print(f"  First entry: {align[0]}")
    
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 60)
print("DEBUG COMPLETE - NO SEGFAULT!")
print("=" * 60)
