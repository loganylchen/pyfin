#!/usr/bin/env python3
"""
Quick test to verify the test data files load correctly.
Run this to check the test infrastructure before running full tests.
"""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Test data paths
TEST_DATA_DIR = PROJECT_ROOT / "tests" / "testdata"


def test_paths_exist():
    """Verify test data files exist."""
    print("\n=== Test Data Files ===")
    
    files = {
        "RNA004 POD5": TEST_DATA_DIR / "RNA004.test.pod5",
        "RNA004 FASTQ": TEST_DATA_DIR / "RNA004.test.fq.gz",
        "RNA004 BAM": TEST_DATA_DIR / "RNA004.test.bam",
        "RNA004 f5c TSV": TEST_DATA_DIR / "RNA004.test.tsv.gz",
        "Reference FASTA": TEST_DATA_DIR / "test.fa",
        "BLOW5": TEST_DATA_DIR / "RNA004.test.blow5",
    }
    
    all_exist = True
    for name, path in files.items():
        exists = path.exists()
        size = path.stat().st_size / 1024 / 1024 if exists else 0
        status = f"✓ {size:.1f} MB" if exists else "✗ MISSING"
        print(f"  {name}: {status}")
        all_exist = all_exist and exists
    
    return all_exist


def test_f5c_tsv_format():
    """Verify f5c TSV file can be parsed."""
    import gzip
    
    tsv_path = TEST_DATA_DIR / "RNA004.test.tsv.gz"
    if not tsv_path.exists():
        print("✗ f5c TSV file not found")
        return False
    
    print("\n=== f5c TSV Format Check ===")
    
    with gzip.open(tsv_path, 'rt') as f:
        # Read header
        header = f.readline().strip().split('\t')
        print(f"  Columns: {len(header)}")
        print(f"  Header: {header[:8]}...")
        
        # Read first data line
        first_line = f.readline().strip().split('\t')
        print(f"  Sample read_id: {first_line[3][:20]}...")
        print(f"  Sample reference: {first_line[0]}")
        
        # Sample first 10000 lines only (file is too large to fully count)
        read_ids = {first_line[3]}
        refs = {first_line[0]}
        sample_count = 2  # header + first line already read
        
        for i, line in enumerate(f):
            if i >= 10000:
                break
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                read_ids.add(parts[3])
                refs.add(parts[0])
            sample_count += 1
    
    print(f"\n  Sampled lines: {sample_count:,} (first 10K only)")
    print(f"  Unique reads (in sample): {len(read_ids)}")
    print(f"  Unique refs (in sample): {len(refs)}")
    print(f"  References: {sorted(refs)[:5]}...")
    
    return True


def test_pod5_loads():
    """Verify POD5 file can be loaded."""
    try:
        import pod5
    except ImportError:
        print("\n=== POD5 Check ===")
        print("  ✗ pod5 not installed")
        return False
    
    pod5_path = TEST_DATA_DIR / "RNA004.test.pod5"
    if not pod5_path.exists():
        print("  ✗ POD5 file not found")
        return False
    
    print("\n=== POD5 Check ===")
    
    read_count = 0
    with pod5.Reader(str(pod5_path)) as reader:
        for read in reader.reads():
            read_count += 1
            if read_count == 1:
                signal = read.signal_pa
                sample_rate = read.run_info.sample_rate
                print(f"  First read: {read.read_id}")
                print(f"  Signal length: {len(signal):,} samples")
                print(f"  Sample rate: {sample_rate} Hz")
                print(f"  Duration: {len(signal)/sample_rate:.2f} s")
    
    print(f"  Total reads: {read_count}")
    return True


def test_bam_loads():
    """Verify BAM file can be loaded."""
    try:
        import pysam
    except ImportError:
        print("\n=== BAM Check ===")
        print("  ✗ pysam not installed")
        return False
    
    bam_path = TEST_DATA_DIR / "RNA004.test.bam"
    if not bam_path.exists():
        print("  ✗ BAM file not found")
        return False
    
    print("\n=== BAM Check ===")
    
    read_count = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            read_count += 1
            if read_count == 1:
                print(f"  First read: {read.query_name}")
                print(f"  Reference: {read.reference_name}")
                print(f"  Position: {read.reference_start}-{read.reference_end}")
    
    print(f"  Total alignments: {read_count}")
    return True


def test_reference_loads():
    """Verify reference FASTA can be loaded."""
    fasta_path = TEST_DATA_DIR / "test.fa"
    if not fasta_path.exists():
        print("\n=== Reference FASTA Check ===")
        print("  ✗ FASTA file not found")
        return False
    
    print("\n=== Reference FASTA Check ===")
    
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    print(f"  Total sequences: {len(sequences)}")
    print(f"  Sequence names: {list(sequences.keys())[:5]}...")
    total_bp = sum(len(s) for s in sequences.values())
    print(f"  Total bp: {total_bp:,}")
    
    return True


if __name__ == "__main__":
    print("=" * 60)
    print("PyFIN Test Data Verification")
    print("=" * 60)
    
    results = []
    results.append(("File paths", test_paths_exist()))
    results.append(("f5c TSV format", test_f5c_tsv_format()))
    results.append(("POD5 loading", test_pod5_loads()))
    results.append(("BAM loading", test_bam_loads()))
    results.append(("Reference FASTA", test_reference_loads()))
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {name}: {status}")
    
    all_passed = all(r[1] for r in results)
    print("\n" + ("All tests passed!" if all_passed else "Some tests failed."))
    sys.exit(0 if all_passed else 1)
