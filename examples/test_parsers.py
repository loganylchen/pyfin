#!/usr/bin/env python3
"""
Example usage of format parsers in the fin package
"""

import sys
from pathlib import Path

# Add parent directory to path for importing fin
sys.path.insert(0, str(Path(__file__).parent.parent))

from fin.io import Fast5Reader, Slow5Reader, Pod5Reader, BamReader, GTFReader, BEDReader


def test_fast5(file_path):
    """Test Fast5Reader"""
    print(f"\n=== Testing Fast5Reader: {file_path} ===")
    try:
        with Fast5Reader(file_path) as reader:
            print(f"Read IDs: {len(reader.read_ids[:5])}... (showing first 5)")
            if reader.read_ids:
                first_read = reader.read_ids[0]
                signal, metadata = reader.get_raw_signal(first_read)
                print(f"Example signal length: {len(signal)}")
                print(f"Metadata: {metadata}")
    except Exception as e:
        print(f"Error: {e}")


def test_slow5(file_path):
    """Test Slow5Reader"""
    print(f"\n=== Testing Slow5Reader: {file_path} ===")
    try:
        with Slow5Reader(file_path) as reader:
            read_ids = reader.read_ids
            print(f"Read IDs: {len(read_ids[:5])}... (showing first 5)")
            if read_ids:
                first_read = read_ids[0]
                signal, metadata = reader.get_signal(first_read)
                print(f"Example signal length: {len(signal)}")
                print(f"Metadata keys: {list(metadata.keys())}")
    except Exception as e:
        print(f"Error: {e}")


def test_pod5(file_path):
    """Test Pod5Reader"""
    print(f"\n=== Testing Pod5Reader: {file_path} ===")
    try:
        with Pod5Reader(file_path) as reader:
            print(f"Total reads: {reader.read_count}")
            read_ids = reader.read_ids
            print(f"Read IDs: {len(read_ids[:5])}... (showing first 5)")
            if read_ids:
                first_read = read_ids[0]
                signal, metadata = reader.get_signal(first_read)
                print(f"Example signal length: {len(signal)}")
                print(f"Metadata keys: {list(metadata.keys())}")
    except Exception as e:
        print(f"Error: {e}")


def test_bam(file_path):
    """Test BamReader"""
    print(f"\n=== Testing BamReader: {file_path} ===")
    try:
        with BamReader(file_path) as reader:
            stats = reader.get_file_stats()
            print(f"File stats: {stats}")

            if 'nreferences' in stats and stats['nreferences'] > 0:
                print(f"References: {stats.get('references', [])[:3]}... (showing first 3)")
    except Exception as e:
        print(f"Error: {e}")


def test_gtf(file_path):
    """Test GTFReader"""
    print(f"\n=== Testing GTFReader: {file_path} ===")
    try:
        with GTFReader(file_path) as reader:
            reader.parse()
            stats = reader.get_gene_stats()
            print(f"Gene stats: {stats}")
    except Exception as e:
        print(f"Error: {e}")


def test_bed(file_path):
    """Test BEDReader"""
    print(f"\n=== Testing BEDReader: {file_path} ===")
    try:
        with BEDReader(file_path) as reader:
            features = list(reader.iterate_features())
            print(f"Total features: {len(features)}")
            if features:
                feature = features[0]
                print(f"First feature: {feature}")
                print(f"Feature length: {feature.length}")
    except Exception as e:
        print(f"Error: {e}")


def main():
    """Main test function"""
    print("fin Package Format Parsers Test")
    print("=" * 50)

    # Note: These are example paths - replace with actual files to test
    test_files = {
        'fast5': 'test_data/sample.fast5',
        'slow5': 'test_data/sample.blow5',
        'pod5': 'test_data/sample.pod5',
        'bam': 'test_data/sample.bam',
        'gtf': 'test_data/sample.gtf',
        'bed': 'test_data/sample.bed',
    }

    for format_type, file_path in test_files.items():
        if Path(file_path).exists():
            if format_type == 'fast5':
                test_fast5(file_path)
            elif format_type == 'slow5':
                test_slow5(file_path)
            elif format_type == 'pod5':
                test_pod5(file_path)
            elif format_type == 'bam':
                test_bam(file_path)
            elif format_type == 'gtf':
                test_gtf(file_path)
            elif format_type == 'bed':
                test_bed(file_path)
        else:
            print(f"\n=== Skipping {format_type}: File not found: {file_path} ===")
            print("Tip: Add test files to test_data/ directory or update paths in test_parsers.py")

    print("\n" + "=" * 50)
    print("Test completed!")


if __name__ == "__main__":
    main()
