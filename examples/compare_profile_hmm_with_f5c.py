#!/usr/bin/env python3
"""
Compare PyFin Profile HMM eventalign with f5c eventalign output

Usage:
    # First run f5c:
    f5c eventalign -b reads.bam -g ref.fa -r reads.fastq > f5c_output.tsv

    # Then compare:
    python compare_profile_hmm_with_f5c.py f5c_output.tsv reads.fastq ref.fa

This will:
1. Parse f5c eventalign TSV output
2. Run PyFin profile_hmm_eventalign on the same data
3. Compare the alignments
"""

import sys
import numpy as np
from pathlib import Path

# Add parent directory to path for development
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from fin._f5c import profile_hmm_eventalign


def parse_f5c_eventalign(tsv_file):
    """Parse f5c eventalign TSV output"""
    import csv

    alignments = []
    with open(tsv_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            alignments.append(
                {
                    "contig": row["contig"],
                    "position": int(row["position"]),
                    "reference_kmer": row["reference_kmer"],
                    "read_index": int(row["read_index"]),
                    "strand": row["strand"],
                    "event_index": int(row["event_index"]),
                    "event_level_mean": float(row["event_level_mean"]),
                    "event_stdv": float(row["event_stdv"]),
                    "event_length": float(row["event_length"]),
                    "model_kmer": row["model_kmer"],
                    "model_mean": float(row["model_mean"]),
                    "model_stdv": float(row["model_stdv"]),
                    "standardized_level": float(row["standardized_level"]),
                    "start_idx": int(row.get("start_idx", 0)),
                    "end_idx": int(row.get("end_idx", 0)),
                }
            )

    return alignments


def compare_alignments(f5c_aln, pyfin_aln):
    """Compare f5c and PyFin alignments"""

    print("=" * 70)
    print("Alignment Comparison: f5c vs PyFin")
    print("=" * 70)

    print(f"\nf5c alignment records: {len(f5c_aln)}")
    print(f"PyFin alignment records: {len(pyfin_aln)}")

    # Compare first N records
    n_compare = min(20, len(f5c_aln), len(pyfin_aln))

    print(f"\nFirst {n_compare} records comparison:")
    print("-" * 70)
    print(
        f"{'Idx':<4} {'f5c Pos':<8} {'PyFin Pos':<10} {'f5c Kmer':<10} {'PyFin Kmer':<11} {'Match':<5}"
    )
    print("-" * 70)

    matches = 0
    for i in range(n_compare):
        f5c = f5c_aln[i]
        pyfin = pyfin_aln[i]

        pos_match = f5c["position"] == pyfin["ref_position"]
        kmer_match = f5c["reference_kmer"] == pyfin["ref_kmer"]
        match = pos_match and kmer_match

        if match:
            matches += 1

        print(
            f"{i:<4} {f5c['position']:<8} {pyfin['ref_position']:<10} "
            f"{f5c['reference_kmer']:<10} {pyfin['ref_kmer']:<11} "
            f"{'✓' if match else '✗':<5}"
        )

    print("-" * 70)
    print(f"Matches: {matches}/{n_compare} ({100*matches/n_compare:.1f}%)")

    # Compare event indices
    print(f"\nEvent index comparison:")
    f5c_events = [a["event_index"] for a in f5c_aln[:n_compare]]
    pyfin_events = [a["event_idx"] for a in pyfin_aln[:n_compare]]

    event_matches = sum(1 for f, p in zip(f5c_events, pyfin_events) if f == p)
    print(
        f"  Event index matches: {event_matches}/{n_compare} ({100*event_matches/n_compare:.1f}%)"
    )

    # Compare model means
    print(f"\nModel mean comparison (first 5):")
    for i in range(min(5, n_compare)):
        f5c = f5c_aln[i]
        pyfin = pyfin_aln[i]
        diff = abs(f5c["model_mean"] - pyfin["model_mean"])
        print(
            f"  {i}: f5c={f5c['model_mean']:.3f}, pyfin={pyfin['model_mean']:.3f}, diff={diff:.6f}"
        )


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    f5c_tsv = sys.argv[1]
    fastq_file = sys.argv[2]
    ref_fasta = sys.argv[3]

    # Parse f5c output
    print("Loading f5c eventalign output...")
    f5c_alignments = parse_f5c_eventalign(f5c_tsv)
    print(f"Loaded {len(f5c_alignments)} f5c alignment records")

    # TODO: Load reads and reference from files
    # For now, just show how to use profile_hmm_eventalign
    print("\nTo run PyFin profile_hmm_eventalign:")
    print("```python")
    print("from fin import _eventalign")
    print("import numpy as np")
    print()
    print("# Load your raw signal as float32 numpy array")
    print("signal = np.array([...], dtype=np.float32)")
    print()
    print("# Load reference sequence")
    print("sequence = 'ACGT...'")
    print()
    print("# Run Profile HMM eventalign")
    print("result = _eventalign.profile_hmm_eventalign(")
    print("    raw_signal=signal,")
    print("    sequence=sequence,")
    print("    is_rna=1,  # or 0 for DNA")
    print("    kmer_size=5,  # or 9 for RNA004")
    print("    events_per_base=3.0")
    print(")")
    print()
    print("# Access alignment")
    print("for aln in result['alignment']:")
    print(
        "    print(f\"{aln['ref_position']}: {aln['ref_kmer']} -> event {aln['event_idx']} (state={aln['hmm_state']})\")"
    )
    print("```")


if __name__ == "__main__":
    main()
