"""
Sequence utility functions.

This module provides basic sequence manipulation utilities.
"""

# Simple reverse complement implementation to avoid circular imports
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
              'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
              'N': 'N', 'n': 'n'}

def reverse_complement(seq: str) -> str:
    """
    Compute reverse complement of DNA/RNA sequence.

    Args:
        seq: DNA/RNA sequence string

    Returns:
        Reverse complement string
    """
    # Use .get() to handle unknown characters gracefully
    return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))


__all__ = [
    "reverse_complement",
]
