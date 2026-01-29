#!/usr/bin/env python3
"""
Unit tests for fin.utils.sequences module.

Tests the sequence utility functions:
- reverse_complement: DNA/RNA reverse complement calculation

Run with:
    pytest tests/unit/test_sequences.py -v
"""

import pytest
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.utils.sequences import reverse_complement, COMPLEMENT


class TestComplement:
    """Tests for the COMPLEMENT dictionary."""

    def test_complement_contains_standard_bases(self):
        """Test that COMPLEMENT contains all standard DNA bases."""
        assert 'A' in COMPLEMENT
        assert 'T' in COMPLEMENT
        assert 'C' in COMPLEMENT
        assert 'G' in COMPLEMENT
        assert 'N' in COMPLEMENT

    def test_complement_contains_lowercase_bases(self):
        """Test that COMPLEMENT contains lowercase DNA bases."""
        assert 'a' in COMPLEMENT
        assert 't' in COMPLEMENT
        assert 'c' in COMPLEMENT
        assert 'g' in COMPLEMENT
        assert 'n' in COMPLEMENT

    def test_complement_pairs_uppercase(self):
        """Test that uppercase complement pairs are correct."""
        assert COMPLEMENT['A'] == 'T'
        assert COMPLEMENT['T'] == 'A'
        assert COMPLEMENT['C'] == 'G'
        assert COMPLEMENT['G'] == 'C'
        assert COMPLEMENT['N'] == 'N'

    def test_complement_pairs_lowercase(self):
        """Test that lowercase complement pairs are correct."""
        assert COMPLEMENT['a'] == 't'
        assert COMPLEMENT['t'] == 'a'
        assert COMPLEMENT['c'] == 'g'
        assert COMPLEMENT['g'] == 'c'
        assert COMPLEMENT['n'] == 'n'


class TestReverseComplement:
    """Tests for the reverse_complement function."""

    # Basic functionality tests
    def test_reverse_complement_empty_string(self):
        """Test reverse complement of empty string returns empty string."""
        assert reverse_complement("") == ""

    def test_reverse_complement_single_base_A(self):
        """Test reverse complement of single 'A' returns 'T'."""
        assert reverse_complement("A") == "T"

    def test_reverse_complement_single_base_T(self):
        """Test reverse complement of single 'T' returns 'A'."""
        assert reverse_complement("T") == "A"

    def test_reverse_complement_single_base_C(self):
        """Test reverse complement of single 'C' returns 'G'."""
        assert reverse_complement("C") == "G"

    def test_reverse_complement_single_base_G(self):
        """Test reverse complement of single 'G' returns 'C'."""
        assert reverse_complement("G") == "C"

    def test_reverse_complement_single_base_N(self):
        """Test reverse complement of single 'N' returns 'N'."""
        assert reverse_complement("N") == "N"

    # Standard DNA sequences
    def test_reverse_complement_simple_sequence(self):
        """Test reverse complement of a simple sequence."""
        # ATCG -> complement: TAGC -> reverse: CGAT
        assert reverse_complement("ATCG") == "CGAT"

    def test_reverse_complement_all_same_base(self):
        """Test reverse complement of sequence with all same base."""
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("TTTT") == "AAAA"
        assert reverse_complement("CCCC") == "GGGG"
        assert reverse_complement("GGGG") == "CCCC"

    def test_reverse_complement_palindrome(self):
        """Test reverse complement of a palindromic sequence."""
        # GCGC -> complement: CGCG -> reverse: GCGC
        assert reverse_complement("GCGC") == "GCGC"

    def test_reverse_complement_longer_sequence(self):
        """Test reverse complement of a longer sequence."""
        seq = "ACGTACGTACGT"
        expected = "ACGTACGTACGT"
        assert reverse_complement(seq) == expected

    # Case handling
    def test_reverse_complement_lowercase(self):
        """Test reverse complement preserves lowercase."""
        assert reverse_complement("atcg") == "cgat"

    def test_reverse_complement_mixed_case(self):
        """Test reverse complement with mixed case sequence."""
        assert reverse_complement("AtCg") == "cGaT"

    # Sequences with N
    def test_reverse_complement_with_N(self):
        """Test reverse complement handles N bases."""
        assert reverse_complement("ACNGT") == "ACNGT"
        assert reverse_complement("NNNNN") == "NNNNN"

    # Inverse property (running twice should return original)
    def test_reverse_complement_is_involution(self):
        """Test that applying reverse_complement twice returns original."""
        sequences = ["ATCG", "ACGTACGT", "NNNNN", "AtCgN", ""]
        for seq in sequences:
            assert reverse_complement(reverse_complement(seq)) == seq

    # Edge cases
    def test_reverse_complement_unknown_character(self):
        """Test that unknown characters are passed through unchanged."""
        # The function uses .get(base, base) so unknown chars are kept
        assert reverse_complement("ACXGT") == "ACXGT"

    def test_reverse_complement_with_numbers(self):
        """Test behavior with non-DNA characters (numbers)."""
        # Numbers should be passed through unchanged but the sequence is reversed
        # "A1T" -> complement: "T1A" -> reverse: "A1T" (since 1 stays as 1)
        assert reverse_complement("A1T") == "A1T"

    # Real biological sequence test
    def test_reverse_complement_biological_primer(self):
        """Test reverse complement of a realistic primer sequence."""
        # Forward primer
        forward = "ATGGCGATCGATCGATCG"
        # Expected reverse complement
        expected = "CGATCGATCGATCGCCAT"
        assert reverse_complement(forward) == expected


class TestReverseComplementPerformance:
    """Performance-related tests for reverse_complement."""

    def test_reverse_complement_long_sequence(self):
        """Test reverse complement handles long sequences."""
        # Create a long sequence
        long_seq = "ATCG" * 10000  # 40,000 bp
        result = reverse_complement(long_seq)
        assert len(result) == len(long_seq)
        # Verify first and last bases are correct
        assert result[-1] == "T"  # complement of first 'A'
        assert result[0] == "C"   # complement of last 'G'


# Pytest fixtures for parameterized testing
@pytest.mark.parametrize("input_seq,expected", [
    ("", ""),
    ("A", "T"),
    ("T", "A"),
    ("C", "G"),
    ("G", "C"),
    ("ATCG", "CGAT"),
    ("atcg", "cgat"),
    ("AAAA", "TTTT"),
    ("GCGC", "GCGC"),  # palindrome
])
def test_reverse_complement_parametrized(input_seq, expected):
    """Parameterized tests for reverse_complement function."""
    assert reverse_complement(input_seq) == expected


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
