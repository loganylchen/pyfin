"""Tests for intron chain extraction and read grouping."""

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.candidates.intron_chains import (
    extract_intron_chain,
    group_reads_by_three_prime_and_intron_chain,
    gtf_transcript_to_intron_chain,
    pick_representative_read,
)


class TestExtractIntronChain:
    """Tests for CIGAR -> intron chain extraction."""

    def test_single_exon_no_introns(self):
        """A simple match with no N ops should produce empty intron chain."""
        # 100M = 100bp match
        cigar = [(0, 100)]
        chain = extract_intron_chain(cigar, reference_start=1000)
        assert chain.introns == ()
        assert chain.is_single_exon
        assert chain.num_introns == 0

    def test_single_intron(self):
        """50M 500N 50M should produce one intron."""
        cigar = [(0, 50), (3, 500), (0, 50)]
        chain = extract_intron_chain(cigar, reference_start=1000)
        assert chain.introns == ((1050, 1550),)
        assert chain.num_introns == 1
        assert not chain.is_single_exon

    def test_two_introns(self):
        """30M 200N 40M 300N 30M should produce two introns."""
        cigar = [(0, 30), (3, 200), (0, 40), (3, 300), (0, 30)]
        chain = extract_intron_chain(cigar, reference_start=500)
        assert chain.introns == ((530, 730), (770, 1070))
        assert chain.num_introns == 2

    def test_soft_clips_ignored(self):
        """Soft clips (op=4) don't consume reference."""
        cigar = [(4, 10), (0, 50), (3, 100), (0, 50), (4, 5)]
        chain = extract_intron_chain(cigar, reference_start=0)
        assert chain.introns == ((50, 150),)

    def test_insertions_ignored(self):
        """Insertions (op=1) don't consume reference."""
        cigar = [(0, 30), (1, 5), (0, 20), (3, 100), (0, 50)]
        chain = extract_intron_chain(cigar, reference_start=0)
        # M consumes 30 ref, I doesn't, M consumes 20 ref -> intron at 50
        assert chain.introns == ((50, 150),)

    def test_deletions_consume_reference(self):
        """Deletions (op=2) consume reference."""
        cigar = [(0, 30), (2, 10), (0, 20), (3, 100), (0, 50)]
        chain = extract_intron_chain(cigar, reference_start=0)
        # M=30, D=10, M=20 -> ref_pos = 60 at N start
        assert chain.introns == ((60, 160),)

    def test_hashable_for_grouping(self):
        """IntronChain should be usable as dict key."""
        chain1 = extract_intron_chain([(0, 50), (3, 100), (0, 50)], 0)
        chain2 = extract_intron_chain([(0, 50), (3, 100), (0, 50)], 0)
        chain3 = extract_intron_chain([(0, 50), (3, 200), (0, 50)], 0)

        assert chain1 == chain2
        assert chain1 != chain3
        assert hash(chain1) == hash(chain2)

        d = {chain1: "a", chain3: "b"}
        assert d[chain2] == "a"


class TestGtfTranscriptToIntronChain:
    """Tests for GTF transcript -> intron chain conversion."""

    def test_single_exon_transcript(self):
        """Single exon transcript has no introns."""

        class MockTranscript:
            exons = [(100, 500)]

        chain = gtf_transcript_to_intron_chain(MockTranscript())
        assert chain.introns == ()

    def test_two_exon_transcript(self):
        """Two exons with gap = one intron."""

        class MockTranscript:
            exons = [(100, 300), (500, 700)]

        chain = gtf_transcript_to_intron_chain(MockTranscript())
        assert chain.introns == ((300, 500),)

    def test_three_exon_transcript(self):
        """Three exons = two introns."""

        class MockTranscript:
            exons = [(100, 200), (400, 500), (700, 800)]

        chain = gtf_transcript_to_intron_chain(MockTranscript())
        assert chain.introns == ((200, 400), (500, 700))

    def test_unsorted_exons(self):
        """Exons should be sorted before computing gaps."""

        class MockTranscript:
            exons = [(500, 700), (100, 300)]

        chain = gtf_transcript_to_intron_chain(MockTranscript())
        assert chain.introns == ((300, 500),)


class TestGroupReadsByThreePrimeAndIntronChain:
    """Tests for read grouping."""

    def _make_read(self, name, ref_start, cigar, ref_end, is_forward=True):
        return {
            "query_name": name,
            "reference_start": ref_start,
            "reference_end": ref_end,
            "cigartuples": cigar,
            "is_forward": is_forward,
            "query_length": ref_end - ref_start,
        }

    def test_same_chain_close_three_prime(self):
        """Reads with same intron chain and close 3' should group together."""
        reads = [
            self._make_read("r1", 100, [(0, 50), (3, 100), (0, 50)], 300),
            self._make_read("r2", 100, [(0, 50), (3, 100), (0, 50)], 305),
        ]
        groups = group_reads_by_three_prime_and_intron_chain(reads, "+", threshold=24)
        assert len(groups) == 1
        key = list(groups.keys())[0]
        assert len(groups[key]) == 2

    def test_different_chains_separate_groups(self):
        """Reads with different intron chains should be in separate groups."""
        reads = [
            self._make_read("r1", 100, [(0, 50), (3, 100), (0, 50)], 300),
            self._make_read("r2", 100, [(0, 50), (3, 200), (0, 50)], 300),
        ]
        groups = group_reads_by_three_prime_and_intron_chain(reads, "+", threshold=24)
        assert len(groups) == 2

    def test_distant_three_prime_separate_groups(self):
        """Reads with same chain but distant 3' should be in separate groups."""
        reads = [
            self._make_read("r1", 100, [(0, 100)], 200),
            self._make_read("r2", 100, [(0, 100)], 300),
        ]
        groups = group_reads_by_three_prime_and_intron_chain(reads, "+", threshold=24)
        assert len(groups) == 2


class TestPickRepresentativeRead:
    """Tests for representative read selection."""

    def test_picks_longest(self):
        reads = [
            {"query_name": "short", "query_length": 100},
            {"query_name": "long", "query_length": 500},
            {"query_name": "medium", "query_length": 250},
        ]
        rep = pick_representative_read(reads)
        assert rep["query_name"] == "long"
