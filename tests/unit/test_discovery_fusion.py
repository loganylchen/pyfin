"""Tests for merge_fusion_candidates (US-006 / AC-05)."""

from __future__ import annotations



from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.candidates.discovery import merge_fusion_candidates
from fin.io.interval_manager import GenomicInterval


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _interval() -> GenomicInterval:
    return GenomicInterval(chrom="chr1", start=1000, end=5000, strand="+")


def _chain() -> IntronChain:
    return IntronChain(introns=((1200, 1400), (2000, 2200)))


def _gtf_candidate(candidate_id: str, start: int = 1000, end: int = 5000) -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=candidate_id,
        intron_chain=_chain(),
        three_prime_pos=end,
        sequence="ACGT",
        source="gtf",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=start,
        end=end,
    )


def _fusion_candidate(candidate_id: str) -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=candidate_id,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=9000,
        sequence="TTTT",
        source="fusion",
        supporting_read_ids={"read1"},
        chrom="chr1",
        strand="+",
        start=8000,
        end=10000,
        fusion_junction=(8500, 8600),
    )


def _candidate_set(gtf_ids: list[str]) -> CandidateSet:
    return CandidateSet(
        interval=_interval(),
        candidates=[_gtf_candidate(cid) for cid in gtf_ids],
        read_ids={"r1", "r2"},
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_merge_basic():
    """3 GTF candidates + 1 fusion → result has 4 candidates with unique IDs."""
    cs = _candidate_set(["tx1", "tx2", "tx3"])
    fc = _fusion_candidate("fus1")

    result = merge_fusion_candidates(cs, [fc])

    assert result.num_candidates == 4
    ids = result.candidate_ids()
    assert len(ids) == len(set(ids)), "Duplicate IDs after merge"
    assert "fus1" in ids


def test_merge_preserves_original():
    """Original GTF candidates are unchanged (same objects, same order)."""
    gtf_cands = [_gtf_candidate("tx1"), _gtf_candidate("tx2"), _gtf_candidate("tx3")]
    cs = CandidateSet(interval=_interval(), candidates=list(gtf_cands), read_ids=set())

    merge_fusion_candidates(cs, [_fusion_candidate("fus1")])

    # First three entries are the original objects in original order
    assert cs.candidates[:3] == gtf_cands


def test_merge_deduplicates():
    """If fusion candidate_id already present, it is not added again."""
    cs = _candidate_set(["tx1", "tx2"])
    # Add a fusion candidate whose ID collides with tx1
    duplicate = _fusion_candidate("tx1")

    result = merge_fusion_candidates(cs, [duplicate])

    assert result.num_candidates == 2
    ids = result.candidate_ids()
    assert ids.count("tx1") == 1


def test_merge_empty_fusion_list():
    """Passing [] returns the candidate_set unchanged."""
    cs = _candidate_set(["tx1", "tx2", "tx3"])
    original_len = cs.num_candidates

    result = merge_fusion_candidates(cs, [])

    assert result is cs
    assert result.num_candidates == original_len


def test_merge_preserves_fusion_source():
    """Fusion candidates in result still have source='fusion'."""
    cs = _candidate_set(["tx1"])
    fc = _fusion_candidate("fus1")

    result = merge_fusion_candidates(cs, [fc])

    fus = result.get_candidate("fus1")
    assert fus is not None
    assert fus.source == "fusion"
