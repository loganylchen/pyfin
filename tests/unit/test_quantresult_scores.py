import pytest

from fin.analysis.quantification import QuantResult, aggregate_across_intervals


def test_quantresult_defaults():
    qr = QuantResult(
        candidate_id="tx1",
        abundance=1.0,
        confidence=0.9,
        num_assigned_reads=3,
        source="gtf",
    )
    assert qr.coherence_score == 0.0
    assert qr.discrimination_score == 0.0
    assert qr.combined_score == 0.0
    assert qr.breakpoint_left is None
    assert qr.breakpoint_right is None


def test_quantresult_with_scores():
    qr = QuantResult(
        candidate_id="tx1",
        abundance=1.0,
        confidence=0.9,
        num_assigned_reads=3,
        source="gtf",
        coherence_score=0.85,
        discrimination_score=0.7,
        combined_score=0.77,
    )
    assert qr.coherence_score == 0.85
    assert qr.discrimination_score == 0.7
    assert qr.combined_score == 0.77


def test_quantresult_with_breakpoints():
    qr = QuantResult(
        candidate_id="fusion_1",
        abundance=5.0,
        confidence=0.8,
        num_assigned_reads=5,
        source="fusion",
        breakpoint_left=("chr1", 1000, "+"),
        breakpoint_right=("chr2", 5000, "+"),
    )
    assert qr.breakpoint_left == ("chr1", 1000, "+")
    assert qr.breakpoint_right == ("chr2", 5000, "+")


def test_aggregate_preserves_scores():
    qr1 = QuantResult(
        candidate_id="tx1",
        abundance=2.0, confidence=0.9, num_assigned_reads=2, source="gtf",
        coherence_score=0.5, discrimination_score=0.6, combined_score=0.55,
    )
    qr2 = QuantResult(
        candidate_id="tx1",
        abundance=3.0, confidence=0.8, num_assigned_reads=3, source="gtf",
        coherence_score=0.7, discrimination_score=0.4, combined_score=0.6,
    )
    agg = aggregate_across_intervals([[qr1], [qr2]])
    assert "tx1" in agg
    r = agg["tx1"]
    assert r.abundance == 5.0
    # Read-count-weighted average: (0.5*2 + 0.7*3) / 5 = 0.62
    assert r.coherence_score == pytest.approx(0.62)
    # (0.6*2 + 0.4*3) / 5 = 0.48
    assert r.discrimination_score == pytest.approx(0.48)
    # (0.55*2 + 0.6*3) / 5 = 0.58
    assert r.combined_score == pytest.approx(0.58)


def test_aggregate_zero_reads_falls_back_to_unit_weight():
    """When all intervals have zero assigned reads, scores degrade to a simple mean."""
    qr1 = QuantResult(
        candidate_id="tx1",
        abundance=0.0, confidence=0.0, num_assigned_reads=0, source="gtf",
        coherence_score=0.4, discrimination_score=0.6, combined_score=0.5,
    )
    qr2 = QuantResult(
        candidate_id="tx1",
        abundance=0.0, confidence=0.0, num_assigned_reads=0, source="gtf",
        coherence_score=0.8, discrimination_score=0.2, combined_score=0.4,
    )
    agg = aggregate_across_intervals([[qr1], [qr2]])
    r = agg["tx1"]
    # Weight=1.0 fallback for both -> simple mean.
    assert r.coherence_score == pytest.approx(0.6)
    assert r.discrimination_score == pytest.approx(0.4)
    assert r.combined_score == pytest.approx(0.45)
    # Confidence has no fallback -> stays 0.0 when no reads.
    assert r.confidence == 0.0


def test_aggregate_mixed_zero_and_nonzero_reads():
    """Intervals with zero reads use weight=1.0; positive-read intervals weight by count."""
    qr1 = QuantResult(
        candidate_id="tx1",
        abundance=0.0, confidence=0.0, num_assigned_reads=0, source="gtf",
        coherence_score=0.0, discrimination_score=0.0, combined_score=0.0,
    )
    qr2 = QuantResult(
        candidate_id="tx1",
        abundance=4.0, confidence=0.7, num_assigned_reads=4, source="gtf",
        coherence_score=1.0, discrimination_score=1.0, combined_score=1.0,
    )
    agg = aggregate_across_intervals([[qr1], [qr2]])
    r = agg["tx1"]
    # Total weight = 1.0 + 4.0 = 5.0; (0*1 + 1.0*4) / 5 = 0.8
    assert r.coherence_score == pytest.approx(0.8)
    assert r.discrimination_score == pytest.approx(0.8)
    assert r.combined_score == pytest.approx(0.8)


def test_aggregate_preserves_breakpoints():
    qr = QuantResult(
        candidate_id="fusion_1",
        abundance=3.0, confidence=0.7, num_assigned_reads=3, source="fusion",
        breakpoint_left=("chrA", 100, "+"),
        breakpoint_right=("chrB", 500, "-"),
    )
    agg = aggregate_across_intervals([[qr]])
    r = agg["fusion_1"]
    assert r.breakpoint_left == ("chrA", 100, "+")
    assert r.breakpoint_right == ("chrB", 500, "-")
