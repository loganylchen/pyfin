"""Unit tests for the soft-mass / hard-read ratio filter.

Covers ``fin.analysis.quantification.soft_mass_ratio_drops`` — drops a NOVEL
multi-exon candidate whose EM soft abundance (R.sum) is inflated far above its
hard argmax read count, the signature of a wobble shadow that borrows fractional
soft mass from a high-abundance structural near-copy.
"""
from __future__ import annotations

from fin.analysis.quantification import QuantResult, soft_mass_ratio_drops


def _qr(
    cid,
    abundance,
    num_reads,
    *,
    source="novel",
    chrom="chr1",
    strand="+",
    start=0,
    end=1000,
    exons=None,
):
    if exons is None:
        exons = ((start, 400), (600, end))  # 2 exons => multi-exon
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=1.0,
        num_assigned_reads=num_reads,
        source=source,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        exons=exons,
    )


def _drops(results, thr):
    return soft_mass_ratio_drops({qr.candidate_id: qr for qr in results}, thr)


def test_inflated_soft_mass_shadow_dropped():
    # ratio = 187 / 61 = 3.07 >= 2.0 -> dropped (the 32%-abundance shadow case).
    shadow = _qr("shadow", 187.0, 61)
    assert _drops([shadow], 2.0) == {"shadow"}


def test_honest_isoform_kept():
    # ratio = 24 / 24 = 1.0 < 2.0 -> kept (soft mass tracks hard reads).
    real = _qr("real", 24.0, 24)
    assert _drops([real], 2.0) == set()


def test_zero_hard_reads_always_dropped():
    # A candidate living entirely on borrowed soft mass (0 hard reads).
    ghost = _qr("ghost", 12.0, 0)
    assert _drops([ghost], 2.0) == {"ghost"}


def test_boundary_ratio_at_threshold_dropped():
    # ratio == threshold is a drop (>= comparison).
    edge = _qr("edge", 4.0, 2)  # 4/2 = 2.0
    assert _drops([edge], 2.0) == {"edge"}


def test_just_below_threshold_kept():
    near = _qr("near", 3.8, 2)  # 1.9 < 2.0
    assert _drops([near], 2.0) == set()


def test_gtf_source_exempt():
    # A GTF passthrough with an inflated ratio is never dropped.
    gtf = _qr("gtf", 200.0, 10, source="gtf")
    assert _drops([gtf], 2.0) == set()


def test_fusion_source_exempt():
    fus = _qr("fus", 200.0, 10, source="fusion")
    assert _drops([fus], 2.0) == set()


def test_mono_exon_exempt():
    # Single-exon candidate (no internal junction) is exempt.
    mono = _qr("mono", 200.0, 10, exons=((0, 1000),))
    assert _drops([mono], 2.0) == set()


def test_threshold_zero_disables():
    shadow = _qr("shadow", 187.0, 61)
    assert _drops([shadow], 0.0) == set()


def test_mixed_locus_only_shadows_dropped():
    real = _qr("real", 100.0, 101)       # ratio ~0.99
    shadow = _qr("shadow", 50.0, 5)      # ratio 10.0
    gtf = _qr("anchor", 600.0, 12, source="gtf")  # exempt
    assert _drops([real, shadow, gtf], 2.0) == {"shadow"}
