"""Unit tests for the minimum-isoform-fraction (locus-relative abundance) filter.

Covers ``fin.analysis.quantification.isoform_fraction_drops`` — the production
port of the ``benchmarks/diag_fp_gtf_features.py`` relabund lever (Cufflinks
``--min-isoform-fraction`` / StringTie ``-f`` minor-isoform suppression).
"""
from __future__ import annotations

import pytest

from fin.analysis.quantification import QuantResult, isoform_fraction_drops


def _qr(
    cid,
    abundance,
    *,
    source="novel",
    chrom="chr1",
    strand="+",
    start=0,
    end=1000,
    exons=None,
    family_id=None,
):
    if exons is None:
        exons = ((start, 400), (600, end))  # 2 exons => multi-exon
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=1.0,
        num_assigned_reads=int(abundance),
        source=source,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        exons=exons,
        family_id=family_id,
    )


def _drops(results, thr, locus="family"):
    return isoform_fraction_drops(
        {qr.candidate_id: qr for qr in results}, thr, locus=locus
    )


def test_family_mode_keeps_overlap_from_another_family():
    dominant = _qr("dom", 100.0, family_id="fam_a")
    minor = _qr("min", 0.5, family_id="fam_b")
    assert _drops([dominant, minor], 0.01) == set()


def test_family_mode_drops_minor_inside_same_family():
    dominant = _qr("dom", 100.0, family_id="fam_a")
    minor = _qr("min", 0.5, family_id="fam_a")
    assert _drops([dominant, minor], 0.01) == {"min"}


def test_attached_gtf_contributes_to_family_maximum():
    gtf = _qr("gtf", 100.0, source="gtf", family_id="fam_a")
    minor = _qr("min", 0.5, family_id="fam_a")
    assert _drops([gtf, minor], 0.01) == {"min"}


def test_gtf_from_another_family_does_not_suppress_novel():
    gtf = _qr("gtf", 100.0, source="gtf", family_id="fam_a")
    minor = _qr("min", 0.5, family_id="fam_b")
    assert _drops([gtf, minor], 0.01) == set()


def test_overlap_mode_restores_historical_denominator():
    dominant = _qr("dom", 100.0, family_id="fam_a")
    minor = _qr("min", 0.5, family_id="fam_b")
    assert _drops([dominant, minor], 0.01, locus="overlap") == {"min"}


def test_invalid_locus_mode_is_rejected():
    with pytest.raises(ValueError, match="locus"):
        _drops([_qr("a", 1.0)], 0.01, locus="gene")


def test_minor_overlapping_isoform_below_fraction_is_dropped():
    dominant = _qr("dom", 100.0)
    minor = _qr("min", 0.5)  # 0.5 / 100 = 0.005 < 0.01
    assert _drops([dominant, minor], 0.01) == {"min"}


def test_minor_at_or_above_fraction_is_kept():
    dominant = _qr("dom", 100.0)
    minor = _qr("min", 2.0)  # 2 / 100 = 0.02 >= 0.01
    assert _drops([dominant, minor], 0.01) == set()


def test_dominant_isoform_never_dropped():
    # The dominant isoform has relabund == 1.0 by construction.
    dominant = _qr("dom", 100.0)
    minor = _qr("min", 0.5)
    drops = _drops([dominant, minor], 0.01)
    assert "dom" not in drops


def test_gtf_and_fusion_sources_are_exempt():
    dominant = _qr("dom", 100.0)
    gtf_minor = _qr("gtf_min", 0.1, source="gtf")
    fusion_minor = _qr("fus_min", 0.1, source="fusion")
    assert _drops([dominant, gtf_minor, fusion_minor], 0.01) == set()


def test_mono_exon_novel_is_exempt():
    dominant = _qr("dom", 100.0)
    mono = _qr("mono", 0.1, exons=((0, 1000),))  # single exon
    assert _drops([dominant, mono], 0.01) == set()


def test_non_overlapping_candidates_are_each_their_own_locus():
    # Different chromosome => no overlap => each is locus-dominant (relabund 1).
    a = _qr("a", 100.0, chrom="chr1")
    b = _qr("b", 0.5, chrom="chr2")
    assert _drops([a, b], 0.01) == set()


def test_opposite_strand_does_not_form_a_locus():
    a = _qr("a", 100.0, strand="+")
    b = _qr("b", 0.5, strand="-")
    assert _drops([a, b], 0.01) == set()


def test_zero_fraction_disables_filter():
    dominant = _qr("dom", 100.0)
    minor = _qr("min", 0.0001)
    assert _drops([dominant, minor], 0.0) == set()


def test_higher_threshold_drops_more():
    dominant = _qr("dom", 100.0)
    minor = _qr("min", 5.0)  # relabund 0.05
    assert _drops([dominant, minor], 0.01) == set()
    assert _drops([dominant, minor], 0.10) == {"min"}


def test_min_fraction_above_one_never_drops_dominant():
    # Defensive: a mis-configured fraction > 1.0 must NOT drop the locus
    # dominant (relabund == 1.0). Non-dominant minors may still be dropped.
    dominant = _qr("dom", 100.0)
    minor = _qr("min", 50.0)
    drops = _drops([dominant, minor], 1.5)
    assert "dom" not in drops


def test_three_isoform_locus_uses_max_as_denominator():
    dominant = _qr("dom", 100.0)
    mid = _qr("mid", 50.0)  # 0.5 -> kept
    minor = _qr("min", 0.5)  # 0.005 -> dropped
    assert _drops([dominant, mid, minor], 0.01) == {"min"}


# --- all-source denominator (the wobble-shadow fix) -----------------------
def test_novel_shadow_dropped_against_high_gtf_at_locus():
    # The locus dominant is a GTF transcript (read support in the hundreds).
    # The novel shadow only competes with another novel; without the GTF in the
    # denominator its relabund (4/8 = 0.5) would keep it. With GTF included
    # (4/800 = 0.005) it is correctly dropped.
    gtf = _qr("gtf_dom", 800.0, source="gtf")
    novel_strong = _qr("nov_strong", 40.0)   # 40/800 = 0.05 >= 0.01 -> kept
    novel_shadow = _qr("nov_shadow", 4.0)     # 4/800 = 0.005 < 0.01 -> dropped
    drops = _drops([gtf, novel_strong, novel_shadow], 0.01)
    assert "nov_shadow" in drops          # 0.005 < 0.01 vs GTF max
    assert "nov_strong" not in drops      # 0.05 >= 0.01 -> kept
    assert "gtf_dom" not in drops         # GTF never dropped


def test_gtf_max_keeps_novel_at_or_above_fraction():
    # A novel at exactly the fraction vs the GTF max survives (>= keeps).
    gtf = _qr("gtf_dom", 800.0, source="gtf")
    novel = _qr("nov", 8.0)               # 8/800 = 0.01 -> kept (>=)
    assert _drops([gtf, novel], 0.01) == set()


def test_gtf_only_locus_partner_still_drops_lone_novel_shadow():
    # Even if the novel is the ONLY novel at the locus, a high GTF partner can
    # now drop it (previously a lone novel was always locus-max -> never dropped).
    gtf = _qr("gtf_dom", 500.0, source="gtf")
    lone_shadow = _qr("shadow", 2.0)      # 2/500 = 0.004 < 0.01
    assert _drops([gtf, lone_shadow], 0.01) == {"shadow"}


def test_low_novel_with_no_overlapping_gtf_is_kept():
    # Recall-safety: a low-abundance novel with no overlapping higher transcript
    # is locus-dominant and must NOT be dropped.
    lone = _qr("lone", 3.0)               # only candidate at its locus
    assert _drops([lone], 0.01) == set()
