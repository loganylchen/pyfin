"""Tests for finalized-model read-supported junction correction."""
from __future__ import annotations

from collections import Counter

from fin.analysis.quantification import QuantResult
from fin.pipeline.junction_snap import snap_quant_results


def _qr(
    cid: str,
    exons,
    *,
    source: str = "novel",
    abundance: float = 1.0,
    confidence: float = 0.5,
    assigned=(),
    max_r: float = 0.5,
    family_id=None,
) -> QuantResult:
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=confidence,
        num_assigned_reads=len(assigned),
        source=source,
        chrom="chr1",
        strand="+",
        start=exons[0][0],
        end=exons[-1][1],
        exons=tuple(exons),
        assigned_read_ids=tuple(assigned),
        max_R=max_r,
        family_id=family_id,
    )


def test_snap_uses_zero_based_introns_and_merges_mass():
    weak = _qr(
        "weak",
        ((100, 200), (300, 400)),
        abundance=1.0,
        confidence=0.5,
        assigned=("r1",),
        max_r=0.6,
        family_id="fam_b",
    )
    strong = _qr(
        "strong",
        ((100, 203), (297, 400)),
        abundance=3.0,
        confidence=1.0,
        assigned=("r2",),
        max_r=1.0,
        family_id="fam_a",
    )
    guided = _qr(
        "guided",
        ((100, 200), (300, 400)),
        source="gtf",
        abundance=2.0,
        assigned=("g1",),
    )
    guided_duplicate = _qr(
        "guided_duplicate",
        ((100, 200), (300, 400)),
        source="gtf",
        abundance=1.0,
        assigned=("g2",),
    )
    observed = {
        ("chr1", "+"): Counter({(200, 300): 1, (203, 297): 6})
    }

    out, snapped, merged = snap_quant_results(
        {
            qr.candidate_id: qr
            for qr in (weak, strong, guided, guided_duplicate)
        },
        observed,
        tolerance=6,
        min_support=2,
        min_ratio=2.0,
    )

    assert snapped == 1
    assert merged == 1
    assert set(out) == {"strong", "guided", "guided_duplicate"}
    corrected = out["strong"]
    assert corrected.exons == ((100, 203), (297, 400))
    assert corrected.abundance == 4.0
    assert corrected.assigned_read_ids == ("r1", "r2")
    assert corrected.num_assigned_reads == 2
    assert corrected.confidence == 0.75
    assert corrected.max_R == 1.0
    assert corrected.family_id == "fam_a"
    assert out["guided"].exons == ((100, 200), (300, 400))
    assert out["guided_duplicate"].exons == ((100, 200), (300, 400))


def test_snap_can_return_absorbed_candidate_redirects():
    weak = _qr(
        "weak", ((100, 200), (300, 400)), abundance=1.0, assigned=("r1",)
    )
    strong = _qr(
        "strong", ((100, 203), (297, 400)), abundance=3.0, assigned=("r2",)
    )
    observed = {
        ("chr1", "+"): Counter({(200, 300): 1, (203, 297): 6})
    }
    out, snapped, merged, redirects = snap_quant_results(
        {"weak": weak, "strong": strong},
        observed,
        tolerance=6,
        min_support=2,
        min_ratio=2.0,
        return_redirects=True,
    )
    assert set(out) == {"strong"}
    assert snapped == 1
    assert merged == 1
    assert redirects == {"weak": "strong"}


def test_snap_reuses_canonical_gate_for_target_junction():
    candidate = _qr("candidate", ((100, 200), (300, 400)))
    observed = {
        ("chr1", "+"): Counter({(200, 300): 1, (203, 297): 6})
    }
    genome = list("A" * 500)
    genome[203:205] = "GT"
    genome[295:297] = "AG"

    out, snapped, _ = snap_quant_results(
        {candidate.candidate_id: candidate},
        observed,
        tolerance=6,
        min_support=2,
        min_ratio=2.0,
        genome_fasta={"chr1": "".join(genome)},
        canonical_motifs=("GT-AG",),
        require_canonical=True,
    )
    assert snapped == 1
    assert out["candidate"].exons == ((100, 203), (297, 400))

    out, snapped, _ = snap_quant_results(
        {candidate.candidate_id: candidate},
        observed,
        tolerance=6,
        min_support=2,
        min_ratio=2.0,
        genome_fasta={"chr1": "A" * 500},
        canonical_motifs=("GT-AG",),
        require_canonical=True,
    )
    assert snapped == 0
    assert out["candidate"].exons == candidate.exons


def test_snap_requires_strictly_stronger_target_support():
    candidate = _qr("candidate", ((100, 200), (300, 400)))
    observed = {
        ("chr1", "+"): Counter({(200, 300): 1, (203, 297): 2})
    }

    out, snapped, merged = snap_quant_results(
        {candidate.candidate_id: candidate},
        observed,
        tolerance=6,
        min_support=2,
        min_ratio=2.0,
    )

    assert out["candidate"].exons == candidate.exons
    assert snapped == 0
    assert merged == 0
