from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.candidates.discovery import _collapse_candidates


def _mk_fusion(cid, junction_pos):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=junction_pos,
        sequence="ACGT" * 100,
        source="fusion",
        supporting_read_ids=set(),
        chrom="chr1::chr2",
        strand="+",
        start=0,
        end=400,
        fusion_junction=(junction_pos, 5000),
        breakpoint_left=("chr1", junction_pos, "+"),
        breakpoint_right=("chr2", 5000, "+"),
    )


def _mk_gtf(cid, three_prime, introns=()):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=introns),
        three_prime_pos=three_prime,
        sequence="ACGT" * 100,
        source="gtf",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=400,
    )


def test_transcript_candidate_has_fusion_fields():
    c = _mk_fusion("fusion_abc", 1000)
    assert c.fusion_junction == (1000, 5000)
    assert c.breakpoint_left == ("chr1", 1000, "+")
    assert c.breakpoint_right == ("chr2", 5000, "+")
    assert c.source == "fusion"


def test_default_fusion_fields_none_for_gtf():
    c = _mk_gtf("tx1", 500)
    assert c.fusion_junction is None
    assert c.breakpoint_left is None
    assert c.breakpoint_right is None


def test_fusion_candidates_not_collapsed():
    # Two DIFFERENT fusions with identical empty intron_chain and close 3' positions
    # Under OLD (pre-fix) logic they would collapse. With guard they must not.
    f1 = _mk_fusion("fusion_1", 1000)
    f2 = _mk_fusion("fusion_2", 1010)  # within threshold of f1's three_prime_pos
    collapsed = _collapse_candidates([f1, f2], threshold=24)
    assert len(collapsed) == 2
    ids = {c.candidate_id for c in collapsed}
    assert ids == {"fusion_1", "fusion_2"}


def test_gtf_candidates_still_collapse_when_identical():
    # Two GTF candidates with same intron chain and close 3' position DO collapse
    g1 = _mk_gtf("tx_a", 1000, introns=((100, 200),))
    g2 = _mk_gtf("tx_b", 1010, introns=((100, 200),))
    collapsed = _collapse_candidates([g1, g2], threshold=24)
    # Still collapses under existing rule (same chain, within threshold)
    assert len(collapsed) == 1


def test_mixed_fusion_and_gtf_preserved():
    g = _mk_gtf("tx_x", 2000, introns=((500, 600),))
    f = _mk_fusion("fusion_x", 1000)
    collapsed = _collapse_candidates([g, f], threshold=24)
    assert len(collapsed) == 2
    sources = {c.source for c in collapsed}
    assert sources == {"gtf", "fusion"}
