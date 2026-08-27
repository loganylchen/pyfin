"""Functional smoke test of the clustering="families" adapter in _chain_cluster_candidates.

Exercises the new path end-to-end at the discovery level (no BAM I/O): cluster_families + collapse
+ mono-bucket routing → the shared emission loop → a CandidateSet with candidates / clusters /
shadows. Confirms the 5' degradation ladder folds into the full-length and mono reads still emit.
"""
from fin.candidates.discovery import _chain_cluster_candidates
from fin.io.interval_manager import GenomicInterval

GENOME = "ACGT" * 3000  # 12 kb dummy reference (coords below stay < len)


def _read(qname, ref_start, blocks):
    """blocks = list of (match_len, intron_len); last intron_len=0 ends the read.
    Builds cigartuples (op 0 = M, op 3 = N) and reference_start/_end."""
    cig = []
    pos = ref_start
    for mlen, ilen in blocks:
        cig.append((0, mlen)); pos += mlen
        if ilen:
            cig.append((3, ilen)); pos += ilen
    return {"query_name": qname, "reference_start": ref_start, "reference_end": pos,
            "cigartuples": cig, "is_forward": True}


def test_families_folds_ladder_and_emits_mono():
    # full chain ((1100,2000),(3000,4000)): two reads.
    full = [(100, 900), (1000, 1000), (100, 0)]        # 1000..1100 M, 1100..2000 N, ... etc.
    r1 = _read("f1", 1000, full)
    r2 = _read("f2", 1000, full)
    assert r1["reference_end"] == 4100
    # 5' truncation: chain ((3000,4000),) only -> exact suffix sub-chain of the full.
    trunc = [(100, 1000), (100, 0)]                    # 2900..3000 M, 3000..4000 N, 4000..4100 M
    r3 = _read("t1", 2900, trunc)
    # a single-exon (mono) read -> holding bucket.
    mono = {"query_name": "m1", "reference_start": 5000, "reference_end": 5300,
            "cigartuples": [(0, 300)], "is_forward": True}

    reads = [r1, r2, r3, mono]
    ivl = GenomicInterval(chrom="chr1", start=900, end=5400, strand="+")
    read_ids = {r["query_name"] for r in reads}

    cs = _chain_cluster_candidates(
        ivl, reads, read_ids, {}, None, GENOME, "+", 1,
        wobble_bp=6, cassette_max_exon_bp=70, clustering="families")

    novel = [c for c in cs.candidates if c.source == "novel"]
    by_chain = {tuple(c.intron_chain.introns): c for c in novel}
    full_chain = ((1100, 2000), (3000, 4000))
    assert full_chain in by_chain                                    # full-length emitted
    full_cand = by_chain[full_chain]
    assert full_cand.family_id is not None
    # the 5' truncation read folded into the full-length candidate (reads pooled).
    assert full_cand.supporting_read_ids == {"f1", "f2", "t1"}
    # the folded sub-chain is recorded as a shadow for downstream 5'-TSS recovery.
    shadows = cs.shadows.get(full_cand.candidate_id, [])
    assert any(sc == ((3000, 4000),) for sc, _ in shadows)
    # the truncation did NOT emit its own standalone candidate.
    assert ((3000, 4000),) not in by_chain
    # the mono read emitted a family-less single-exon candidate.
    mono_candidates = [
        c for c in novel if not c.intron_chain.introns and "m1" in c.supporting_read_ids
    ]
    assert len(mono_candidates) == 1
    assert mono_candidates[0].family_id is None


def test_related_variants_share_stable_family_id():
    exact = _read("exact", 1000, [(100, 900), (1000, 1000), (100, 0)])
    wobble = _read("wobble", 1000, [(102, 898), (1000, 1000), (100, 0)])
    unrelated = _read("other", 5900, [(100, 500), (100, 0)])
    reads = [exact, wobble, unrelated]
    ivl = GenomicInterval(chrom="chr1", start=900, end=6700, strand="+")

    cs = _chain_cluster_candidates(
        ivl, reads, {r["query_name"] for r in reads}, {}, None, GENOME, "+", 1,
        wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
    )
    by_chain = {tuple(c.intron_chain.introns): c for c in cs.candidates}
    exact_id = by_chain[((1100, 2000), (3000, 4000))].family_id
    wobble_id = by_chain[((1102, 2000), (3000, 4000))].family_id
    other_id = by_chain[((6000, 6500),)].family_id
    assert exact_id is not None and exact_id == wobble_id
    assert other_id is not None and other_id != exact_id

    reversed_cs = _chain_cluster_candidates(
        ivl, list(reversed(reads)), {r["query_name"] for r in reads}, {}, None,
        GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70,
        clustering="families",
    )
    reversed_ids = {
        tuple(c.intron_chain.introns): c.family_id for c in reversed_cs.candidates
    }
    assert reversed_ids == {chain: cand.family_id for chain, cand in by_chain.items()}


class _Transcript:
    def __init__(self, transcript_id, exons):
        self.transcript_id = transcript_id
        self.exons = list(exons)
        self.strand = "+"
        self.start = exons[0][0]
        self.end = exons[-1][1]

    def sort_features(self):
        self.exons.sort()

    def get_spliced_sequence(self, _genome):
        return "A" * sum(end - start for start, end in self.exons)


class _GTFReader:
    def __init__(self, transcripts):
        self.transcripts = transcripts

    def get_transcripts_in_region(self, _chrom, _start, _end):
        return list(self.transcripts)


def test_attached_gtf_and_novel_share_family_id():
    exact = _read("exact", 1000, [(100, 900), (1000, 1000), (100, 0)])
    wobble = _read("wobble", 1000, [(102, 898), (1000, 1000), (100, 0)])
    gtf = _Transcript("gtf_exact", [(1000, 1100), (2000, 3000), (4000, 4100)])
    ivl = GenomicInterval(chrom="chr1", start=900, end=4200, strand="+")

    cs = _chain_cluster_candidates(
        ivl, [exact, wobble], {"exact", "wobble"}, {}, _GTFReader([gtf]),
        GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70,
        clustering="families",
    )
    gtf_candidate = next(c for c in cs.candidates if c.source == "gtf")
    novel_candidate = next(c for c in cs.candidates if c.source == "novel")
    assert gtf_candidate.family_id is not None
    assert gtf_candidate.family_id == novel_candidate.family_id


def test_families_fold_monoexon_folds_contained_mono():
    # with fold_monoexon on (the non-m2_em path), a mono read wholly inside a multi member's exon
    # folds into it instead of emitting a standalone mono candidate.
    full = [(100, 900), (1000, 1000), (100, 0)]        # chain ((1100,2000),(3000,4000)), span 1000-4100
    r1 = _read("f1", 1000, full)
    mono = {"query_name": "m1", "reference_start": 2100, "reference_end": 2400,  # inside exon [2000,3000]
            "cigartuples": [(0, 300)], "is_forward": True}
    reads = [r1, mono]
    ivl = GenomicInterval(chrom="chr1", start=900, end=4200, strand="+")
    cs = _chain_cluster_candidates(
        ivl, reads, {"f1", "m1"}, {}, None, GENOME, "+", 1,
        wobble_bp=6, cassette_max_exon_bp=70, fold_monoexon_contained=True, clustering="families")
    novel = [c for c in cs.candidates if c.source == "novel"]
    assert not any(not c.intron_chain.introns for c in novel)      # NO standalone mono candidate
    multi = [c for c in novel if c.intron_chain.introns]
    assert len(multi) == 1 and multi[0].supporting_read_ids == {"f1", "m1"}   # mono folded in
