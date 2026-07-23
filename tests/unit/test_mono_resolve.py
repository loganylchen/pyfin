"""Post-EM mono-exon resolution (mono_resolve_drops).

Each surviving single-exon candidate's reads are re-resolved against the surviving
multi-exon candidates by strict strand-aware exonic containment (guard 2), only against
candidates NOT already dropped (guard 1); multi-cover folds into the highest-EM-abundance
candidate (guard 3); a mono candidate under the uncovered-read floor is dropped.
"""
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.analysis.quantification import QuantResult
from fin.scoring.m2_junction_nll import mono_resolve_drops, _mono_read_in_exon


def _cand(cid, introns, start, end, strand="+"):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end, sequence="A", source="novel",
        supporting_read_ids=set(), chrom="chr1", strand=strand, start=start, end=end)


def _qr(cid, reads, abundance=None):
    reads = tuple(sorted(reads))
    return QuantResult(
        candidate_id=cid,
        abundance=float(len(reads)) if abundance is None else abundance,
        confidence=1.0, num_assigned_reads=len(reads), source="novel",
        assigned_read_ids=reads)


# multi M: exons [1000,1500) and [2000,3000); intron (1500,2000)
def _setup():
    M = _cand("M", [(1500, 2000)], 1000, 3000)
    MONO = _cand("MO", [], 1000, 3000)
    return M, MONO


class TestMonoReadInExon:
    def test_in_exon(self):
        assert _mono_read_in_exon(1100, 1200, ((1500, 2000),), 1000, 3000)
        assert _mono_read_in_exon(2100, 2200, ((1500, 2000),), 1000, 3000)

    def test_in_intron_or_crossing_false(self):
        assert not _mono_read_in_exon(1600, 1900, ((1500, 2000),), 1000, 3000)  # intron
        assert not _mono_read_in_exon(1400, 2100, ((1500, 2000),), 1000, 3000)  # crosses

    def test_outside_footprint_false(self):
        assert not _mono_read_in_exon(4000, 4100, ((1500, 2000),), 1000, 3000)


class TestMonoResolve:
    def test_covered_reads_fold_into_multi(self):
        M, MONO = self._setup2()
        qm, qmo = _qr("M", ["a"]), _qr("MO", ["r1", "r2", "r3"])
        spans = {"r1": (1100, 1200), "r2": (2100, 2200), "r3": (1600, 1900)}  # r3 intron
        drops = mono_resolve_drops(
            [M, MONO], [qm, qmo], [10.0, 0.0], spans,
            exclude=set(), slop_bp=10, min_reads=1)
        assert drops == set()                       # MONO keeps r3 (>= 1)
        assert set(qm.assigned_read_ids) == {"a", "r1", "r2"}
        assert qm.abundance == 3.0                  # 1 + 2 folded
        assert set(qmo.assigned_read_ids) == {"r3"}
        assert qmo.abundance == 1.0

    def test_mono_dropped_under_min_reads(self):
        M, MONO = self._setup2()
        qm, qmo = _qr("M", []), _qr("MO", ["r1", "r2", "r3"])
        spans = {"r1": (1100, 1200), "r2": (2100, 2200), "r3": (1600, 1900)}
        drops = mono_resolve_drops(
            [M, MONO], [qm, qmo], [10.0, 0.0], spans,
            exclude=set(), slop_bp=10, min_reads=2)
        assert drops == {1}                         # MONO keeps only r3 (< 2)

    def test_multi_cover_picks_highest_abundance(self):
        M1 = _cand("M1", [(1500, 2000)], 1000, 3000)
        M2 = _cand("M2", [(1500, 2000)], 1000, 3000)
        MONO = _cand("MO", [], 1000, 3000)
        q1, q2, qmo = _qr("M1", []), _qr("M2", []), _qr("MO", ["r1", "u"])
        spans = {"r1": (1100, 1200), "u": (4000, 4100)}   # r1 covered by both, u outside
        mono_resolve_drops([M1, M2, MONO], [q1, q2, qmo], [10.0, 5.0, 0.0], spans,
                           exclude=set(), slop_bp=10, min_reads=1)
        assert "r1" in q1.assigned_read_ids          # higher abundance wins
        assert "r1" not in q2.assigned_read_ids

    def test_guard1_excluded_multi_not_targeted(self):
        M, MONO = self._setup2()
        qm, qmo = _qr("M", []), _qr("MO", ["r1", "r2"])
        spans = {"r1": (1100, 1200), "r2": (2100, 2200)}
        drops = mono_resolve_drops(
            [M, MONO], [qm, qmo], [10.0, 0.0], spans,
            exclude={0}, slop_bp=10, min_reads=1)   # M (col 0) already dropped
        assert set(qmo.assigned_read_ids) == {"r1", "r2"}   # nothing folded
        assert qm.assigned_read_ids == ()

    def test_guard2_intron_crossing_and_strand(self):
        M, MONO = self._setup2()
        MONO_minus = _cand("MO", [], 1000, 3000, strand="-")   # opposite strand
        qm, qmo = _qr("M", []), _qr("MO", ["r1"])
        spans = {"r1": (1100, 1200)}                # in exon but strand mismatches M(+)
        mono_resolve_drops([M, MONO_minus], [qm, qmo], [10.0, 0.0], spans,
                           exclude=set(), slop_bp=10, min_reads=1)
        assert set(qmo.assigned_read_ids) == {"r1"}         # strand guard keeps it

    def test_read_without_span_stays(self):
        M, MONO = self._setup2()
        qm, qmo = _qr("M", []), _qr("MO", ["r1", "nospan"])
        spans = {"r1": (1100, 1200)}
        mono_resolve_drops([M, MONO], [qm, qmo], [10.0, 0.0], spans,
                           exclude=set(), slop_bp=10, min_reads=1)
        assert "nospan" in qmo.assigned_read_ids
        assert "r1" in qm.assigned_read_ids

    def _setup2(self):
        return _setup()
