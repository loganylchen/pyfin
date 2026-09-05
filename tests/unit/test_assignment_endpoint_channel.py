"""The 3'/polyA channel inside the assignment posterior (spec A8 / the `tes`
assignment channel). A polyA-PASS read whose 3' end lies past an annotated
model's claimed TES must prefer the model that ends where it does; reads
without a PASS call abstain; the row's best cell stays at 0."""
from __future__ import annotations

import numpy as np
import pytest

from fin.analysis.sample_nuisance import EndpointModel
from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.candidates.discovery import _chain_cluster_candidates
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.assignment import MISSING, apply_endpoint_likelihood


def _cand(cid, start, end, strand="+", source="novel"):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=()), three_prime_pos=end,
        sequence="ACGT" * 30, source=source, supporting_read_ids=set(), chrom="c",
        strand=strand, start=start, end=end)


def _ctx(polya, ends):
    return {"endpoint": EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.04, n_fit=500),
            "polya": polya, "ends": ends}


class TestApplyEndpointLikelihood:
    def test_readthrough_prefers_the_longer_model(self):
        trunc = _cand("anno_trunc", 1000, 1900, source="gtf")   # claims TES 1900
        full = _cand("novel_full", 1000, 2000)                  # ends at 2000
        d = np.full((2, 2), MISSING); d[:, :] = 0.0             # both tied, flat
        ties = {0: [0, 1], 1: [0, 1]}
        ctx = _ctx({"r_pass": (40.0, "PASS"), "r_fail": (5.0, "FAIL")},
                   {"r_pass": ("c", "+", 2000), "r_fail": ("c", "+", 2000)})
        n = apply_endpoint_likelihood(d, ties, ["r_pass", "r_fail"], [trunc, full], ctx)
        assert n == 1
        assert d[0, 1] == 0.0 and d[0, 0] > 5.0        # PASS read: 100 bp past trunc -> penalised
        assert d[1, 0] == 0.0 and d[1, 1] == 0.0        # FAIL read: untouched (abstains)

    def test_minus_strand(self):
        trunc = _cand("anno_trunc", 1100, 2000, strand="-", source="gtf")   # TES = start 1100
        full = _cand("novel_full", 1000, 2000, strand="-")                  # TES 1000
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = _ctx({"r": (40.0, "PASS")}, {"r": ("c", "-", 1000)})
        apply_endpoint_likelihood(d, ties, ["r"], [trunc, full], ctx)
        assert d[0, 1] == 0.0 and d[0, 0] > 5.0

    def test_same_tes_is_neutral_and_missing_cells_untouched(self):
        a = _cand("a", 1000, 2000); b = _cand("b", 1000, 2000); c = _cand("c", 1000, 2500)
        d = np.array([[0.0, 0.0, MISSING]]); ties = {0: [0, 1]}
        ctx = _ctx({"r": (40.0, "PASS")}, {"r": ("c", "+", 2000)})
        apply_endpoint_likelihood(d, ties, ["r"], [a, b, c], ctx)
        assert d[0, 0] == 0.0 and d[0, 1] == 0.0 and d[0, 2] == MISSING

    def test_unavailable_model_is_a_no_op(self):
        a = _cand("a", 1000, 2000); b = _cand("b", 1000, 1900)
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = {"endpoint": EndpointModel(1.0, 80.0, 0.04, n_fit=1), "polya": {"r": (40.0, "PASS")},
               "ends": {"r": ("c", "+", 2000)}}
        assert apply_endpoint_likelihood(d, ties, ["r"], [a, b], ctx) == 0
        assert (d == 0).all()


GENOME = "ACGT" * 3000


class _Tx:
    def __init__(self, tid, exons):
        self.transcript_id, self.exons, self.strand = tid, list(exons), "+"
        self.start, self.end = exons[0][0], exons[-1][1]

    def sort_features(self):
        self.exons.sort()

    def get_spliced_sequence(self, _g):
        return "A" * sum(e - s for s, e in self.exons)


class _Reader:
    def __init__(self, txs):
        self.txs = txs

    def get_transcripts_in_region(self, *_):
        return list(self.txs)


class TestMonoReadExtentHypothesis:
    def _reads(self):
        return [{"query_name": f"m{i}", "reference_start": 5000, "reference_end": 5400,
                 "cigartuples": [(0, 400)], "is_forward": True} for i in range(6)]

    def test_default_absorbs_reads_into_the_overlapping_annotation(self):
        reads = self._reads()
        cs = _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=4900, end=5500, strand="+"), reads,
            {r["query_name"] for r in reads}, {}, _Reader([_Tx("anno_mono", [(5000, 5300)])]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families")
        assert [c.source for c in cs.candidates] == ["gtf"]

    def test_channel_keeps_the_read_extent_when_the_reads_contradict_the_annotation(self):
        reads = self._reads()
        cs = _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=4900, end=5500, strand="+"), reads,
            {r["query_name"] for r in reads}, {}, _Reader([_Tx("anno_mono", [(5000, 5300)])]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=lambda gc, rids, alts=(): {"tes": sorted(rids)})
        ids = sorted(c.candidate_id for c in cs.candidates)
        assert ids == ["anno_mono", "anno_mono.ext"]
        sib = next(c for c in cs.candidates if c.extent_of == "anno_mono")
        assert sib.source == "gtf"                                 # the annotation's trust
        assert (sib.start, sib.end) == (5000, 5400)               # the reads' extent
        assert len(cs.clusters) == 1 and len(cs.clusters[0]) == 2   # same locus, both compete

    def test_consistent_reads_do_not_spawn_a_duplicate(self):
        reads = self._reads()
        cs = _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=4900, end=5500, strand="+"), reads,
            {r["query_name"] for r in reads}, {}, _Reader([_Tx("anno_mono", [(5000, 5300)])]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=lambda gc, rids, alts=(): {})
        assert [c.source for c in cs.candidates] == ["gtf"]


class TestMonoEveryOverlappingAnnotationIsExamined:
    def test_short_model_is_repaired_even_when_the_long_one_wins_overlap(self):
        # a=[1000,1900) and b=[1000,5000); reads at 2000 and 5000. The pooled
        # span attaches the group to b, yet a is the contradicted model.
        reads = [{"query_name": f"s{i}", "reference_start": 1000, "reference_end": 2000,
                  "cigartuples": [(0, 1000)], "is_forward": True} for i in range(8)]
        reads += [{"query_name": f"l{i}", "reference_start": 1000, "reference_end": 5000,
                   "cigartuples": [(0, 4000)], "is_forward": True} for i in range(8)]
        seen = {}
        def hook(gc, rids, alts=()):
            seen[gc.candidate_id] = sorted(a.candidate_id for a in alts)
            if gc.candidate_id == "a":
                return {"tes": [r for r in rids if r.startswith("s")]}
            return {}
        cs = _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=900, end=5100, strand="+"), reads,
            {r["query_name"] for r in reads}, {},
            _Reader([_Tx("a", [(1000, 1900)]), _Tx("b", [(1000, 5000)])]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=hook)
        assert seen == {"a": ["b"], "b": ["a"]}
        sib = next(c for c in cs.candidates if c.extent_of == "a")
        assert (sib.start, sib.end) == (1000, 2000)               # the short reads' TES, not 5000
        assert not any(c.extent_of == "b" for c in cs.candidates)


class TestMultiExonReadExtentHypothesis:
    """The c_trunc failure mode: an annotated key whose reads all carry a WOBBLE
    of its chain (minimap2 shifts non-canonical SIRV junctions, and the canonical
    gate later removes those read variants) and whose 3' end is 100 bp short of
    where every read ends. The reads' extent must be a separate hypothesis with
    the ANNOTATED chain, so it -- not a wobble sibling that happens to carry the
    full span -- wins the reads in the posterior. The evidence for the
    contradiction is the family's wobble-equivalent reads, not the key's own."""

    def _reads(self, n=8, junc=(1200, 1500)):
        rs, re = 1000, 2000
        return [{"query_name": f"m{i}", "reference_start": rs, "reference_end": re,
                 "cigartuples": [(0, junc[0] - rs), (3, junc[1] - junc[0]), (0, re - junc[1])],
                 "is_forward": True} for i in range(n)]

    def _run(self, hook, anno_exons):
        reads = self._reads()
        return _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=900, end=2100, strand="+"), reads,
            {r["query_name"] for r in reads}, {}, _Reader([_Tx("anno", anno_exons)]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=hook)

    def test_contradicted_truncated_annotation_gets_a_read_extent_sibling(self):
        seen = {}
        def hook(gc, rids, alts=()):
            seen[gc.candidate_id] = set(rids)
            return {"tes": sorted(rids)}
        # annotation: same chain (wobbled by 2 bp -> folded under the annotated key),
        # 3' end 1900 while every read ends at 2000
        cs = self._run(hook, [(1000, 1202), (1500, 1900)])
        by_id = {c.candidate_id: c for c in cs.candidates}
        assert seen["anno"] == {f"m{i}" for i in range(8)}          # wobble reads are the evidence
        sib = by_id["anno.ext"]
        assert sib.source == "gtf" and sib.extent_of == "anno"
        assert tuple(sib.intron_chain.introns) == ((1202, 1500),)  # the annotation's chain
        assert (sib.start, sib.end) == (1000, 2000)                # 3' from the reads, 5' annotated
        assert sib.family_id == by_id["anno"].family_id            # same family -> compete
        assert len(cs.clusters) == 1 and set(cs.clusters[0]) >= {"anno", "anno.ext"}

    def test_same_chain_alternatives_are_offered_to_the_hook(self):
        got = {}
        def hook(gc, rids, alts=()):
            got[gc.candidate_id] = sorted(a.candidate_id for a in alts)
            return {}
        reads = self._reads()
        _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=900, end=2100, strand="+"), reads,
            {r["query_name"] for r in reads}, {},
            _Reader([_Tx("short", [(1000, 1202), (1500, 1900)]), _Tx("long", [(1000, 1202), (1500, 2000)])]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=hook)
        assert got == {"short": ["long"], "long": ["short"]}

    def test_sibling_id_never_shadows_a_real_annotation_id(self):
        reads = self._reads()
        reader = _Reader([_Tx("anno", [(1000, 1202), (1500, 1900)]),
                          _Tx("anno.ext", [(1000, 1202), (1500, 2000)])])
        reader.transcripts = {"anno": None, "anno.ext": None}
        cs = _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=900, end=2100, strand="+"), reads,
            {r["query_name"] for r in reads}, {}, reader,
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=lambda gc, rids, alts=(): {"tes": sorted(rids)} if gc.candidate_id == "anno" else {})
        ids = sorted(c.candidate_id for c in cs.candidates if c.source == "gtf")
        assert ids == ["anno", "anno.ext", "anno.ext_"]
        assert next(c for c in cs.candidates if c.candidate_id == "anno.ext_").extent_of == "anno"

    def test_degraded_reads_keep_the_annotated_5p_end(self):
        # Reads start INSIDE exon 2 (5' degradation) and only the 3' end is
        # contradicted: the sibling must keep the annotated TSS, so its span
        # still contains the first intron (a min(read start) would not).
        reads = [{"query_name": f"m{i}", "reference_start": 1550, "reference_end": 2000,
                  "cigartuples": [(0, 450)], "is_forward": True} for i in range(4)]
        reads += self._reads(n=2)                                 # two full-length reads
        cs = _chain_cluster_candidates(
            GenomicInterval(chrom="chr1", start=900, end=2100, strand="+"), reads,
            {r["query_name"] for r in reads}, {}, _Reader([_Tx("anno", [(1000, 1200), (1500, 1900)])]),
            GENOME, "+", 1, wobble_bp=6, cassette_max_exon_bp=70, clustering="families",
            read_extent_contradicted=lambda gc, rids, alts=(): {"tes": sorted(rids)})
        sib = next(c for c in cs.candidates if c.extent_of == "anno")
        assert (sib.start, sib.end) == (1000, 2000)
        assert sib.start < sib.intron_chain.introns[0][0]

    def test_only_the_contradicted_end_moves(self):
        from fin.candidates.discovery import _read_extent
        plus = _cand("a", 1000, 1900, source="gtf"); minus = _cand("b", 1000, 1900, strand="-", source="gtf")
        sp = {"r": (900, 2000)}
        assert _read_extent(plus, {"tes": ["r"]}, sp) == (1000, 2000)
        assert _read_extent(plus, {"tss": ["r"]}, sp) == (900, 1900)
        assert _read_extent(minus, {"tes": ["r"]}, sp) == (900, 1900)   # TES is start on '-'
        assert _read_extent(minus, {"tss": ["r"]}, sp) == (1000, 2000)
        assert _read_extent(plus, {"tes": ["r"]}, {"r": (1100, 1850)}) == (1000, 1900)  # never shrinks

    def test_only_contradicting_reads_define_the_end(self):
        # annotated APA: reads explained by the long model "b" (ending 5000)
        # must not set the repaired TES of the short model "a" (reads at 2000).
        from fin.candidates.discovery import _read_extent
        a = _cand("a", 1000, 1900, source="gtf")
        spans = {"s1": (1000, 2000), "s2": (1000, 2001), "l1": (1000, 5000)}
        assert _read_extent(a, {"tes": ["s1", "s2"]}, spans) == (1000, 2001)

    def test_uncontradicted_annotation_stays_alone(self):
        cs = self._run(lambda gc, rids, alts=(): {}, [(1000, 1202), (1500, 1900)])
        assert not any(c.extent_of for c in cs.candidates)

    def test_hook_off_is_the_baseline(self):
        a = self._run(None, [(1000, 1202), (1500, 1900)])
        b = self._run(lambda gc, rids, alts=(): {}, [(1000, 1202), (1500, 1900)])
        assert [(c.candidate_id, c.start, c.end) for c in a.candidates] == \
               [(c.candidate_id, c.start, c.end) for c in b.candidates]


class TestStartTerm:
    def test_upstream_start_penalises_the_truncated_model(self):
        from fin.analysis.sample_nuisance import StartModel
        trunc = _cand("anno_trunc5p", 1100, 2000, source="gtf")   # claims TSS 1100
        full = _cand("novel_full", 1000, 2000)                   # starts at 1000
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = {"endpoint": None, "polya": {}, "ends": {},
               "start": StartModel(b_narrow=0.5, b_broad=500.0, eta=0.05, n_fit=1000),
               "starts": {"r": ("c", "+", 1000)}}                 # read begins at 1000
        n = apply_endpoint_likelihood(d, ties, ["r"], [trunc, full], ctx)
        assert n == 1
        assert d[0, 1] == 0.0 and d[0, 0] > 5.0        # 100 bp upstream of trunc's TSS

    def test_downstream_start_abstains(self):
        from fin.analysis.sample_nuisance import StartModel
        a = _cand("a", 1000, 2000); b = _cand("b", 1100, 2000)
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = {"endpoint": None, "polya": {}, "ends": {},
               "start": StartModel(0.5, 500.0, 0.05, n_fit=1000),
               "starts": {"r": ("c", "+", 1500)}}                 # degraded: inside both
        apply_endpoint_likelihood(d, ties, ["r"], [a, b], ctx)
        assert d[0, 0] == 0.0 and d[0, 1] == 0.0

    def test_minus_strand_upstream(self):
        from fin.analysis.sample_nuisance import StartModel
        trunc = _cand("t", 1000, 1900, strand="-", source="gtf")   # TSS = end 1900
        full = _cand("f", 1000, 2000, strand="-")                  # TSS 2000
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = {"endpoint": None, "polya": {}, "ends": {},
               "start": StartModel(0.5, 500.0, 0.05, n_fit=1000),
               "starts": {"r": ("c", "-", 2000)}}
        apply_endpoint_likelihood(d, ties, ["r"], [trunc, full], ctx)
        assert d[0, 1] == 0.0 and d[0, 0] > 5.0


class TestExtentPairResolution:
    """One model, two extent estimates: the posterior picks one, the other folds."""

    def _qo(self, mass_anno, mass_ext):
        from types import SimpleNamespace
        from fin.analysis.quantification import QuantResult
        anno = _cand("anno", 1000, 1900, source="gtf")
        ext = _cand("anno.ext", 1000, 2000, source="gtf"); ext.extent_of = "anno"
        other = _cand("novel_x", 3000, 4000)
        R = np.zeros((mass_anno + mass_ext + 1, 3))
        R[:mass_anno, 0] = 1.0; R[mass_anno:mass_anno + mass_ext, 1] = 1.0; R[-1, 2] = 1.0
        def qr(c, j):
            rows = [i for i in range(R.shape[0]) if R[i, j] == 1.0]
            return QuantResult(candidate_id=c.candidate_id, abundance=float(len(rows)), confidence=1.0,
                               num_assigned_reads=len(rows), source=c.source, assigned_read_ids=tuple(f"r{i}" for i in rows),
                               max_R=1.0 if rows else 0.0)
        cands = [anno, ext, other]
        return cands, R, [qr(c, j) for j, c in enumerate(cands)]

    def test_read_extent_wins_and_absorbs(self):
        from fin.pipeline.selection import resolve_extent_pairs
        cands, R, qrs = self._qo(2, 5)
        drop, out = set(), []
        assert resolve_extent_pairs(cands, R, {}, {}, qrs, drop, out)[0] == 1
        assert drop == {0} and out[0].action == "fold" and out[0].fold_into == "anno.ext"
        assert qrs[1].num_assigned_reads == 7 and qrs[1].abundance == 7.0
        assert qrs[2].num_assigned_reads == 1                      # unrelated column untouched

    def test_confidence_is_read_weighted_after_pooling(self):
        from fin.pipeline.selection import resolve_extent_pairs
        cands, R, qrs = self._qo(2, 5)
        qrs[0].confidence, qrs[1].confidence = 1.0, 0.3
        resolve_extent_pairs(cands, R, {}, {}, qrs, set(), [])
        assert qrs[1].confidence == pytest.approx((0.3 * 5 + 1.0 * 2) / 7)

    def test_tie_keeps_the_annotation(self):
        from fin.pipeline.selection import resolve_extent_pairs
        cands, R, qrs = self._qo(3, 3)
        drop, out = set(), []
        resolve_extent_pairs(cands, R, {}, {}, qrs, drop, out)
        assert drop == {1} and out[0].fold_into == "anno" and qrs[0].num_assigned_reads == 6

    def test_no_siblings_is_a_no_op(self):
        from fin.pipeline.selection import resolve_extent_pairs
        cands, R, qrs = self._qo(1, 1)
        cands[1].extent_of = None
        drop, out = set(), []
        assert resolve_extent_pairs(cands, R, {}, {}, qrs, drop, out)[0] == 0 and not drop and not out


class TestContradictionAlternatives:
    """runner._read_extent_contradicted: a read explained by another annotated
    model's end (annotated APA / alternative TSS) contradicts nothing."""

    def _runner(self, tmp_path):
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        cfg = PipelineConfig(bam_path="/fake.bam", work_dir=str(tmp_path),
                             assignment_channels="ends", signal_path="/fake.blow5")
        r = PipelineRunner(cfg)
        polya = {f"r{i}": [40.0, "PASS"] for i in range(120)}
        ends = {f"r{i}": ["c", "+", 2000 + (i % 3) - 1, 1000 + (i % 3) - 1] for i in range(120)}
        (tmp_path / r._SAMPLE_EVIDENCE_NAME).write_text(json.dumps(
            {"polya": polya, "ends": ends, "version": "sample-evidence-v2"}))
        r._annotated_tes = lambda: {("c", "+"): [2000]}
        r._annotated_tss = lambda: {("c", "+"): [1000]}
        return r

    def test_readthrough_explained_by_annotated_apa_is_not_a_contradiction(self, tmp_path):
        r = self._runner(tmp_path)
        hook = r._read_extent_contradicted()
        short = _cand("short", 1000, 1900, source="gtf")
        long_ = _cand("long", 1000, 2000, source="gtf")
        reads = [f"r{i}" for i in range(120)]          # all end at ~2000
        assert sorted(hook(short, reads)) == ["tes"] and len(hook(short, reads)["tes"]) == 120
        assert not hook(short, reads, [long_])            # the long model explains them
        assert not hook(long_, reads, [short])            # consistent with its own end

    def test_stale_work_dir_evidence_is_rebuilt_not_trusted(self, tmp_path):
        import json
        r = self._runner(tmp_path)
        path = tmp_path / r._SAMPLE_EVIDENCE_NAME
        assert json.loads(path.read_text()).get("identity") is None
        # a file without the input identity is not accepted as-is: it is rebuilt
        r.config.signal_path = None
        r._prepare_sample_evidence()
        assert json.loads(path.read_text())["identity"]["bam"].endswith("fake.bam")


class TestPooledMassVisibleToLaterGates:
    def test_resolve_returns_pooled_matrices(self):
        from fin.pipeline.selection import resolve_extent_pairs
        from fin.analysis.quantification import QuantResult
        anno = _cand("anno", 1000, 1900, source="gtf")
        ext = _cand("anno.ext", 1000, 2000, source="gtf"); ext.extent_of = "anno"
        R = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
        qrs = [QuantResult("anno", 2.0, 1.0, 2, "gtf", assigned_read_ids=("r0", "r1")),
               QuantResult("anno.ext", 3.0, 1.0, 3, "gtf", assigned_read_ids=("r2", "r3", "r4"))]
        ties = {0: [0, 1], 2: [1], 4: [0]}
        nlls = {0: {0: 0.5, 1: 0.2}, 2: {1: 0.1}}
        n, R2, t2, n2 = resolve_extent_pairs([anno, ext], R, ties, nlls, qrs, set(), [])
        assert n == 1
        assert R2.sum(axis=0).tolist() == [0.0, 5.0]           # winner carries all mass
        assert R.sum(axis=0).tolist() == [2.0, 3.0]            # input untouched
        assert t2 == {0: [1], 2: [1], 4: [1]}
        assert n2 == {0: {1: 0.2}, 2: {1: 0.1}}


class TestConfigGuard:
    def test_channels_require_m2_em(self, tmp_path):
        from fin.pipeline.config import PipelineConfig
        bam, sig = tmp_path / "r.bam", tmp_path / "r.blow5"
        bam.write_bytes(b""); sig.write_bytes(b"")
        with pytest.raises(ValueError, match="m2_em"):
            PipelineConfig(bam_path=str(bam), signal_path=str(sig),
                           assignment_channels="tes", quant_mode="argmax").validate()
        PipelineConfig(bam_path=str(bam), signal_path=str(sig),
                       assignment_channels="tes", quant_mode="m2_em").validate()


class TestFoldedLoserIsInvisibleDownstream:
    """After pair resolution the loser must not feed the opt-in decision
    layers (containment evidence / junction election), nor the aggregate."""

    def test_containment_layer_receives_live_models_only(self, tmp_path, monkeypatch):
        from types import SimpleNamespace
        import fin.analysis.containment_evidence as ce
        from fin.analysis.quantification import QuantResult
        from fin.pipeline import selection as S
        anno = _cand("anno", 1100, 2000, source="gtf")
        ext = _cand("anno.ext", 1000, 2000, source="gtf"); ext.extent_of = "anno"
        R = np.array([[1.0, 0.0]] * 6 + [[0.0, 1.0]] * 4)
        qrs = [QuantResult("anno", 6.0, 1.0, 6, "gtf", assigned_read_ids=tuple(f"r{i}" for i in range(6))),
               QuantResult("anno.ext", 4.0, 1.0, 4, "gtf", assigned_read_ids=tuple(f"r{i}" for i in range(6, 10)))]
        qo = SimpleNamespace(cand_list=[anno, ext], R=R, ties_by_read={}, nlls_by_read={},
                             quant_results=qrs, kept_read_ids=[f"r{i}" for i in range(10)])
        seen = {}
        def fake_eval(cands, quant, spans, interval_name=""):
            seen["ids"] = [c.candidate_id for c in cands]
            return []
        monkeypatch.setattr(ce, "evaluate_containment", fake_eval)
        cfg = SimpleNamespace(containment_evidence_mode="audit", junction_election_mode="off",
                              m2_cluster_recheck=False, novel_junction_min_reads=0,
                              guided_junction_min_reads=0, m2_support_gate=False,
                              containment_cluster=False, containment_collapse=False,
                              mono_resolve_post_em=False, work_dir=str(tmp_path))
        cs = SimpleNamespace(read_spans={})
        iv = GenomicInterval(chrom="c", start=0, end=3000, strand="+")
        res, out = S.select_m2_interval(cfg, qo, cs, iv, lambda i: {})
        assert seen["ids"] == ["anno"]                         # the folded sibling is gone
        assert [q.candidate_id for q in res] == ["anno"] and res[0].num_assigned_reads == 10

    def test_unify_relabels_the_loser_before_aggregation(self):
        from fin.analysis.quantification import QuantResult, aggregate_across_intervals
        from fin.pipeline.selection import unify_extent_pairs
        # interval A kept the annotation (mass 2), interval B its sibling (mass 5)
        a = QuantResult("a", 2.0, 1.0, 2, "gtf", chrom="c", strand="+", start=1000, end=1900,
                        exons=((1000, 1900),), assigned_read_ids=("l1", "l2"))
        e = QuantResult("a.ext", 5.0, 1.0, 5, "gtf", chrom="c", strand="+", start=1000, end=2000,
                        exons=((1000, 2000),), assigned_read_ids=tuple(f"r{i}" for i in range(5)))
        e.extent_of = "a"
        results = [[a], [e]]
        redirects, out = unify_extent_pairs(results)
        assert redirects == {"a": "a.ext"} and out[0].action == "fold" and out[0].fold_into == "a.ext"
        assert a.candidate_id == "a.ext" and a.extent_of == "a"      # relabelled in place
        agg = aggregate_across_intervals(results)
        assert set(agg) == {"a.ext"}
        assert agg["a.ext"].num_assigned_reads == 7 and agg["a.ext"].abundance == pytest.approx(7.0)
        assert (agg["a.ext"].start, agg["a.ext"].end) == (1000, 2000)

    def test_unify_shared_reads_are_deduplicated_by_the_aggregator(self):
        from fin.analysis.quantification import QuantResult, aggregate_across_intervals
        from fin.pipeline.selection import unify_extent_pairs
        a = QuantResult("a", 2.0, 1.0, 2, "gtf", chrom="c", strand="+", start=1000, end=1900,
                        exons=((1000, 1900),), assigned_read_ids=("r1", "r2"))
        e = QuantResult("a.ext", 3.0, 1.0, 3, "gtf", chrom="c", strand="+", start=1000, end=2000,
                        exons=((1000, 2000),), assigned_read_ids=("r1", "r2", "r3"))
        e.extent_of = "a"
        results = [[a], [e]]
        unify_extent_pairs(results)
        agg = aggregate_across_intervals(results)["a.ext"]
        assert agg.num_assigned_reads == 3 and agg.abundance == pytest.approx(3.0)

    def test_unify_annotation_wins_and_keeps_its_geometry(self):
        from fin.analysis.quantification import QuantResult, aggregate_across_intervals
        from fin.pipeline.selection import unify_extent_pairs
        a = QuantResult("a", 6.0, 1.0, 6, "gtf", chrom="c", strand="+", start=1000, end=1900,
                        exons=((1000, 1900),), assigned_read_ids=tuple(f"l{i}" for i in range(6)))
        e = QuantResult("a.ext", 2.0, 1.0, 2, "gtf", chrom="c", strand="+", start=1000, end=2000,
                        exons=((1000, 2000),), assigned_read_ids=("r1", "r2"))
        e.extent_of = "a"
        results = [[e], [a]]                                          # sibling first
        redirects, _ = unify_extent_pairs(results)
        assert redirects == {"a.ext": "a"} and e.candidate_id == "a" and e.extent_of is None
        agg = aggregate_across_intervals(results)["a"]
        assert (agg.start, agg.end, agg.exons) == (1000, 1900, ((1000, 1900),))
        assert agg.num_assigned_reads == 8

    def test_unify_is_a_no_op_without_split_pairs(self):
        from fin.analysis.quantification import QuantResult
        from fin.pipeline.selection import unify_extent_pairs
        x = QuantResult("x", 1.0, 1.0, 1, "novel")
        assert unify_extent_pairs([[x]]) == ({}, []) and x.candidate_id == "x"

    def test_global_fold_is_a_refit_redirect(self):
        from fin.analysis.abundance_refit import ResponsibilityLedger, refit_survivor_abundance
        from fin.analysis.quantification import QuantResult
        # interval A kept "a" (2 reads), interval B kept "a.ext" (3 reads); the
        # global fold sent "a" -> "a.ext": refit must route A's reads there.
        la = ResponsibilityLedger(interval_key="A", input_read_ids=("l1", "l2"),
                                  weights={"l1": {"a": 1.0}, "l2": {"a": 1.0}}, redirects={}, forced={})
        lb = ResponsibilityLedger(interval_key="B", input_read_ids=("r1", "r2", "r3"),
                                  weights={r: {"a.ext": 1.0} for r in ("r1", "r2", "r3")}, redirects={}, forced={})
        e = QuantResult("a.ext", 5.0, 1.0, 5, "gtf", exons=((0, 100),))
        out, diag = refit_survivor_abundance({"a.ext": e}, [la, lb], snap_redirects={"a": "a.ext"})
        assert out["a.ext"].abundance == pytest.approx(5.0)
        out2, _ = refit_survivor_abundance({"a.ext": QuantResult("a.ext", 5.0, 1.0, 5, "gtf", exons=((0, 100),))}, [la, lb])
        assert out2["a.ext"].abundance == pytest.approx(3.0)      # without the redirect the reads orphan

    def test_sibling_geometry_across_intervals_is_order_independent(self):
        from fin.analysis.quantification import QuantResult, aggregate_across_intervals
        def mk(start, end, rid):
            q = QuantResult("a.ext", 1.0, 1.0, 1, "gtf", chrom="c", strand="+", start=start, end=end,
                            exons=((start, 1200), (1500, end)), assigned_read_ids=(rid,))
            q.extent_of = "a"
            return q
        left, right = mk(900, 5000, "l"), mk(1000, 5100, "r")
        g1 = aggregate_across_intervals([[left], [right]])["a.ext"]
        g2 = aggregate_across_intervals([[right], [left]])["a.ext"]
        for g in (g1, g2):
            assert (g.start, g.end) == (900, 5100) and g.exons == ((900, 1200), (1500, 5100))


class TestKernelAvailability:
    def test_bimodal_toy_sample_abstains(self):
        from fin.analysis.sample_nuisance import fit_endpoint_model
        m = fit_endpoint_model([0] * 8 + [100] * 8)
        assert not m.available                    # 16 ends, half 100 bp off: not a jitter kernel

    def test_sparse_sample_abstains_and_calibrated_sample_is_available(self):
        import random
        from fin.analysis.sample_nuisance import MIN_KERNEL_FIT, fit_endpoint_model
        rng = random.Random(1)
        d = [int(round(rng.expovariate(1.0))) * rng.choice([-1, 1]) for _ in range(400)]
        d += [int(round(rng.expovariate(1 / 80.0))) * rng.choice([-1, 1]) for _ in range(16)]
        assert fit_endpoint_model(d).available
        assert not fit_endpoint_model(d[:MIN_KERNEL_FIT - 1]).available
