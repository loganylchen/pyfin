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


def _tiny_bam(path, reads=((1000, False), (1200, True))):
    """Two primary alignments on contig "c" (+ and -), indexed."""
    import pysam
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 10000, "SN": "c"}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i, (start, rev) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"; a.flag = 16 if rev else 0
            a.reference_id = 0; a.reference_start = start; a.mapping_quality = 60
            a.cigarstring = "500M"; a.query_sequence = "A" * 500
            a.query_qualities = pysam.qualitystring_to_array("I" * 500)
            fh.write(a)
    pysam.index(str(path))
    return path


class TestContradictionAlternatives:
    """runner._read_extent_contradicted: a read explained by another annotated
    model's end (annotated APA / alternative TSS) contradicts nothing."""

    def _runner(self, tmp_path):
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        cfg = PipelineConfig(bam_path=str(_tiny_bam(tmp_path / "r.bam")),
                             work_dir=str(tmp_path),
                             assignment_channels="ends", signal_path="/fake.blow5")
        r = PipelineRunner(cfg)
        polya = {f"r{i}": [40.0, "PASS"] for i in range(120)}
        ends = {f"r{i}": ["c", "+", 2000 + (i % 3) - 1, 1000 + (i % 3) - 1] for i in range(120)}
        (r._evidence_path()).write_text(json.dumps(
            {"polya": polya, "ends": ends, "version": "sample-evidence-v3",
             "identity": r._evidence_identity(), "polya_status": "ok",
             "n_ends": len(ends), "n_polya": len(polya)}))
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
        path = r._evidence_path()
        stale = json.loads(path.read_text())
        stale["identity"]["bam"] = "/some/other.bam:1:2:3:edge:0000000000000000"
        path.write_text(json.dumps(stale))
        # the evidence no longer belongs to this run's inputs: this run must not
        # consume it, and must produce its own
        r.config.signal_path = None
        r.config.sample_evidence_path = ""
        r._invalidate_sample_caches()
        if hasattr(r, "_evidence_identity_cache"):
            del r._evidence_identity_cache
        r._prepare_sample_evidence()
        d = json.loads(r._evidence_path().read_text())
        assert "r.bam:" in d["identity"]["bam"]
        assert d["version"] == "sample-evidence-v3" and d["polya_status"] == "no_signal"
        assert set(d["ends"]) == {"r0", "r1"}          # built from the real BAM

    def test_regenerated_bam_at_the_same_path_is_a_cache_miss(self, tmp_path):
        import json
        import time
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                             assignment_channels="ends", quant_mode="m2_em")
        r = PipelineRunner(cfg)
        r._prepare_sample_evidence()
        first = json.loads(r._evidence_path().read_text())
        time.sleep(0.01)
        _tiny_bam(tmp_path / "r.bam", reads=((4000, False),))     # same path, new content
        PipelineRunner(cfg)._prepare_sample_evidence()
        second = json.loads(PipelineRunner(cfg)._evidence_path().read_text())
        assert second["identity"] != first["identity"]
        assert set(second["ends"]) == {"r0"} and second["ends"]["r0"][2] == 4500

    def test_failed_polya_is_not_published_as_final(self, tmp_path, monkeypatch):
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        import fin.scoring.polya as polya_mod
        sig = tmp_path / "s.blow5"; sig.write_bytes(b"x")
        bam = _tiny_bam(tmp_path / "r.bam")
        cache = tmp_path / "cache"; cache.mkdir()
        monkeypatch.setenv("PYFIN_SAMPLE_EVIDENCE_CACHE", str(cache))
        monkeypatch.setattr(polya_mod, "compute_polya", lambda *a, **k: {})
        cfg = PipelineConfig(bam_path=str(bam), signal_path=str(sig),
                             work_dir=str(tmp_path / "w"),
                             assignment_channels="ends", quant_mode="m2_em")
        PipelineRunner(cfg)._prepare_sample_evidence()
        d = json.loads(PipelineRunner(cfg)._evidence_path().read_text())
        assert d["polya_status"] == "unavailable" and d["ends"]        # ends still usable
        assert not list(cache.glob("sample_evidence_*.json"))          # never cached as final
        # ... and a later run retries instead of inheriting the abstention
        calls = []
        monkeypatch.setattr(polya_mod, "compute_polya",
                            lambda seqs, *a, **k: calls.append(1) or {r: (40.0, "PASS") for r in seqs})
        PipelineRunner(cfg)._prepare_sample_evidence()
        d2 = json.loads(PipelineRunner(cfg)._evidence_path().read_text())
        assert calls and d2["polya_status"] == "ok" and d2["polya"]

    def test_malformed_records_are_rejected_with_matching_headers(self, tmp_path):
        """Record-level validation, exercised PAST the header checks: a payload
        with the right version, identity and status must still be rejected when
        a record has the wrong shape -- otherwise it silently disables a channel
        or fabricates a coordinate."""
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        work = tmp_path / "w"; work.mkdir()
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(work),
                             assignment_channels="ends", quant_mode="m2_em")
        r = PipelineRunner(cfg)
        ident = r._evidence_identity()
        path = r._evidence_path()

        def payload(ends, polya=None):
            polya = {} if polya is None else polya
            return json.dumps({"version": "sample-evidence-v3", "identity": ident,
                               "polya_status": "no_signal", "ends": ends,
                               "polya": polya,
                               "n_ends": len(ends), "n_polya": len(polya)})

        good = {"r0": ["c", "+", 1500, 1000]}
        path.write_text(payload(good))
        assert PipelineRunner._read_sample_evidence(path, ident, "no_signal") is not None
        for ends, polya in [
            ({"r0": ["c", "+", 1500]}, None),                 # no 5' field
            ({"r0": "c+12"}, None),                           # string, not a record
            ({"r0": [["c"], "+", 1500, 1000]}, None),         # chrom not a string
            ({"r0": ["c", "?", 1500, 1000]}, None),           # bad strand
            ({"r0": ["c", "+", 1500, float("inf")]}, None),   # non-finite coordinate
            ({"r0": ["c", "+", 1500, 1000.5]}, None),         # fractional coordinate
            ({"r0": ["c", "+", 1500, 10 ** 400]}, None),      # outside the double range
            (good, {"r0": ["bad", "PASS"]}),                  # polyA length not a number
            (good, {"r0": ["PASS"]}),                         # polyA record too short
            (good, {"r0": [10 ** 400, "PASS"]}),              # length outside the double range
        ]:
            path.write_text(payload(ends, polya))
            assert PipelineRunner._read_sample_evidence(path, ident, "no_signal") is None, (ends, polya)

    def test_evidence_from_another_sample_is_not_consumed(self, tmp_path):
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        other = _tiny_bam(tmp_path / "other.bam", reads=((9000, False),))
        work = tmp_path / "w"; work.mkdir()
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(work),
                             assignment_channels="ends", quant_mode="m2_em")
        foreign = PipelineRunner(PipelineConfig(
            bam_path=str(other), work_dir=str(work),
            assignment_channels="ends", quant_mode="m2_em"))
        foreign._prepare_sample_evidence()                    # writes OTHER sample's evidence
        d = json.loads(foreign._evidence_path().read_text())
        assert set(d["ends"]) == {"r0"} and d["ends"]["r0"][2] == 9500
        # the evidence is addressed by input identity, so the two runs cannot
        # collide in a shared work dir, and this consumer sees no evidence at
        # all rather than the other sample's
        mine = PipelineRunner(cfg)
        assert mine._evidence_path() != foreign._evidence_path()
        assert mine._sample_evidence() is None
        mine._prepare_sample_evidence()                       # both files coexist
        assert foreign._evidence_path().exists() and mine._evidence_path().exists()
        assert set(mine._sample_evidence()[1]) == {"r0", "r1"}
        assert set(json.loads(foreign._evidence_path().read_text())["ends"]) == {"r0"}

    def test_malformed_evidence_is_rebuilt_not_trusted(self, tmp_path):
        import dataclasses
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        work = tmp_path / "w"; work.mkdir()
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(work),
                             assignment_channels="ends", quant_mode="m2_em")
        path = PipelineRunner(cfg)._evidence_path()
        for payload in ("[]", "{}", '{"version": "sample-evidence-v3", "ends": {"r0": []}}',
                        '{"version": "sample-evidence-v3", "ends": {}}', "not json"):
            path.write_text(payload)
            assert PipelineRunner._read_sample_evidence(path) is None, payload
            r = PipelineRunner(dataclasses.replace(cfg))
            r._prepare_sample_evidence()                       # never crashes
            # the corrupt file is NOT overwritten (another run may be reading
            # it); this run works from its own copy
            assert path.read_text() == payload
            assert "_run" in r._evidence_path().name
            assert r._sample_evidence()[1]["r0"][2] == 1500
            r.cleanup()

    def test_fingerprint_sees_an_equal_metadata_replacement(self, tmp_path):
        # size and mtime alone cannot distinguish a file replaced in place; the
        # head/tail digest can, which is what keeps derived evidence honest.
        import os
        from fin.pipeline.runner import input_fingerprint
        a = tmp_path / "in.bin"
        a.write_bytes(b"A" * 4096)
        st = os.stat(a)
        first = input_fingerprint(a)
        a.write_bytes(b"B" * 4096)                      # same size, different content
        os.utime(a, ns=(st.st_atime_ns, st.st_mtime_ns))   # same mtime
        second = input_fingerprint(a)
        assert os.stat(a).st_size == st.st_size
        assert os.stat(a).st_mtime_ns == st.st_mtime_ns
        assert second != first
        assert input_fingerprint(a) == second and input_fingerprint("") == ""

    def test_interior_edit_needs_the_strict_hash(self, monkeypatch, tmp_path):
        """The documented bound: the fast fingerprint samples the edges, so a
        same-length interior edit with restored timestamps can hide from it;
        PYFIN_EVIDENCE_STRICT_HASH=1 hashes the whole file and cannot."""
        import os
        from fin.pipeline.runner import input_fingerprint
        a = tmp_path / "big.bin"
        a.write_bytes(b"A" * 196608)
        st = os.stat(a)
        fast_before = input_fingerprint(a)
        monkeypatch.setenv("PYFIN_EVIDENCE_STRICT_HASH", "1")
        strict_before = input_fingerprint(a)
        monkeypatch.delenv("PYFIN_EVIDENCE_STRICT_HASH")
        with a.open("r+b") as fh:                       # same length, interior byte
            fh.seek(80000); fh.write(b"B")
        os.utime(a, ns=(st.st_atime_ns, st.st_mtime_ns))
        assert os.stat(a).st_size == st.st_size and os.stat(a).st_mtime_ns == st.st_mtime_ns
        monkeypatch.setenv("PYFIN_EVIDENCE_STRICT_HASH", "1")
        assert input_fingerprint(a) != strict_before     # strict mode always sees it
        monkeypatch.delenv("PYFIN_EVIDENCE_STRICT_HASH")
        assert strict_before != fast_before              # the two modes never collide

    def test_regenerated_bam_is_a_cache_miss(self, tmp_path):
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                             assignment_channels="ends", quant_mode="m2_em")
        PipelineRunner(cfg)._prepare_sample_evidence()
        old_path = PipelineRunner(cfg)._evidence_path()
        first = json.loads(old_path.read_text())
        _tiny_bam(tmp_path / "r.bam", reads=((7000, False), (7200, True)))
        r = PipelineRunner(cfg)
        r._reset_input_snapshot()                    # what run() does at its boundary
        r._prepare_sample_evidence()
        assert r._evidence_path() != old_path        # different input -> different file
        second = json.loads(r._evidence_path().read_text())
        assert second["identity"] != first["identity"]
        assert second["ends"]["r0"][2] == 7500
        assert r._sample_evidence()[1]["r0"][2] == 7500   # never the stale endpoints


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


class TestProductionDefaults:
    def test_real_drna_profile_enables_ends_and_tss_require(self):
        from fin.pipeline.config import PIPELINE_PROFILES
        for name in ("real-drna", "real-drna-precision"):
            assert PIPELINE_PROFILES[name]["assignment_channels"] == "ends"
            assert PIPELINE_PROFILES[name]["folded_tss_recovery_mode"] == "require"
        assert "assignment_channels" not in PIPELINE_PROFILES["sirv"]   # benchmark working point unchanged

    def test_channels_without_signal_abstain_instead_of_failing(self, tmp_path, caplog):
        from fin.pipeline.config import PipelineConfig
        bam = tmp_path / "r.bam"; bam.write_bytes(b"")
        gtf = tmp_path / "a.gtf"; gtf.write_text("")
        # ends + annotation: the 3' term abstains, the 5' term still applies
        PipelineConfig(bam_path=str(bam), gtf_path=str(gtf),
                       assignment_channels="ends", quant_mode="m2_em").validate()
        assert any("the 5' start channel still applies" in r.getMessage()
                   for r in caplog.records)
        caplog.clear()
        # no annotation (or tes-only): there is no endpoint evidence at all
        PipelineConfig(bam_path=str(bam), assignment_channels="ends",
                       quant_mode="m2_em").validate()
        assert any("no endpoint evidence" in r.getMessage() for r in caplog.records)

    def test_profile_implied_settings_degrade_under_another_quantifier(self, tmp_path, caplog):
        from fin.pipeline.config import PipelineConfig
        bam = tmp_path / "x.bam"; bam.touch()
        c = PipelineConfig.from_profile("real-drna", bam_path=str(bam), quant_mode="argmax")
        c.validate()
        assert c.assignment_channels == "none" and c.folded_tss_recovery_mode == "off"
        assert "assignment_channels=ends disabled" in caplog.text
        with pytest.raises(ValueError, match="m2_em"):     # explicit choice still errors
            PipelineConfig.from_profile("real-drna", bam_path=str(bam), quant_mode="argmax",
                                        assignment_channels="ends").validate()

    def test_ends_are_collected_without_signal(self, tmp_path):
        """The 5' channel needs only BAM read starts: the prepass must collect
        the strand-aware ends even when no signal file is configured, so the
        promoted default degrades to "3' abstains", not "both channels dead"."""
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                             assignment_channels="ends", quant_mode="m2_em")
        r = PipelineRunner(cfg)
        r._prepare_sample_evidence()
        d = json.loads(r._evidence_path().read_text())
        assert d["polya"] == {}                                  # no signal -> every call MISSING
        assert d["ends"]["r0"] == ["c", "+", 1500, 1000]          # 3' = end, 5' = start on +
        assert d["ends"]["r1"] == ["c", "-", 1200, 1700]          # mirrored on -
        polya_map, read_ends, read_starts = r._sample_evidence()
        assert read_starts["r0"] == ("c", "+", 1000) and read_starts["r1"] == ("c", "-", 1700)


class TestEvidenceLifecycle:
    """The evidence file is an artifact of ONE run over ONE input set."""

    def _cfg(self, tmp_path, bam, **kw):
        from fin.pipeline.config import PipelineConfig
        return PipelineConfig(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                              assignment_channels="ends", quant_mode="m2_em", **kw)

    def test_failed_polya_never_occupies_the_shared_name(self, tmp_path, monkeypatch):
        """Two runs over the SAME inputs, one whose polyA pass fails: the failed
        one must not own the identity-addressed file, or the successful run's
        workers would read an abstaining payload."""
        import json
        import fin.scoring.polya as polya_mod
        from fin.pipeline.runner import PipelineRunner
        sig = tmp_path / "s.blow5"; sig.write_bytes(b"x")
        bam = _tiny_bam(tmp_path / "r.bam")
        monkeypatch.setattr(polya_mod, "compute_polya", lambda *a, **k: {})
        bad = PipelineRunner(self._cfg(tmp_path, bam, signal_path=str(sig)))
        bad._prepare_sample_evidence()
        bad_path = bad._evidence_path()
        assert json.loads(bad_path.read_text())["polya_status"] == "unavailable"
        assert f"_run" in bad_path.name                          # run-local, not shared
        assert bad._sample_evidence()[1]                         # its own ends are usable

        monkeypatch.setattr(polya_mod, "compute_polya",
                            lambda seqs, *a, **k: {r: (40.0, "PASS") for r in seqs})
        cfg_ok = self._cfg(tmp_path, bam, signal_path=str(sig))
        good = PipelineRunner(cfg_ok)
        good._prepare_sample_evidence()               # the failed run did not satisfy it
        good_path = good._evidence_path()
        assert good_path != bad_path and "_run" not in good_path.name
        assert json.loads(good_path.read_text())["polya"]
        assert json.loads(bad_path.read_text())["polya"] == {}   # untouched
        assert PipelineRunner(cfg_ok)._sample_evidence()[0]      # consumers see the calls

    def test_late_failure_leaves_published_evidence_intact(self, tmp_path, monkeypatch):
        """Success first, failure second: the failed writer must not touch the
        artifact the successful run's workers are already reading."""
        import json
        import fin.scoring.polya as polya_mod
        from fin.pipeline.runner import PipelineRunner
        sig = tmp_path / "s.blow5"; sig.write_bytes(b"x")
        bam = _tiny_bam(tmp_path / "r.bam")
        cfg_ok = self._cfg(tmp_path, bam, signal_path=str(sig))
        monkeypatch.setattr(polya_mod, "compute_polya",
                            lambda seqs, *a, **k: {r: (40.0, "PASS") for r in seqs})
        good = PipelineRunner(cfg_ok); good._prepare_sample_evidence()
        shared_path = good._evidence_path()
        before = shared_path.read_text()
        # a concurrent run over the same inputs that cannot call polyA; it must
        # not find the published artifact acceptable for itself NOR replace it
        monkeypatch.setattr(polya_mod, "compute_polya", lambda *a, **k: {})
        late = PipelineRunner(self._cfg(tmp_path, bam, signal_path=str(sig)))
        late.config.sample_evidence_path = ""
        late._prepare_sample_evidence()
        assert shared_path.read_text() == before                 # byte-identical
        assert json.loads(shared_path.read_text())["polya"]

    def test_broken_calls_are_dropped_before_publication(self, tmp_path, monkeypatch):
        """A NaN length is not a failed call, it is a broken one: dropping it at
        the source keeps the artifact parseable, so the read ends survive."""
        import json
        import fin.scoring.polya as polya_mod
        from fin.pipeline.runner import PipelineRunner
        sig = tmp_path / "s.blow5"; sig.write_bytes(b"x")
        bam = _tiny_bam(tmp_path / "r.bam")
        monkeypatch.setattr(polya_mod, "compute_polya", lambda seqs, *a, **k: {
            "r0": (float("nan"), "FAIL"), "r1": (40.0, "PASS")})
        cfg = self._cfg(tmp_path, bam, signal_path=str(sig))
        r = PipelineRunner(cfg); r._prepare_sample_evidence()
        d = json.loads(r._evidence_path().read_text())
        assert set(d["polya"]) == {"r1"}
        # an incomplete producer result is this run's evidence, never the
        # artifact another run inherits for these inputs
        assert d["polya_status"] == "partial" and "_run" in r._evidence_path().name
        polya_map, read_ends, _ = r._sample_evidence()            # artifact stays usable
        assert set(read_ends) == {"r0", "r1"} and set(polya_map) == {"r1"}
        r.cleanup()

    def test_reset_drops_signal_and_annotation_state(self, tmp_path):
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam))
        closed = []

        class _Reader:
            def close(self):
                closed.append(True)

        r._annotated_tss_cache = {("c", "+"): [1]}
        r._annotated_tes_cache = {("c", "+"): [2]}
        r._jc_ctx = ("aligner", "krill", lambda rid: None)   # the real shape: a closure
        r._jc_slow5 = _Reader()                              # the handle it closes over
        r._reset_input_snapshot()
        assert not hasattr(r, "_annotated_tss_cache") and not hasattr(r, "_annotated_tes_cache")
        assert not hasattr(r, "_jc_ctx") and closed == [True]

    def test_run_resets_the_snapshot(self, tmp_path, monkeypatch):
        # the reset must be wired into run(), not only available as a method
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam))
        called = []
        monkeypatch.setattr(PipelineRunner, "_reset_input_snapshot",
                            lambda self: called.append(True))
        monkeypatch.setattr("fin.pipeline.runner.generate_isolated_intervals",
                            lambda *a, **k: {"intervals": [], "fusion_reads": []})
        r.run()
        assert called == [True]

    def test_reused_runner_refingerprints_its_inputs(self, tmp_path, monkeypatch):
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam))
        r._prepare_sample_evidence()
        first = r._evidence_path()
        _tiny_bam(tmp_path / "r.bam", reads=((7000, False),))     # same path, new content
        r._reset_input_snapshot()                                 # what run() does
        r._prepare_sample_evidence()
        assert r._evidence_path() != first
        assert r._sample_evidence()[1]["r0"][2] == 7500

    def test_workers_read_the_recorded_path_not_their_own_guess(self, tmp_path):
        """A worker must consume the file its RUN prepared. Pointing the
        recorded path at a different (valid) artifact proves it is honoured
        rather than re-derived from the worker's own fingerprint."""
        import dataclasses
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        other = _tiny_bam(tmp_path / "other.bam", reads=((9000, False),))
        parent = PipelineRunner(self._cfg(tmp_path, bam))
        parent._prepare_sample_evidence()
        assert parent.config.sample_evidence_path                 # recorded for workers
        elsewhere = PipelineRunner(self._cfg(tmp_path, other))
        elsewhere._prepare_sample_evidence()
        worker_cfg = dataclasses.replace(
            parent.config, sample_evidence_path=str(elsewhere._evidence_path()))
        worker = PipelineRunner(worker_cfg)
        assert worker._evidence_path() == elsewhere._evidence_path()
        assert set(worker._sample_evidence()[1]) == {"r0"}        # the recorded artifact
        # and the ordinary case: same config -> the parent's own file
        plain = PipelineRunner(dataclasses.replace(parent.config))
        assert plain._evidence_path() == parent._evidence_path()
        assert set(plain._sample_evidence()[1]) == {"r0", "r1"}

    def test_run_local_artifact_is_removed_on_cleanup(self, tmp_path, monkeypatch):
        import fin.scoring.polya as polya_mod
        from fin.pipeline.runner import PipelineRunner
        sig = tmp_path / "s.blow5"; sig.write_bytes(b"x")
        bam = _tiny_bam(tmp_path / "r.bam")
        monkeypatch.setattr(polya_mod, "compute_polya", lambda *a, **k: {})
        r = PipelineRunner(self._cfg(tmp_path, bam, signal_path=str(sig)))
        r._prepare_sample_evidence()
        path = r._evidence_path()
        assert "_run" in path.name and path.exists()
        r.cleanup()
        assert not path.exists()                                  # disposable, not left behind

    def test_missing_prepared_evidence_is_reported(self, tmp_path, caplog):
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam))
        r._prepare_sample_evidence()
        r._evidence_path().unlink()
        r._invalidate_sample_caches()
        assert r._sample_evidence() is None
        assert any("is missing" in rec.getMessage() for rec in caplog.records)

    def test_retry_on_one_runner_uses_the_new_evidence(self, tmp_path, monkeypatch):
        """A retry must not adopt the artifact of its own earlier failed
        attempt: each run-local publication takes a fresh token."""
        import json
        import fin.scoring.polya as polya_mod
        from fin.pipeline.runner import PipelineRunner
        sig = tmp_path / "s.blow5"; sig.write_bytes(b"x")
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam, signal_path=str(sig)))
        monkeypatch.setattr(polya_mod, "compute_polya", lambda *a, **k: {})
        r._prepare_sample_evidence()                       # attempt 1: unavailable
        first = r._evidence_path()
        assert json.loads(first.read_text())["polya"] == {}
        # attempt 2 on the SAME runner: one good call, one broken -> partial
        monkeypatch.setattr(polya_mod, "compute_polya", lambda seqs, *a, **k: {
            "r0": (float("nan"), "FAIL"), "r1": (40.0, "PASS")})
        r.config.sample_evidence_path = ""
        r._invalidate_sample_caches()
        r._prepare_sample_evidence()
        second = r._evidence_path()
        assert second != first
        assert set(json.loads(second.read_text())["polya"]) == {"r1"}
        assert set(r._sample_evidence()[0]) == {"r1"}      # the NEW calls, not the old
        r.cleanup()
        assert not first.exists() and not second.exists()  # both attempts cleaned up

    def test_annotation_replaced_during_parse_is_refused(self, tmp_path, monkeypatch):
        import pytest
        from fin.pipeline.runner import PipelineRunner
        gtf = tmp_path / "a.gtf"; gtf.write_text("chr1\tx\ttranscript\t1\t2\t.\t+\t.\t\n")
        bam = _tiny_bam(tmp_path / "r.bam")
        cfg = self._cfg(tmp_path, bam)
        cfg.gtf_path = str(gtf)
        r = PipelineRunner(cfg)
        import fin.io.io_gtf as io_gtf

        real_parse = io_gtf.GTFReader.parse

        def parse_then_swap(self_):
            out = real_parse(self_)
            gtf.write_text("chr1\tx\ttranscript\t9\t9\t.\t+\t.\t\n")   # replaced under us
            return out

        monkeypatch.setattr(io_gtf.GTFReader, "parse", parse_then_swap)
        with pytest.raises(RuntimeError, match="changed while it was being parsed"):
            r._open_gtf_reader()

    def test_parseable_subset_is_rejected(self, tmp_path):
        """A file whose records were removed but whose header still matches is
        a damaged artifact, not a sample with fewer reads."""
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam))
        r._prepare_sample_evidence()
        path = r._evidence_path()
        full = json.loads(path.read_text())
        assert full["n_ends"] == 2
        for damaged in ({}, {"r0": full["ends"]["r0"]}):
            d = dict(full, ends=damaged)                  # n_ends still says 2
            path.write_text(json.dumps(d))
            assert PipelineRunner._read_sample_evidence(path) is None
        path.write_text(json.dumps(full))
        assert PipelineRunner._read_sample_evidence(path) is not None

    def test_input_replaced_mid_scan_is_refused(self, tmp_path, monkeypatch):
        """The identity is taken before the scan; if the BAM is swapped while it
        is being read, the evidence would mix two files under one name."""
        import pytest
        import fin.pipeline.runner as R
        from fin.pipeline.runner import PipelineRunner
        bam = _tiny_bam(tmp_path / "r.bam")
        r = PipelineRunner(self._cfg(tmp_path, bam))
        r._evidence_identity()                            # fingerprint taken
        _tiny_bam(tmp_path / "r.bam", reads=((5000, False),))   # replaced underneath
        with pytest.raises(RuntimeError, match="changed while"):
            r._prepare_sample_evidence()

    def test_failed_annotation_reopen_is_retried(self, tmp_path, monkeypatch):
        from fin.pipeline.runner import PipelineRunner
        gtf = tmp_path / "a.gtf"; gtf.write_text("")
        bam = _tiny_bam(tmp_path / "r.bam")
        cfg = self._cfg(tmp_path, bam)
        cfg.gtf_path = str(gtf)
        r = PipelineRunner(cfg)
        calls = []

        def flaky(self_):
            calls.append(1)
            if len(calls) == 1:
                raise OSError("transient")
            self_._gtf_fingerprint = R_fp(self_.config.gtf_path)
            return "reader"

        import fin.pipeline.runner as R
        R_fp = R.input_fingerprint
        monkeypatch.setattr(PipelineRunner, "_open_gtf_reader", flaky)
        try:
            r._reset_input_snapshot()
        except OSError:
            pass
        assert r._gtf_reader is None and r._gtf_fingerprint is None
        r._reset_input_snapshot()                          # unchanged bytes -> retried
        assert r._gtf_reader == "reader" and len(calls) == 2


class TestFivePrimeClusterRequirement:
    """A read starting upstream of an annotated TSS is evidence about that
    transcript's start only when other molecules end there too -- the 5'
    analogue of the polyA-PASS requirement on the 3' side. Measured on human
    dRNA: 60% of upstream-starting reads are alone at their position."""

    def _hook(self, tmp_path, n_upstream, min_reads=3):
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        cfg = PipelineConfig(bam_path=str(_tiny_bam(tmp_path / "r.bam")),
                             work_dir=str(tmp_path / "w"), signal_path="/fake.blow5",
                             assignment_channels="ends", quant_mode="m2_em",
                             tss_contradiction_min_reads=min_reads)
        r = PipelineRunner(cfg)
        polya = {f"b{i}": [40.0, "PASS"] for i in range(120)}
        # 120 background reads at the annotated TSS/TES, plus n_upstream reads
        # that start 500 bp upstream of it (all at the same position)
        ends = {f"b{i}": ["c", "+", 3000 + (i % 3) - 1, 2000 + (i % 3) - 1]
                for i in range(120)}
        for j in range(n_upstream):
            ends[f"u{j}"] = ["c", "+", 3000, 1500]
            polya[f"u{j}"] = [40.0, "PASS"]
        (tmp_path / "w").mkdir(parents=True, exist_ok=True)
        path = r._evidence_path()
        path.write_text(json.dumps(
            {"polya": polya, "ends": ends, "version": "sample-evidence-v3",
             "identity": r._evidence_identity(), "polya_status": "ok",
             "n_ends": len(ends), "n_polya": len(polya)}))
        r._annotated_tes = lambda: {("c", "+"): [3000]}
        r._annotated_tss = lambda: {("c", "+"): [2000]}
        return r._read_extent_contradicted()

    def test_a_lone_upstream_read_does_not_contradict(self, tmp_path):
        hook = self._hook(tmp_path, n_upstream=1)
        model = _cand("anno", 2000, 3000, source="gtf")
        assert "tss" not in hook(model, ["u0"] + [f"b{i}" for i in range(120)])

    def test_a_cluster_of_upstream_reads_does(self, tmp_path):
        hook = self._hook(tmp_path, n_upstream=4)
        model = _cand("anno", 2000, 3000, source="gtf")
        which = hook(model, [f"u{j}" for j in range(4)] + [f"b{i}" for i in range(120)])
        assert sorted(which.get("tss", [])) == ["u0", "u1", "u2", "u3"]

    def test_threshold_one_restores_the_old_behaviour(self, tmp_path):
        hook = self._hook(tmp_path, n_upstream=1, min_reads=1)
        model = _cand("anno", 2000, 3000, source="gtf")
        assert hook(model, ["u0"] + [f"b{i}" for i in range(120)]).get("tss") == ["u0"]

    def test_start_support_counts_within_tolerance(self):
        from fin.analysis.sample_nuisance import (
            TSS_CLUSTER_TOL_BP, sample_nuisance_evidence)
        starts = {f"r{i}": ("c", "+", 1000 + i) for i in range(5)}
        starts["far"] = ("c", "+", 5000)
        starts["other"] = ("c", "-", 1000)
        ev = sample_nuisance_evidence({}, {}, read_starts=starts,
                                      claimed_tss={("c", "+"): [900]})
        assert ev.start_support("c", "+", 1002) == 5          # 1000..1004 within 10
        assert ev.start_support("c", "+", 5000) == 1
        assert ev.start_support("c", "-", 1000) == 1          # strand-separated
        assert ev.start_support("c", "+", 1002, tol=1) == 3
        assert ev.start_support("zz", "+", 1) == 0
        assert TSS_CLUSTER_TOL_BP == 10


class TestUncontradictedAnnotationExemption:
    """The support gate asks a read to DISCRIMINATE a candidate. A 5'-degraded
    fragment shared by several annotated isoforms cannot, so on a real
    transcriptome the gate drops genuine transcripts. The exemption keeps its
    anti-corruption purpose by asking whether the reads CONTRADICT the
    annotation's junctions -- a jittered junction sits beside an observed one,
    a merely unobserved junction does not."""

    def _cands(self):
        from fin.candidates.dataclasses import IntronChain
        def mk(cid, introns, source):
            c = _cand(cid, 1000, 5000, source=source)
            c.intron_chain = IntronChain(introns=introns)
            return c
        return [
            mk("real_anno", ((1200, 1500), (2000, 2400)), "gtf"),   # unobserved but not contradicted
            mk("jittered_anno", ((1203, 1500), (2000, 2400)), "gtf"),  # 3 bp off an observed junction
            mk("novel_x", ((1200, 1500),), "novel"),                 # novel: never exempt
            mk("mono_anno", (), "gtf"),                              # mono: no junctions
        ]

    def test_only_uncontradicted_annotated_chains_are_exempt(self):
        from fin.scoring.m2_junction_nll import uncontradicted_gtf_columns
        observed = {"+": {(1200, 1500): 7}}          # the real junction was seen
        ex = uncontradicted_gtf_columns(self._cands(), observed, wobble_bp=6)
        assert ex == {0}                              # jittered is contradicted; novel/mono excluded

    def test_containment_only_drops_candidates_with_a_wobble_twin(self):
        """A candidate with a same-intron-count twin COULD have been decided by
        a read; only a structurally undiscriminable one is exempt."""
        from fin.candidates.dataclasses import IntronChain
        from fin.scoring.m2_junction_nll import uncontradicted_gtf_columns
        cands = self._cands()
        # nothing observed nearby -> both annotated chains are uncontradicted,
        # but they are wobble twins of each other (3 bp apart)
        ex_all = uncontradicted_gtf_columns(cands, {"+": {(9000, 9500): 3}})
        assert ex_all == {0, 1}
        assert uncontradicted_gtf_columns(
            cands, {"+": {(9000, 9500): 3}}, containment_only=True) == set()
        # a nested pair (different intron counts) has no twin -> exempt
        nested = [cands[0]]
        c = _cand("longer_anno", 1000, 5000, source="gtf")
        c.intron_chain = IntronChain(introns=((1200, 1500), (2000, 2400), (3000, 3400)))
        nested.append(c)
        assert uncontradicted_gtf_columns(
            nested, {"+": {(9000, 9500): 3}}, containment_only=True) == {0, 1}

    def test_silence_is_not_contradiction(self):
        from fin.scoring.m2_junction_nll import uncontradicted_gtf_columns
        # the reads produced a junction far away: says nothing about either model
        ex = uncontradicted_gtf_columns(self._cands(), {"+": {(9000, 9500): 3}})
        assert ex == {0, 1}

    def test_no_registry_exempts_nothing(self):
        from fin.scoring.m2_junction_nll import uncontradicted_gtf_columns
        assert uncontradicted_gtf_columns(self._cands(), None) == set()
        assert uncontradicted_gtf_columns(self._cands(), {}) == set()

    def test_exemption_is_off_by_default(self):
        from fin.pipeline.config import PIPELINE_PROFILES, PipelineConfig
        assert PipelineConfig(bam_path="/x.bam").m2_support_gate_gtf_exempt is False
        for name in PIPELINE_PROFILES:
            assert "m2_support_gate_gtf_exempt" not in PIPELINE_PROFILES[name]


class TestExemptionNeedsMass:
    def test_exempt_candidate_must_still_carry_reads(self):
        """The exemption says no read COULD have singled the candidate out, not
        that no read supports it: an EM-massless annotation stays gated."""
        from types import SimpleNamespace
        import numpy as np
        from fin.analysis.quantification import QuantResult
        from fin.candidates.dataclasses import IntronChain
        from fin.pipeline import selection as S

        def mk(cid, introns):
            c = _cand(cid, 1000, 5000, source="gtf")
            c.intron_chain = IntronChain(introns=introns)
            return c
        # two annotated isoforms of different intron counts (a nested contest),
        # neither contradicted; one carries 3 reads of mass, the other 0.5
        cands = [mk("rich", ((1200, 1500), (2000, 2400))), mk("poor", ((1200, 1500),))]
        R = np.zeros((4, 2)); R[:3, 0] = 1.0; R[3, 1] = 0.5
        qrs = [QuantResult("rich", 3.0, 1.0, 3, "gtf"), QuantResult("poor", 0.5, 1.0, 1, "gtf")]
        qo = SimpleNamespace(cand_list=cands, R=R, ties_by_read={0: [0, 1]},
                             nlls_by_read={}, quant_results=qrs,
                             kept_read_ids=[f"r{i}" for i in range(4)])
        cfg = SimpleNamespace(
            m2_support_gate=True, m2_support_gate_tie=True,
            m2_support_gate_gtf_exempt=True,
            m2_support_gate_gtf_exempt_containment_only=True,
            m2_support_gate_gtf_exempt_min_mass=2.0,
            chain_cluster_wobble_bp=6, containment_evidence_mode="off",
            junction_election_mode="off", m2_cluster_recheck=False,
            novel_junction_min_reads=0, guided_junction_min_reads=0,
            containment_cluster=False, containment_collapse=False,
            mono_resolve_post_em=False, complete_chain_rescue_reads=0,
            work_dir="/tmp")
        cs = SimpleNamespace(read_spans={})
        iv = GenomicInterval(chrom="c", start=0, end=9000, strand="+")
        res, _ = S.select_m2_interval(cfg, qo, cs, iv, lambda i: {"+": {(9500, 9600): 4}})
        kept = {q.candidate_id for q in res}
        assert "rich" in kept and "poor" not in kept


class TestExemptionIgnoresDroppedColumns:
    """A column another gate already folded away is not a rival that could have
    discriminated anything: it must neither earn an exemption nor block one."""

    def test_masked_columns_neither_earn_nor_block(self):
        from fin.candidates.dataclasses import IntronChain
        from fin.scoring.m2_junction_nll import uncontradicted_gtf_columns
        def mk(cid, introns):
            c = _cand(cid, 1000, 5000, source="gtf")
            c.intron_chain = IntronChain(introns=introns)
            return c
        survivor = mk("anno", ((1200, 1500), (2000, 2400)))
        sibling = mk("anno.ext", ((1200, 1500), (2000, 2400)))   # identical chain
        obs = {"+": {(9000, 9500): 3}}                            # silence nearby
        # with the sibling present it is a twin and blocks the exemption
        assert uncontradicted_gtf_columns([survivor, sibling], obs,
                                          containment_only=True) == set()
        # masked (already folded): the survivor is exempt, the mask earns nothing
        assert uncontradicted_gtf_columns([survivor, None], obs,
                                          containment_only=True) == {0}
        assert uncontradicted_gtf_columns([None, None], obs) == set()


class TestPolyaQcAdmission:
    """``polya_end_accept_adapter``: an ADAPTER call means krill LOCATED the
    polyA region but its LENGTH is adapter-contaminated. The 3' channel scores
    the read's genomic end and never the length, so the verdict alone decides.
    Default OFF -- the default path must behave exactly as PASS-only did."""

    def test_adapter_abstains_by_default_and_scores_when_admitted(self):
        trunc = _cand("anno_trunc", 1000, 1900, source="gtf")
        full = _cand("novel_full", 1000, 2000)
        polya = {"r": (None, "ADAPTER")}
        ends = {"r": ("c", "+", 2000)}

        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        assert apply_endpoint_likelihood(d, ties, ["r"], [trunc, full], _ctx(polya, ends)) == 0
        assert (d == 0).all()

        d = np.zeros((1, 2))
        ctx = _ctx(polya, ends); ctx["accept_adapter"] = True
        assert apply_endpoint_likelihood(d, ties, ["r"], [trunc, full], ctx) == 1
        assert d[0, 1] == 0.0 and d[0, 0] > 5.0

    def test_admission_does_not_reach_the_other_qc_classes(self):
        """NOREGION lands a median 334 bp from the true end on human dRNA and
        READ_FAILED is krill declining to process the read; neither becomes
        evidence when ADAPTER is admitted."""
        trunc = _cand("anno_trunc", 1000, 1900, source="gtf")
        full = _cand("novel_full", 1000, 2000)
        reads = ["r_noregion", "r_failed", "r_none"]
        polya = {"r_noregion": (None, "NOREGION"), "r_failed": (None, "READ_FAILED"),
                 "r_none": None}
        ends = {r: ("c", "+", 2000) for r in reads}
        d = np.zeros((3, 2)); ties = {i: [0, 1] for i in range(3)}
        ctx = _ctx(polya, ends); ctx["accept_adapter"] = True
        assert apply_endpoint_likelihood(d, ties, reads, [trunc, full], ctx) == 0
        assert (d == 0).all()

    def test_a_pass_verdict_without_a_usable_length_is_still_not_evidence(self):
        """The PASS branch stays byte-identical to ``polya_state``'s PASS, so
        admitting ADAPTER cannot widen PASS by a side effect."""
        trunc = _cand("anno_trunc", 1000, 1900, source="gtf")
        full = _cand("novel_full", 1000, 2000)
        d = np.zeros((2, 2)); ties = {0: [0, 1], 1: [0, 1]}
        ctx = _ctx({"r_zero": (0.0, "PASS"), "r_nolen": (None, "PASS")},
                   {"r_zero": ("c", "+", 2000), "r_nolen": ("c", "+", 2000)})
        ctx["accept_adapter"] = True
        assert apply_endpoint_likelihood(d, ties, ["r_zero", "r_nolen"],
                                         [trunc, full], ctx) == 0
        assert (d == 0).all()


class TestQcConditionalChannel:
    """With per-class kernels the channel scores every read krill gave a verdict
    for, under its OWN class's kernel, and abstains for a class whose kernel the
    sample could not fit."""

    def _cands(self):
        return [_cand("anno_trunc", 1000, 1900, source="gtf"), _cand("novel_full", 1000, 2000)]

    def _ctx(self, polya, ends, by_qc):
        return {"endpoint": None, "polya": polya, "ends": ends,
                "qc_conditional": True, "endpoint_by_qc": by_qc}

    def test_each_class_scores_under_its_own_kernel(self):
        good = EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.04, n_fit=500)
        noisy = EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.40, n_fit=500)
        reads = ["r_pass", "r_failed"]
        d = np.zeros((2, 2)); ties = {0: [0, 1], 1: [0, 1]}
        ctx = self._ctx({"r_pass": (40.0, "PASS"), "r_failed": (None, "READ_FAILED")},
                        {r: ("c", "+", 2000) for r in reads},
                        {"PASS": good, "READ_FAILED": noisy})
        assert apply_endpoint_likelihood(d, ties, reads, self._cands(), ctx) == 2
        # both penalise the truncated model, the noisy class by less
        assert d[0, 0] > d[1, 0] > 0
        assert d[0, 1] == 0.0 and d[1, 1] == 0.0

    def test_a_class_without_a_usable_kernel_abstains(self):
        good = EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.04, n_fit=500)
        unfit = EndpointModel(b_narrow=4407.0, b_broad=500.0, eta=0.81, n_fit=1702)
        assert not unfit.available
        reads = ["r_noregion", "r_unlisted"]
        d = np.zeros((2, 2)); ties = {0: [0, 1], 1: [0, 1]}
        ctx = self._ctx({"r_noregion": (None, "NOREGION"), "r_unlisted": (None, "OTHER")},
                        {r: ("c", "+", 2000) for r in reads},
                        {"PASS": good, "NOREGION": unfit})
        assert apply_endpoint_likelihood(d, ties, reads, self._cands(), ctx) == 0
        assert (d == 0).all()

    def test_class_kernels_carry_the_channel_when_the_pooled_fit_is_absent(self):
        """The context passes `endpoint: None` in conditional mode whenever the
        pooled fit is unavailable; the class kernels must still switch the 3'
        channel on, or the runner and the nuisance model disagree."""
        good = EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.04, n_fit=500)
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = self._ctx({"r": (40.0, "PASS")}, {"r": ("c", "+", 2000)}, {"PASS": good})
        assert ctx["endpoint"] is None
        assert apply_endpoint_likelihood(d, ties, ["r"], self._cands(), ctx) == 1
        assert d[0, 0] > 5.0 and d[0, 1] == 0.0

    def test_no_class_qualifying_is_not_a_licence_to_use_the_pooled_kernel(self):
        """An empty class map in conditional mode means no class fitted, which is
        abstention -- not a fallback to the pooled model."""
        pooled = EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.04, n_fit=500)
        d = np.zeros((1, 2)); ties = {0: [0, 1]}
        ctx = self._ctx({"r": (40.0, "PASS")}, {"r": ("c", "+", 2000)}, {})
        ctx["endpoint"] = pooled
        assert apply_endpoint_likelihood(d, ties, ["r"], self._cands(), ctx) == 0
        assert (d == 0).all()


class TestDiscoveryIsNotWidenedByTheConditionalMode:
    """`--endpoint-qc-conditional` changes the ASSIGNMENT channel only. The
    discovery contradiction hook tests `narrow_posterior(d) < 0.5`, an ABSOLUTE
    threshold: a class fitted as less reliable has a lower narrow posterior at
    every distance, so it would raise MORE contradictions and be protected by
    FEWER annotated alternatives. Less confidence must not mean more aggression,
    so this hook keeps the pooled kernel and the verdict-name rule."""

    def _hook(self, tmp_path, *, qc_conditional):
        """200 PASS reads AT the annotated TES define the narrow kernel; 30 PASS
        reads 400 bp past it are genuine contradictions in either mode; and 400
        NOREGION reads are shaped so their own class kernel is AVAILABLE and
        high-eta -- otherwise a wrong implementation that selected class kernels
        here would abstain for an unrelated reason and the comparison would pass
        vacuously."""
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner
        tmp_path.mkdir(parents=True, exist_ok=True)
        cfg = PipelineConfig(bam_path=str(_tiny_bam(tmp_path / "r.bam")),
                             work_dir=str(tmp_path / "w"), signal_path="/fake.blow5",
                             assignment_channels="ends", quant_mode="m2_em",
                             endpoint_qc_conditional=qc_conditional)
        r = PipelineRunner(cfg)
        polya, ends = {}, {}
        for i in range(200):
            polya[f"p{i}"] = [40.0, "PASS"]
            ends[f"p{i}"] = ["c", "+", 3000 + (i % 3) - 1, 2000]
        for i in range(30):
            polya[f"q{i}"] = [40.0, "PASS"]
            ends[f"q{i}"] = ["c", "+", 3400 + (i % 3) - 1, 2000]
        for i in range(400):
            polya[f"n{i}"] = [-1.0, "NOREGION"]
            e = (i % 3) - 1 if i % 20 < 9 else ((i * 53) % 401) - 200
            ends[f"n{i}"] = ["c", "+", 3000 + e, 2000]
        (tmp_path / "w").mkdir(parents=True, exist_ok=True)
        r._evidence_path().write_text(json.dumps(
            {"polya": polya, "ends": ends, "version": "sample-evidence-v3",
             "identity": r._evidence_identity(), "polya_status": "ok",
             "n_ends": len(ends), "n_polya": len(polya)}))
        r._annotated_tes = lambda: {("c", "+"): [3000]}
        r._annotated_tss = lambda: {("c", "+"): [2000]}
        return r

    def _reads(self):
        return ([f"p{i}" for i in range(200)] + [f"q{i}" for i in range(30)]
                + [f"n{i}" for i in range(400)])

    def test_the_noregion_class_really_would_have_scored(self, tmp_path):
        """The precondition the comparison below depends on: fitted per class,
        NOREGION is available and low-confidence, so an implementation that let
        class kernels into this hook WOULD produce extra contradictions."""
        from fin.analysis.sample_nuisance import ETA_MAX, fit_from_evidence

        r = self._hook(tmp_path, qc_conditional=True)
        m = fit_from_evidence(r._nuisance_evidence())
        k = m.endpoint_for("NOREGION")
        assert k is not None and k.available and k.eta > ETA_MAX
        # and such a kernel would indeed call the scattered reads contradictions
        assert k.narrow_posterior(150) < 0.5

    def test_the_hook_gives_the_same_verdicts_with_the_mode_on(self, tmp_path):
        model = _cand("anno", 2000, 3000, source="gtf")
        reads = self._reads()
        off = self._hook(tmp_path / "off", qc_conditional=False)
        on = self._hook(tmp_path / "on", qc_conditional=True)
        a = off._read_extent_contradicted()(model, reads)
        b = on._read_extent_contradicted()(model, reads)
        # non-vacuous: the 30 PASS reads past the annotated end DO contradict it,
        # and no NOREGION read is admitted in either mode
        assert a.get("tes") == sorted(f"q{i}" for i in range(30))
        assert a == b

    def test_the_same_holds_when_an_alternative_annotation_protects_reads(self, tmp_path):
        model = _cand("anno", 2000, 3000, source="gtf")
        alt = _cand("alt", 2000, 3400, source="gtf")
        reads = self._reads()
        off = self._hook(tmp_path / "off2", qc_conditional=False)
        on = self._hook(tmp_path / "on2", qc_conditional=True)
        a = off._read_extent_contradicted()(model, reads, [alt])
        b = on._read_extent_contradicted()(model, reads, [alt])
        # the alternative explains exactly the reads that contradicted before
        assert "tes" not in a
        assert a == b

    def test_adapter_admission_does_reach_the_hook(self, tmp_path):
        """The name-based widening is a different thing: one pooled kernel over a
        wider population, no reliability inversion, and it is meant to apply to
        both consumers."""
        import json
        from fin.pipeline.config import PipelineConfig
        from fin.pipeline.runner import PipelineRunner

        def build(accept):
            d = tmp_path / f"ad{int(accept)}"
            d.mkdir(parents=True, exist_ok=True)
            cfg = PipelineConfig(bam_path=str(_tiny_bam(d / "r.bam")),
                                 work_dir=str(d / "w"), signal_path="/fake.blow5",
                                 assignment_channels="ends", quant_mode="m2_em",
                                 polya_end_accept_adapter=accept)
            r = PipelineRunner(cfg)
            polya, ends = {}, {}
            for i in range(200):
                polya[f"p{i}"] = [40.0, "PASS"]
                ends[f"p{i}"] = ["c", "+", 3000 + (i % 3) - 1, 2000]
            for i in range(60):
                polya[f"a{i}"] = [None, "ADAPTER"]
                ends[f"a{i}"] = ["c", "+", 3000 + 60 + (i % 3), 2000]
            (d / "w").mkdir(parents=True, exist_ok=True)
            r._evidence_path().write_text(json.dumps(
                {"polya": polya, "ends": ends, "version": "sample-evidence-v3",
                 "identity": r._evidence_identity(), "polya_status": "ok",
                 "n_ends": len(ends), "n_polya": len(polya)}))
            r._annotated_tes = lambda: {("c", "+"): [3000]}
            r._annotated_tss = lambda: {("c", "+"): [2000]}
            return r._read_extent_contradicted()

        model = _cand("anno", 2000, 3000, source="gtf")
        reads = [f"p{i}" for i in range(200)] + [f"a{i}" for i in range(60)]
        assert not build(False)(model, reads).get("tes")
        assert build(True)(model, reads).get("tes")
