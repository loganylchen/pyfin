"""Tests for finalized-model read-supported junction correction."""
from __future__ import annotations

from collections import Counter
from dataclasses import replace

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
    fulllen_frac: float = -1.0,
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
        fulllen_frac=fulllen_frac,
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


class TestSnapCollisionChainMerge:
    """`snap_merge_chain`: reconcile the models the SNAP itself gave one
    corrected intron chain. Discovery emits two wobble variants of one chain; the
    snap corrects the minority junction onto the consensus; the two then carry
    the same chain, but the default merge key is the full exon tuple, so any
    difference in the terminal coordinates keeps them apart and the chain-based
    scorer credits one and charges the other.

    The rule is deliberately narrow, because one corrected intron chain is NOT
    one transcript -- APA and alternative-TSS isoforms share a chain -- so each
    of the conditions below has a test that it actually blocks the merge."""

    # 0-based half-open introns. (200,400) is the read consensus for the first
    # junction, (200,403) the 3 bp wobble the snap will correct onto it.
    OBS = {("chr1", "+"): Counter({(200, 400): 20, (200, 403): 1,
                                   (500, 700): 20})}

    def _pair(self, *, minority_start=125, family="F1", minority_family=None,
              minority_source="novel"):
        major = _qr("major", [(100, 200), (400, 500), (700, 800)],
                    abundance=10.0, assigned=("r1", "r2", "r3"), family_id=family,
                    fulllen_frac=0.9)
        # same chain after snapping, 25 bp shorter at the 5', one junction wobbled
        minor = _qr("minor", [(minority_start, 200), (403, 500), (700, 800)],
                    abundance=2.0, assigned=("r4",),
                    family_id=minority_family or family,
                    source=minority_source, fulllen_frac=0.4)
        return {"major": major, "minor": minor}

    def _run(self, results, **kw):
        return snap_quant_results(
            results, self.OBS, tolerance=6, min_support=2, min_ratio=2.0,
            return_redirects=True, **kw)

    def test_the_snap_makes_them_identical_but_the_default_keeps_them_apart(self):
        merged, snapped, n, redirects = self._run(self._pair())
        assert snapped == 1                       # the wobble WAS corrected
        assert set(merged) == {"major", "minor"}  # ... and they stay separate
        assert n == 0 and redirects == {}
        assert merged["minor"].exons[1] == (400, 500)   # same chain now

    def test_enabling_the_merge_reconciles_them(self):
        merged, _snapped, n, redirects = self._run(self._pair(), chain_merge=True)
        assert set(merged) == {"major"}
        assert n == 1 and redirects == {"minor": "major"}
        got = merged["major"]
        assert got.abundance == 12.0
        assert set(got.assigned_read_ids) == {"r1", "r2", "r3", "r4"}
        # the survivor keeps its OWN geometry: no new endpoint is invented here
        assert got.exons == ((100, 200), (400, 500), (700, 800))
        assert got.start == 100 and got.end == 800
        # a fraction measured against other ends over a grown read set is stale:
        # the inputs carried real scores (0.9 / 0.4) and the survivor must not
        # report either of them, nor a weighted average of the two
        assert got.fulllen_frac == -1.0

    def test_a_pre_existing_duplicate_the_snap_did_not_create_is_left_alone(self):
        """The provenance condition. Two models already sharing a chain before
        the snap are not this mechanism, and merging them would be an ordinary
        same-chain collapse -- which is exactly what must not happen silently."""
        a = _qr("a", [(100, 200), (400, 500), (700, 800)], family_id="F1")
        b = _qr("b", [(125, 200), (400, 500), (700, 800)], family_id="F1")
        merged, snapped, n, redirects = self._run({"a": a, "b": b}, chain_merge=True)
        assert snapped == 0
        assert set(merged) == {"a", "b"} and n == 0 and redirects == {}

    def test_a_different_three_prime_end_blocks_the_merge(self):
        """An alternative last exon is a different transcript end, not a wobble."""
        p = self._pair()
        p["minor"] = _qr("minor", [(125, 200), (403, 500), (700, 900)],
                         abundance=2.0, family_id="F1")
        merged, _s, n, _r = self._run(p, chain_merge=True)
        assert set(merged) == {"major", "minor"} and n == 0

    def test_a_five_prime_spread_beyond_the_window_blocks_the_merge(self):
        merged, _s, n, _r = self._run(self._pair(minority_start=40),
                                      chain_merge=True)
        assert set(merged) == {"major", "minor"} and n == 0
        merged, _s, n, _r = self._run(self._pair(minority_start=40),
                                      chain_merge=True, chain_merge_window_bp=100)
        assert set(merged) == {"major"} and n == 1

    def test_a_different_discovery_family_blocks_the_merge(self):
        merged, _s, n, _r = self._run(self._pair(minority_family="F2"),
                                      chain_merge=True)
        assert set(merged) == {"major", "minor"} and n == 0

    def test_a_missing_family_blocks_the_merge(self):
        p = self._pair()
        p["minor"] = _qr("minor", [(125, 200), (403, 500), (700, 800)],
                         abundance=2.0, family_id=None)
        merged, _s, n, _r = self._run(p, chain_merge=True)
        assert set(merged) == {"major", "minor"} and n == 0

    def test_an_annotated_model_is_never_absorbed(self):
        """The annotated model carries the CONSENSUS chain, so it is not kept out
        by failing to snap: two novel models around it do snap and would merge,
        and the annotated one still must not be absorbed."""
        anno = _qr("anno", [(100, 200), (400, 500), (700, 800)],
                   source="gtf", abundance=20.0, assigned=("r0",), family_id="F1")
        novel_a = _qr("nov_a", [(110, 200), (400, 500), (700, 800)],
                      abundance=10.0, assigned=("r1", "r2"), family_id="F1")
        novel_b = _qr("nov_b", [(125, 200), (403, 500), (700, 800)],
                      abundance=2.0, assigned=("r3",), family_id="F1")
        merged, snapped, n, redirects = self._run(
            {"anno": anno, "nov_a": novel_a, "nov_b": novel_b}, chain_merge=True)
        assert snapped == 1                       # the novel wobble WAS corrected
        assert n == 1 and redirects == {"nov_b": "nov_a"}
        assert set(merged) == {"anno", "nov_a"}   # the annotation kept its identity
        assert merged["anno"].exons == ((100, 200), (400, 500), (700, 800))

    def test_the_grouping_is_bounded_diameter_not_transitive(self):
        """5' ends at 100/125/150 must not chain into one group through two
        25 bp links: the outer pair is 50 bp apart."""
        base = [(100, 200), (400, 500), (700, 800)]
        rs = {"a": _qr("a", base, abundance=9.0, family_id="F1"),
              "b": _qr("b", [(125, 200), (403, 500), (700, 800)],
                       abundance=5.0, family_id="F1"),
              "c": _qr("c", [(150, 200), (403, 500), (700, 800)],
                       abundance=1.0, family_id="F1")}
        merged, _s, _n, redirects = self._run(rs, chain_merge=True)
        assert len(merged) == 2
        assert redirects == {"b": "a"}
        assert set(merged) == {"a", "c"}

    def test_redirects_compose_through_both_passes(self):
        """A model absorbed by the exon-tuple pass and then by this one must
        redirect to the FINAL survivor, not to an id that no longer exists."""
        base = [(100, 200), (400, 500), (700, 800)]
        rs = {"major": _qr("major", base, abundance=10.0,
                           assigned=("r1", "r2", "r3"), family_id="F1"),
              "dupe": _qr("dupe", [(125, 200), (403, 500), (700, 800)],
                          abundance=4.0, assigned=("r4",), family_id="F1"),
              "twin": _qr("twin", [(125, 200), (403, 500), (700, 800)],
                          abundance=1.0, assigned=("r5",), family_id="F1")}
        merged, _s, n, redirects = self._run(rs, chain_merge=True)
        assert set(merged) == {"major"}
        # `twin` was absorbed by `dupe` in the exon-tuple pass and `dupe` by
        # `major` here: its redirect must name the FINAL survivor
        assert redirects == {"dupe": "major", "twin": "major"}
        assert set(merged["major"].assigned_read_ids) == {"r1", "r2", "r3",
                                                          "r4", "r5"}
        assert n == 2


class TestSnapMergeChainConfigGuards:
    """The mode is only coherent with the snap it reconciles, and it must not be
    able to undo EndpointRefine's endpoint states."""

    def _cfg(self, tmp_path, **kw):
        from fin.pipeline.config import PipelineConfig
        bam = tmp_path / "r.bam"; bam.write_bytes(b"")
        base = dict(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                    snap_merge_chain=True, junction_snap=True)
        base.update(kw)
        cfg = PipelineConfig(**base)
        cfg.validate()
        return cfg

    def test_it_requires_junction_snap(self, tmp_path):
        import pytest
        with pytest.raises(ValueError, match="junction_snap is off"):
            self._cfg(tmp_path, junction_snap=False)

    def test_it_refuses_to_run_beside_endpoint_refine(self, tmp_path):
        import pytest
        with pytest.raises(ValueError, match="endpoint_refine is on"):
            self._cfg(tmp_path, endpoint_refine=True)

    def test_it_is_off_by_default(self, tmp_path):
        from fin.pipeline.config import PipelineConfig
        bam = tmp_path / "r.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(tmp_path / "w"))
        cfg.validate()
        assert cfg.snap_merge_chain is False
        assert cfg.snap_merge_chain_window_bp == 25


class TestSnapCollisionProvenanceIsPerPair:
    """The snap must have created the relationship between EVERY pair being
    collapsed. A wobble model that happens to land in the group must not supply
    the permission to collapse two models that already shared their chain."""

    OBS = {("chr1", "+"): Counter({(200, 400): 20, (200, 403): 1,
                                   (200, 406): 1, (500, 700): 20})}

    def _run(self, results, **kw):
        return snap_quant_results(
            results, self.OBS, tolerance=6, min_support=2, min_ratio=2.0,
            return_redirects=True, chain_merge=True, **kw)

    def test_a_pre_existing_endpoint_pair_is_not_collapsed_by_a_third_wobble(self):
        cons = [(100, 200), (400, 500), (700, 800)]
        rs = {
            # a and b already shared the consensus chain before any snapping
            "a": _qr("a", cons, abundance=10.0, assigned=("r1",), family_id="F1"),
            "b": _qr("b", [(125, 200), (400, 500), (700, 800)],
                     abundance=6.0, assigned=("r2",), family_id="F1"),
            # c is a wobble that snaps onto a's chain and is absorbed by it
            "c": _qr("c", [(100, 200), (403, 500), (700, 800)],
                     abundance=1.0, assigned=("r3",), family_id="F1"),
        }
        merged, snapped, n, redirects = self._run(rs)
        assert snapped == 1
        # the ONE merge is the first pass absorbing c into a; b is untouched
        assert redirects == {"c": "a"} and n == 1
        assert set(merged) == {"a", "b"}

    def test_two_wobbles_of_distinct_originals_still_merge(self):
        """The control for the rule above: when the two survivors' pre-snap
        chains are disjoint, the collision really is snap-created."""
        rs = {
            "a": _qr("a", [(100, 200), (400, 500), (700, 800)],
                     abundance=10.0, assigned=("r1",), family_id="F1"),
            "b": _qr("b", [(125, 200), (403, 500), (700, 800)],
                     abundance=6.0, assigned=("r2",), family_id="F1"),
        }
        merged, _s, n, redirects = self._run(rs)
        assert set(merged) == {"a"} and n == 1 and redirects == {"b": "a"}

    def test_a_first_pass_merge_across_two_families_cannot_launder_them(self):
        """`_merge_group` labels its survivor with the lexicographically first
        non-null family, so a first-pass group spanning F1 and F2 would present
        itself as F1 unless the ORIGINAL families are what is checked."""
        rs = {
            "a": _qr("a", [(100, 200), (400, 500), (700, 800)],
                     abundance=10.0, assigned=("r1",), family_id="F1"),
            # same exon tuple as `a` after snapping -> absorbed by the first pass,
            # bringing a foreign family with it
            "foreign": _qr("foreign", [(100, 200), (403, 500), (700, 800)],
                           abundance=3.0, assigned=("r2",), family_id="F2"),
            # a THIRD pre-snap chain, so the origin sets {P,Q} and {R} are
            # disjoint and the provenance guard cannot be what blocks this --
            # the family check has to be doing the work
            "b": _qr("b", [(125, 200), (406, 500), (700, 800)],
                     abundance=2.0, assigned=("r3",), family_id="F1"),
        }
        merged, _s, n, redirects = self._run(rs)
        # the ONE merge is the first pass absorbing `foreign`; b is untouched
        assert redirects == {"foreign": "a"} and n == 1
        assert set(merged) == {"a", "b"}

    def test_a_member_without_a_family_cannot_borrow_a_siblings(self):
        rs = {
            "a": _qr("a", [(100, 200), (400, 500), (700, 800)],
                     abundance=10.0, assigned=("r1",), family_id="F1"),
            "nofam": _qr("nofam", [(100, 200), (403, 500), (700, 800)],
                         abundance=3.0, assigned=("r2",), family_id=None),
            # a third pre-snap chain, so provenance cannot be the blocker here
            "b": _qr("b", [(125, 200), (406, 500), (700, 800)],
                     abundance=2.0, assigned=("r3",), family_id="F1"),
        }
        merged, _s, n, redirects = self._run(rs)
        assert redirects == {"nofam": "a"} and n == 1
        assert set(merged) == {"a", "b"}

    def test_minus_strand_uses_the_strand_aware_ends(self):
        """On `-` the 3' end is `start` and the 5' end is `end`: a pair sharing a
        3' end and differing 25 bp at the 5' must merge, and a pair differing at
        the 3' end must not."""
        def qr(cid, exons, **kw):
            r = _qr(cid, exons, **kw)
            return replace(r, strand="-")
        obs = {("chr1", "-"): Counter({(200, 400): 20, (200, 403): 1,
                                       (500, 700): 20})}
        same3 = {"a": qr("a", [(100, 200), (400, 500), (700, 800)],
                         abundance=10.0, assigned=("r1",), family_id="F1"),
                 "b": qr("b", [(100, 200), (403, 500), (700, 775)],
                         abundance=2.0, assigned=("r2",), family_id="F1")}
        merged, _s, n, redirects = snap_quant_results(
            same3, obs, tolerance=6, min_support=2, min_ratio=2.0,
            return_redirects=True, chain_merge=True)
        assert set(merged) == {"a"} and n == 1 and redirects == {"b": "a"}

        diff3 = {"a": qr("a", [(100, 200), (400, 500), (700, 800)],
                         abundance=10.0, assigned=("r1",), family_id="F1"),
                 "b": qr("b", [(125, 200), (403, 500), (700, 800)],
                         abundance=2.0, assigned=("r2",), family_id="F1")}
        merged, _s, n, _r = snap_quant_results(
            diff3, obs, tolerance=6, min_support=2, min_ratio=2.0,
            return_redirects=True, chain_merge=True)
        assert set(merged) == {"a", "b"} and n == 0

    def test_a_negative_window_is_rejected(self):
        import pytest
        with pytest.raises(ValueError, match="must be >= 0"):
            self._run({"a": _qr("a", [(100, 200), (400, 500), (700, 800)],
                                family_id="F1"),
                       "b": _qr("b", [(125, 200), (403, 500), (700, 800)],
                                family_id="F1")},
                      chain_merge_window_bp=-1)


class TestSnapMergeChainAdditionalCoverage:
    OBS = {("chr1", "+"): Counter({(200, 400): 30, (200, 403): 1, (200, 406): 1,
                                   (200, 409): 1, (500, 700): 30})}

    def _run(self, results, **kw):
        return snap_quant_results(
            results, self.OBS, tolerance=10, min_support=2, min_ratio=2.0,
            return_redirects=True, chain_merge=True, **kw)

    def test_three_survivors_with_disjoint_origins_merge_into_one(self):
        rs = {"a": _qr("a", [(100, 200), (400, 500), (700, 800)],
                       abundance=9.0, assigned=("r1", "r2"), family_id="F1"),
              "b": _qr("b", [(110, 200), (403, 500), (700, 800)],
                       abundance=5.0, assigned=("r3",), family_id="F1"),
              "c": _qr("c", [(120, 200), (406, 500), (700, 800)],
                       abundance=1.0, assigned=("r4",), family_id="F1")}
        merged, _s, n, redirects = self._run(rs)
        assert set(merged) == {"a"} and n == 2
        assert redirects == {"b": "a", "c": "a"}
        assert set(merged["a"].assigned_read_ids) == {"r1", "r2", "r3", "r4"}

    def test_a_single_exon_model_is_never_bucketed(self):
        mono = _qr("mono", [(100, 800)], abundance=4.0, family_id="F1")
        rs = {"a": _qr("a", [(100, 200), (400, 500), (700, 800)],
                       abundance=9.0, family_id="F1"),
              "b": _qr("b", [(125, 200), (403, 500), (700, 800)],
                       abundance=5.0, family_id="F1"),
              "mono": mono}
        merged, _s, _n, _r = self._run(rs)
        assert "mono" in merged

    def test_an_empty_origin_set_blocks_the_merge(self):
        """`not all(sets)`: no recorded pre-snap chain is UNKNOWN provenance, and
        unknown provenance must not be read as disjoint provenance. Driven
        against the helper directly, because through the public path every
        eligible survivor has its origins recorded."""
        from fin.pipeline.junction_snap import _merge_snap_collisions
        a = _qr("a", [(100, 200), (400, 500), (700, 800)],
                abundance=9.0, assigned=("r1",), family_id="F1")
        b = _qr("b", [(125, 200), (400, 500), (700, 800)],
                abundance=5.0, assigned=("r2",), family_id="F1")
        chain_a = ((200, 403), (500, 700))
        common = dict(members_of={"a": ["a"], "b": ["b"]},
                      origin_families={"a": frozenset({"F1"}),
                                       "b": frozenset({"F1"})},
                      snapped_ids={"a"})
        out, extra, _r = _merge_snap_collisions(
            {"a": a, "b": b}, redirects={},
            origin_chains={"a": frozenset({chain_a}),
                           "b": frozenset({((200, 406), (500, 700))})},
            window=25, **common)
        assert set(out) == {"a"} and extra == 1        # both origins recorded
        out, extra, _r = _merge_snap_collisions(
            {"a": a, "b": b}, redirects={},
            origin_chains={"a": frozenset({chain_a}), "b": frozenset()},
            window=25, **common)
        assert set(out) == {"a", "b"} and extra == 0   # b's provenance unknown


class TestSnapMergeChainRefitGuard:
    """Merging changes every read's responsibility for the merged identity, so
    the mode is only correct when the post-selection refit will recompute
    confidence, max_R and the hard assignments."""

    def _cfg(self, tmp_path, **kw):
        from fin.pipeline.config import PipelineConfig
        tmp_path.mkdir(parents=True, exist_ok=True)
        bam = tmp_path / "r.bam"; bam.write_bytes(b"")
        base = dict(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                    snap_merge_chain=True, junction_snap=True,
                    post_selection_refit=True, quant_mode="m2_em",
                    abundance_feedback=False)
        base.update(kw)
        cfg = PipelineConfig(**base)
        cfg.validate()
        return cfg

    def test_the_refit_conditions_are_each_enforced(self, tmp_path):
        import pytest
        assert self._cfg(tmp_path)                       # the allowed combination
        for i, kw in enumerate((dict(post_selection_refit=False),
                                dict(quant_mode="argmax"),
                                dict(abundance_feedback=True))):
            with pytest.raises(ValueError, match="abundance refit is not in effect"):
                self._cfg(tmp_path / f"c{i}", **kw)

    def test_a_negative_window_is_rejected_by_config_not_only_by_data(self, tmp_path):
        import pytest
        with pytest.raises(ValueError, match="window_bp must be >= 0"):
            self._cfg(tmp_path / "neg", snap_merge_chain_window_bp=-1)


def test_a_negative_window_is_rejected_even_with_no_models():
    """The boundary check must not depend on a bucket having two members."""
    import pytest
    with pytest.raises(ValueError, match="must be >= 0"):
        snap_quant_results({}, {}, tolerance=6, min_support=2, min_ratio=2.0,
                           chain_merge=True, chain_merge_window_bp=-1)


class TestSnapMergeChainProfilePromotion:
    """The `real-drna` profile now carries the mode. A setting the PROFILE
    supplied must step aside with a warning when the rest of the operating point
    cannot support it -- the same contract as folded_tss_recovery_mode -- while
    one the caller asked for explicitly stays an error."""

    def _cfg(self, tmp_path, explicit=(), **kw):
        """A config as the CLI builds one: caller values overlaid by the profile,
        except the fields the caller passed explicitly (`explicit`), which the
        CLI also records in `profile_overrides`."""
        from fin.pipeline.config import PipelineConfig, resolve_profile_values
        tmp_path.mkdir(parents=True, exist_ok=True)
        bam = tmp_path / "r.bam"; bam.write_bytes(b"")
        values = dict(snap_merge_chain=False, junction_snap=False,
                      post_selection_refit=False, quant_mode="m2_em",
                      abundance_feedback=False, assignment_channels="none",
                      folded_tss_recovery_mode="off")
        values.update(kw)
        resolved = resolve_profile_values("real-drna", values,
                                          explicit_fields=explicit)
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(tmp_path / "w"),
                             profile="real-drna",
                             profile_overrides=tuple(sorted(explicit)), **resolved)
        cfg.validate()
        return cfg

    def test_the_profile_turns_it_on(self, tmp_path):
        assert self._cfg(tmp_path).snap_merge_chain is True

    def test_every_incompatibility_makes_the_profile_value_step_aside(self, tmp_path):
        """All three reasons, not just one, and each must WARN rather than raise
        when the value came from the profile."""
        import logging
        cases = (dict(quant_mode="argmax"),
                 dict(abundance_feedback=True),
                 dict(endpoint_refine=True),
                 # junction_snap IS a profile field, so turning it off is an
                 # explicit override -- but snap_merge_chain itself still came
                 # from the profile, so it must step aside rather than raise
                 dict(junction_snap=False, explicit=("junction_snap",)))
        for i, kw in enumerate(cases):
            cfg = self._cfg(tmp_path / f"step{i}", **kw)
            assert cfg.snap_merge_chain is False, kw

    def test_it_steps_aside_when_the_profile_cannot_support_it(self, tmp_path):
        # quant_mode is not a profile field, so a caller asking for argmax keeps
        # it; the profile's refit then never takes effect and the
        # profile-implied merge must disable itself rather than raise
        cfg = self._cfg(tmp_path / "argmax", quant_mode="argmax")
        assert cfg.quant_mode == "argmax"
        assert cfg.snap_merge_chain is False

    def test_an_explicit_request_is_still_an_error(self, tmp_path):
        import pytest
        from fin.pipeline.config import PipelineConfig
        d = tmp_path / "explicit"; d.mkdir(parents=True, exist_ok=True)
        bam = d / "r.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam), work_dir=str(d / "w"),
                             profile="real-drna",
                             profile_overrides=("snap_merge_chain",),
                             snap_merge_chain=True, junction_snap=False)
        with pytest.raises(ValueError, match="junction_snap is off"):
            cfg.validate()
