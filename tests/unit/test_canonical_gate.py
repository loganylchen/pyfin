"""Tests for the production canonical-motif gate (Stage B).

The gate drops NOVEL multi-exon candidates whose intron chain is not
all-canonical. GTF-passthrough (source="gtf") and fusion (source="fusion")
candidates are EXEMPT; single-exon (mono) candidates trivially pass.
"""

from __future__ import annotations

from fin.candidates.dataclasses import (
    CandidateSet,
    IntronChain,
    TranscriptCandidate,
)
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner

# Genome layout (0-based):
#   0-4   "AAAAA"
#   5-6   "GT"      <- canonical donor
#   7-12  "CCCCCC"
#   13-14 "AG"      <- canonical acceptor
#   15-19 "AAAAA"
# Canonical + strand intron = (5, 15): genome[5:7]="GT", genome[13:15]="AG".
GENOME = "AAAAA" + "GT" + "CCCCCC" + "AG" + "AAAAA"
CANON_INTRON = (5, 15)
# Non-canonical: donor genome[6:8]="TC".
NONCANON_INTRON = (6, 15)


def _cand(cid, introns, source, strand="+"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=15,
        sequence="ACGT",
        source=source,
        supporting_read_ids={f"r_{cid}"},
        chrom="chr1",
        strand=strand,
        start=0,
        end=20,
    )


def _make_set(cands):
    return CandidateSet(
        interval=GenomicInterval(chrom="chr1", start=0, end=20, strand="+"),
        candidates=list(cands),
        read_ids={c.candidate_id for c in cands},
    )


def _runner(**cfg_kwargs):
    cfg = PipelineConfig(bam_path="/fake/reads.bam", **cfg_kwargs)
    return PipelineRunner(cfg)


class TestCanonicalGate:
    def test_drops_noncanonical_novel_multiexon(self):
        cs = _make_set([
            _cand("good", [CANON_INTRON], "novel"),
            _cand("bad", [NONCANON_INTRON], "novel"),
        ])
        _runner()._apply_canonical_gate(cs, GENOME)
        ids = cs.candidate_ids()
        assert "good" in ids
        assert "bad" not in ids

    def test_keeps_canonical_novel_multiexon(self):
        cs = _make_set([_cand("good", [CANON_INTRON], "novel")])
        _runner()._apply_canonical_gate(cs, GENOME)
        assert cs.candidate_ids() == ["good"]

    def test_exempts_gtf_passthrough(self):
        # A non-canonical chain survives if it came from the GTF.
        cs = _make_set([_cand("gtf_tx", [NONCANON_INTRON], "gtf")])
        _runner()._apply_canonical_gate(cs, GENOME)
        assert cs.candidate_ids() == ["gtf_tx"]

    def test_exempts_fusion(self):
        cs = _make_set([_cand("fus", [NONCANON_INTRON], "fusion")])
        _runner()._apply_canonical_gate(cs, GENOME)
        assert cs.candidate_ids() == ["fus"]

    def test_mono_novel_trivially_passes(self):
        # Single-exon candidate has no internal junction -> always kept.
        cs = _make_set([_cand("mono", [], "novel")])
        _runner()._apply_canonical_gate(cs, GENOME)
        assert cs.candidate_ids() == ["mono"]

    def test_empty_genome_skips_gate(self):
        # No genome sequence -> cannot evaluate motifs -> keep everything.
        cs = _make_set([_cand("bad", [NONCANON_INTRON], "novel")])
        _runner()._apply_canonical_gate(cs, "")
        assert cs.candidate_ids() == ["bad"]

    def test_read_ids_untouched(self):
        cs = _make_set([
            _cand("good", [CANON_INTRON], "novel"),
            _cand("bad", [NONCANON_INTRON], "novel"),
        ])
        before = set(cs.read_ids)
        _runner()._apply_canonical_gate(cs, GENOME)
        assert cs.read_ids == before

    def test_extended_motifs_keep_gc_ag(self):
        # GC-AG intron: donor genome[s:s+2]="GC". Default motifs include GC-AG.
        genome = "AAAAA" + "GC" + "CCCCCC" + "AG" + "AAAAA"
        cs = _make_set([_cand("gc", [(5, 15)], "novel")])
        _runner()._apply_canonical_gate(cs, genome)
        assert cs.candidate_ids() == ["gc"]

    def test_strict_motifs_drop_gc_ag(self):
        # With strict GT-AG only, a GC-AG novel intron is dropped.
        genome = "AAAAA" + "GC" + "CCCCCC" + "AG" + "AAAAA"
        cs = _make_set([_cand("gc", [(5, 15)], "novel")])
        _runner(canonical_motifs=("GT-AG",))._apply_canonical_gate(cs, genome)
        assert cs.candidate_ids() == []

    def test_minus_strand_canonical(self):
        # - strand intron (s,e): donor = revcomp(genome[e-2:e]),
        # acceptor = revcomp(genome[s:s+2]). For GT-AG we need
        # genome[e-2:e]=revcomp("GT")="AC" and genome[s:s+2]=revcomp("AG")="CT".
        genome = "AAAAA" + "CT" + "CCCCCC" + "AC" + "AAAAA"
        cs = _make_set([_cand("m", [(5, 15)], "novel", strand="-")])
        _runner()._apply_canonical_gate(cs, genome)
        assert cs.candidate_ids() == ["m"]
