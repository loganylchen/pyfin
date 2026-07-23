"""Wiring test for 5'-TSS short-isoform recovery inside ``_quant_cluster``.

The recovery CORE (peak-vs-ramp) is covered by ``test_isoform_recovery.py``.
This test pins the RUNNER wiring: a CandidateSet carrying a folded shadow
(``CandidateSet.shadows``) plus a BAM-span pileup that places a SHARP cluster of
the shadow reads' 5' ends at one position (real TSS) is re-activated as a NEW
novel QuantResult, and the parent's abundance drops by the recovered excess.

Signal-free (``signal_path=""`` -> no krill -> ``m2_resolve`` None): mono-exon
candidates (empty intron chain) so every read/member tie is deterministic and the
recovered short chain's spliced sequence is a plain genome slice.
"""
from __future__ import annotations

import random

import mappy
import pytest

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline import runner as runner_mod
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner


def _rnd(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


# A long random genome chromosome; the parent transcript is a mono-exon slice of
# it, and every read maps to the parent (they share its sequence). The shadow's
# short isoform is a plain sub-slice (empty intron chain -> genome slice).
_RNG = random.Random(2027)
CHROM_SEQ = _rnd(2000, _RNG)
PARENT_START = 300
PARENT_END = 1300
PARENT_SEQ = CHROM_SEQ[PARENT_START:PARENT_END]

TSS = 700  # sharp 5' peak position (genomic) where the short isoform starts


def _cand(cid, start, end):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=end,
        sequence=CHROM_SEQ[start:end],
        source="novel",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=start,
        end=end,
    )


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=2000, strand="+")


class _FakeRead:
    def __init__(self, name, start, end):
        self.query_name = name
        self.reference_start = start
        self.reference_end = end
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False


class _FakeBam:
    def __init__(self, reads):
        self._reads = reads

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def fetch(self, *a, **k):
        return iter(self._reads)


def _patch_bam(monkeypatch, reads):
    monkeypatch.setattr(
        runner_mod.pysam, "AlignmentFile", lambda *a, **k: _FakeBam(reads)
    )


def _runner():
    cfg = PipelineConfig(bam_path="x", quant_mode="cluster")
    r = PipelineRunner(cfg)
    # Provide the chromosome sequence so recovery can stitch the short isoform.
    r._genome_fasta = {"chr1": CHROM_SEQ}
    return r


@pytest.mark.skipif(mappy is None, reason="mappy required for cluster quant")
def test_recovery_activates_short_isoform(monkeypatch):
    parent_id = "c_parent"
    parent = _cand(parent_id, PARENT_START, PARENT_END)

    # Shadow reads: sharp 5'-end pileup exactly at TSS (real TSS peak).
    shadow_reads = [f"tss{i}" for i in range(14)]
    # Degradation reads: 5' ends scattered across the parent's 5' region (ramp),
    # placed in the recovery flanks (|d| in (window, 3*window] = (20, 60]).
    ramp_reads = [f"deg{i}" for i in range(6)]

    reads = []
    read_sequences = {}
    for rid in shadow_reads:
        # short isoform read: genomic [TSS, PARENT_END) sequence.
        reads.append(_FakeRead(rid, TSS, PARENT_END))
        read_sequences[rid] = CHROM_SEQ[TSS:PARENT_END]
    for k, rid in enumerate(ramp_reads):
        p5 = TSS + 40 + k  # in the flank band, not the peak window
        reads.append(_FakeRead(rid, p5, PARENT_END))
        read_sequences[rid] = CHROM_SEQ[p5:PARENT_END]

    all_ids = shadow_reads + ramp_reads
    parent.supporting_read_ids.update(all_ids)

    short_chain = ()  # mono-exon short isoform
    cs = CandidateSet(
        interval=_interval(),
        candidates=[parent],
        read_ids=set(all_ids),
        read_sequences=read_sequences,
        clusters=[[parent_id]],
        shadows={parent_id: [(short_chain, tuple(sorted(shadow_reads)))]},
    )
    _patch_bam(monkeypatch, reads)

    runner = _runner()
    results = runner._quant_cluster(cs, all_ids, _interval())
    out = {r.candidate_id: r for r in results}

    # A NEW novel candidate (not the parent) was emitted for the short isoform.
    novel_ids = [cid for cid in out if cid != parent_id]
    assert len(novel_ids) == 1, f"expected 1 recovered isoform, got {novel_ids}"
    recovered = out[novel_ids[0]]

    # The recovered isoform starts at the TSS and shares the parent's 3' end.
    assert recovered.source == "novel"
    assert recovered.start == TSS
    assert recovered.end == PARENT_END
    assert recovered.abundance > 8.0            # ~ the sharp-peak excess
    assert recovered.num_assigned_reads == len(shadow_reads)

    # Parent mass dropped by the recovered excess (mass conservation).
    total = len(all_ids)
    assert out[parent_id].abundance < total     # excess subtracted off the parent
    assert out[parent_id].abundance >= 0.0
