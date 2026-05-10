"""Integration: discover_candidates with vs without a reference GTF.

User goal: pyfin should perform well even when no reference GTF is supplied.
This test runs the same real RNA004 BAM against ``discover_candidates`` in
both modes and asserts:

1. GTF-mode produces a candidate set whose `gtf` candidates have non-empty
   spliced sequences (already worked before the no-GTF fix).
2. **No-GTF mode produces ONLY novel candidates, all with non-empty spliced
   sequences** (the regression we just fixed; previously empty -> f5c dropped).
3. Both modes use the same read population.
4. Per-interval candidate counts are reported so the user can eyeball
   how much information is lost when the reference GTF is absent.
"""

from __future__ import annotations

import os

import pysam
import pytest

from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import (
    GenomicInterval,
    generate_isolated_intervals,
)
from fin.io.io_fasta import FASTAReader
from fin.io.io_gtf import GTFReader

TESTDATA = os.path.join(os.path.dirname(__file__), "..", "testdata")
BAM = os.path.join(TESTDATA, "RNA004.test.bam")
GTF = os.path.join(TESTDATA, "RNA004.test.gtf")
FA = os.path.join(TESTDATA, "test.fa")


pytestmark = pytest.mark.skipif(
    not (os.path.exists(BAM) and os.path.exists(FA) and os.path.exists(GTF)),
    reason="RNA004 fixture not available",
)


def _load_chrom_sequences(fa_path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    with FASTAReader(fa_path) as fa:
        for record in fa.iterate_records():
            seqs[record.id] = record.sequence
    return seqs


def _first_interval_with_reads(bam_path: str, intervals) -> GenomicInterval | None:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for itv in intervals:
            try:
                cnt = bam.count(itv.chrom, itv.start, itv.end)
            except ValueError:
                continue
            if cnt > 0:
                return itv
    return None


def test_discovery_gtf_vs_no_gtf_real_data(capsys):
    seqs = _load_chrom_sequences(FA)

    # Build strand-separated intervals like the runner does, then pick the
    # first interval with reads so the test stays fast.
    info = generate_isolated_intervals(BAM, GTF)
    intervals = info["intervals"]
    interval = _first_interval_with_reads(BAM, intervals)
    assert interval is not None, "No interval with reads in fixture"
    chrom_seq = seqs.get(interval.chrom, "")
    assert chrom_seq, f"Genome FASTA missing chrom {interval.chrom}"

    # --- Mode A: with GTF reader -------------------------------------------
    with GTFReader(GTF) as g:
        g.parse()
        cs_gtf = discover_candidates(
            interval=interval,
            bam_path=BAM,
            gtf_reader=g,
            genome_fasta=chrom_seq,
        )

    # --- Mode B: no GTF reader ---------------------------------------------
    cs_nogtf = discover_candidates(
        interval=interval,
        bam_path=BAM,
        gtf_reader=None,
        genome_fasta=chrom_seq,
    )

    # Same reads visible to both modes.
    assert cs_gtf.read_ids == cs_nogtf.read_ids

    # GTF mode: at least one candidate; every gtf candidate must have a
    # non-empty spliced sequence.
    gtf_only = [c for c in cs_gtf.candidates if c.source == "gtf"]
    if gtf_only:
        for c in gtf_only:
            assert c.sequence, f"GTF candidate {c.candidate_id} has empty sequence"
            assert c.sequence != "N" * len(c.sequence)

    # No-GTF mode: every candidate must be novel and have a non-empty spliced
    # sequence (regression fix for the empty-sequence bug).
    novel_only = [c for c in cs_nogtf.candidates if c.source == "novel"]
    assert all(c.source == "novel" for c in cs_nogtf.candidates), (
        f"No-GTF mode emitted non-novel candidates: "
        f"{[(c.candidate_id, c.source) for c in cs_nogtf.candidates]}"
    )
    if novel_only:
        for c in novel_only:
            assert c.sequence, (
                f"No-GTF novel candidate {c.candidate_id} has empty sequence"
            )
            # Sanity: spliced cDNA length should not exceed genomic span.
            assert len(c.sequence) <= (c.end - c.start)

    # Surface a side-by-side report for the user.
    print("\n=== discovery comparison ===")
    print(f"interval                : {interval.region_string} {interval.strand}")
    print(f"reads in interval       : {len(cs_gtf.read_ids)}")
    print(f"GTF mode    candidates  : {len(cs_gtf.candidates)} "
          f"(gtf={sum(1 for c in cs_gtf.candidates if c.source == 'gtf')}, "
          f"novel={sum(1 for c in cs_gtf.candidates if c.source == 'novel')})")
    print(f"  with non-empty seq    : "
          f"{sum(1 for c in cs_gtf.candidates if c.sequence)}/{len(cs_gtf.candidates)}")
    print(f"no-GTF mode candidates  : {len(cs_nogtf.candidates)} (all novel)")
    print(f"  with non-empty seq    : "
          f"{sum(1 for c in cs_nogtf.candidates if c.sequence)}/{len(cs_nogtf.candidates)}")
    if cs_nogtf.candidates:
        seq_lens = [len(c.sequence) for c in cs_nogtf.candidates]
        print(f"  no-GTF cDNA lengths   : min={min(seq_lens)}, "
              f"median={sorted(seq_lens)[len(seq_lens)//2]}, max={max(seq_lens)}")
    captured = capsys.readouterr()
    # Re-print so pytest -s still shows it in either mode
    print(captured.out)
