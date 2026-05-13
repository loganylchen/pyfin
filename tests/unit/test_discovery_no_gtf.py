"""Tests for novel candidate discovery without a reference GTF (de novo mode).

When the user runs pyfin without --gtf, discovery must:
1. Produce novel candidates with **non-empty spliced sequences** stitched from
   the genome FASTA + CIGAR-derived intron chain.
2. Apply reverse-complement on the '-' strand.
3. Skip groups where genome bounds are unusable instead of emitting empty
   candidates that f5c will silently drop.
"""

from __future__ import annotations

import os

import pysam

from fin.candidates.dataclasses import IntronChain
from fin.candidates.discovery import (
    _build_spliced_sequence,
    _exons_from_chain,
    _reverse_complement,
    discover_candidates,
)
from fin.io.interval_manager import GenomicInterval


# ---------------------------------------------------------------------------
# Helpers for the BAM fixture
# ---------------------------------------------------------------------------


def _write_bam_with_spliced_reads(tmp_path: str, chrom: str, chrom_len: int,
                                   reads: list[dict]) -> str:
    """Write a tiny BAM with the given reads (cigar, pos, strand)."""
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": chrom, "LN": chrom_len}],
        }
    )
    bam_path = os.path.join(tmp_path, "spliced.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for r in reads:
            seg = pysam.AlignedSegment(header)
            seg.query_name = r["name"]
            # Build a query sequence whose length matches the sum of M ops.
            mlen = sum(length for op, length in r["cigar"] if op in (0, 1, 4, 7, 8))
            seg.query_sequence = "A" * mlen
            seg.flag = 0x10 if r.get("reverse", False) else 0
            seg.reference_id = 0
            seg.reference_start = r["pos"]
            seg.mapping_quality = 60
            seg.cigar = r["cigar"]
            seg.query_qualities = pysam.qualitystring_to_array("I" * mlen)
            bam.write(seg)
    sorted_path = bam_path + ".sorted.bam"
    pysam.sort("-o", sorted_path, bam_path)
    pysam.index(sorted_path)
    return sorted_path


# ---------------------------------------------------------------------------
# Pure-function tests
# ---------------------------------------------------------------------------


class TestExonsFromChain:
    def test_single_exon(self):
        ex = _exons_from_chain(100, 200, IntronChain(introns=()))
        assert ex == [(100, 200)]

    def test_two_exons_one_intron(self):
        chain = IntronChain(introns=((150, 170),))
        assert _exons_from_chain(100, 200, chain) == [(100, 150), (170, 200)]

    def test_three_exons_two_introns(self):
        chain = IntronChain(introns=((120, 140), (160, 180)))
        assert _exons_from_chain(100, 200, chain) == [
            (100, 120),
            (140, 160),
            (180, 200),
        ]


class TestReverseComplement:
    def test_basic(self):
        assert _reverse_complement("ACGT") == "ACGT"
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("GATTACA") == "TGTAATC"

    def test_preserves_n(self):
        assert _reverse_complement("ANNT") == "ANNT"


class TestBuildSplicedSequence:
    def test_single_exon_plus_strand(self):
        # 30-base genome, no introns -> exact slice
        genome = "ACGT" * 8  # 32 bases
        seq = _build_spliced_sequence(genome, 4, 12, IntronChain(introns=()), "+")
        assert seq == genome[4:12]

    def test_spliced_two_exon_plus_strand(self):
        # Build a known-content genome where the intron has different bases
        genome = "A" * 10 + "C" * 10 + "G" * 10 + "T" * 10  # 40 bases
        # Exon1 [5,10) = "AAAAA"; intron [10,20) = "CCCCCCCCCC";
        # Exon2 [20,25) = "GGGGG". Spliced should be "AAAAAGGGGG".
        chain = IntronChain(introns=((10, 20),))
        seq = _build_spliced_sequence(genome, 5, 25, chain, "+")
        assert seq == "AAAAAGGGGG"
        assert "C" not in seq  # intron stripped

    def test_negative_strand_reverse_complements(self):
        genome = "A" * 10 + "C" * 10 + "G" * 10 + "T" * 10
        chain = IntronChain(introns=((10, 20),))
        # Plus strand spliced = "AAAAAGGGGG" -> RC = "CCCCCTTTTT"
        seq = _build_spliced_sequence(genome, 5, 25, chain, "-")
        assert seq == "CCCCCTTTTT"

    def test_empty_genome_returns_empty(self):
        assert _build_spliced_sequence("", 0, 100, IntronChain(introns=()), "+") == ""

    def test_clips_to_chromosome_bounds(self):
        # End beyond chromosome should not crash; just clip.
        genome = "ACGT" * 5  # 20 bases
        seq = _build_spliced_sequence(genome, 10, 100, IntronChain(introns=()), "+")
        assert seq == genome[10:]  # clipped at chrom end


# ---------------------------------------------------------------------------
# End-to-end no-GTF discovery test
# ---------------------------------------------------------------------------


class TestDiscoverCandidatesNoGTF:
    def test_novel_candidates_have_nonempty_sequences(self, tmp_path):
        """gtf_reader=None: every novel candidate must carry a non-empty,
        correctly spliced sequence assembled from the genome FASTA."""
        # Genome with 4 distinct base regions for easy verification:
        # [0,100)=A, [100,200)=C, [200,300)=G, [300,400)=T
        genome = "A" * 100 + "C" * 100 + "G" * 100 + "T" * 100

        # Three reads with the SAME spliced structure: exon[20,100) +
        # intron[100,200) + exon[200,280) on chr1, forward strand.
        # CIGAR per read: 80M 100N 80M
        reads = [
            {"name": f"r{i}", "pos": 20, "reverse": False,
             "cigar": [(0, 80), (3, 100), (0, 80)]}
            for i in range(3)
        ]
        bam_path = _write_bam_with_spliced_reads(
            str(tmp_path), "chr1", chrom_len=400, reads=reads
        )

        interval = GenomicInterval(chrom="chr1", start=0, end=400, strand="+")
        result = discover_candidates(
            interval=interval,
            bam_path=bam_path,
            gtf_reader=None,         # explicitly no GTF
            genome_fasta=genome,
        )

        # All candidates must be novel and have non-empty sequences.
        assert len(result.candidates) >= 1
        for c in result.candidates:
            assert c.source == "novel"
            assert c.sequence, f"Empty sequence for novel candidate {c.candidate_id}"
            # Spliced length = sum of exon lengths, never including intron bases.
            assert "C" not in c.sequence, (
                "Intron bases (C) leaked into spliced sequence: "
                f"{c.sequence[:50]}..."
            )
            # Should be exactly 80 A's + 80 G's = 160 bases.
            assert len(c.sequence) == 160
            assert c.sequence == "A" * 80 + "G" * 80

    def test_junction_snap_consensus(self, tmp_path):
        """With junction_snap_bp=6, noisy donor/acceptor positions across reads
        should collapse to a single consensus chain instead of splitting into
        multiple near-identical novel candidates.
        """
        # Genome: 400 bases of N so no canonical motifs interfere; snap
        # falls back to mode/median consensus.
        genome = "N" * 400

        # Three reads with slightly different intron coordinates (within 6 bp):
        # r1: CIGAR 80M 100N 80M -> intron [100, 200)
        # r2: shifted by +3      -> intron [103, 200)  (5' shift of 3)
        # r3: shifted by +6      -> intron [106, 200)
        reads = [
            {"name": "r1", "pos": 20, "reverse": False,
             "cigar": [(0, 80), (3, 100), (0, 80)]},
            {"name": "r2", "pos": 20, "reverse": False,
             "cigar": [(0, 83), (3, 97), (0, 80)]},
            {"name": "r3", "pos": 20, "reverse": False,
             "cigar": [(0, 86), (3, 94), (0, 80)]},
        ]
        bam_path = _write_bam_with_spliced_reads(
            str(tmp_path), "chr1", chrom_len=400, reads=reads
        )

        interval = GenomicInterval(chrom="chr1", start=0, end=400, strand="+")

        # Without snap, three reads with three different chains -> 3 candidates.
        no_snap = discover_candidates(
            interval=interval,
            bam_path=bam_path,
            gtf_reader=None,
            genome_fasta=genome,
            junction_snap_bp=0,
        )
        novel_no_snap = [c for c in no_snap.candidates if c.source == "novel"]
        assert len(novel_no_snap) == 3

        # With snap_bp=6, all three should collapse to one consensus chain.
        snapped = discover_candidates(
            interval=interval,
            bam_path=bam_path,
            gtf_reader=None,
            genome_fasta=genome,
            junction_snap_bp=6,
            min_junction_support=2,
        )
        novel_snapped = [c for c in snapped.candidates if c.source == "novel"]
        assert len(novel_snapped) == 1
        cand = novel_snapped[0]
        # Donor consensus: cluster {100, 103, 106}; mode tie -> n//2 = index 1
        # of sorted unique [100, 103, 106] = 103.
        # Acceptor cluster {200, 200, 200} -> 200.
        assert cand.intron_chain.introns == ((103, 200),)
        # All three reads should support the snapped candidate.
        assert len(cand.supporting_read_ids) == 3

    def test_novel_candidate_negative_strand_rc(self, tmp_path):
        """Reverse-strand reads must yield reverse-complemented spliced cDNA."""
        genome = "A" * 100 + "C" * 100 + "G" * 100 + "T" * 100
        reads = [
            {"name": f"r{i}", "pos": 20, "reverse": True,
             "cigar": [(0, 80), (3, 100), (0, 80)]}
            for i in range(2)
        ]
        bam_path = _write_bam_with_spliced_reads(
            str(tmp_path), "chr1", chrom_len=400, reads=reads
        )

        interval = GenomicInterval(chrom="chr1", start=0, end=400, strand="-")
        result = discover_candidates(
            interval=interval,
            bam_path=bam_path,
            gtf_reader=None,
            genome_fasta=genome,
        )

        assert len(result.candidates) >= 1
        for c in result.candidates:
            assert c.source == "novel"
            assert c.sequence
            # Plus-strand spliced = AAAAA...GGGGG (80 A + 80 G).
            # Reverse complement = CCCCC...TTTTT (80 C + 80 T).
            assert c.sequence == "C" * 80 + "T" * 80
