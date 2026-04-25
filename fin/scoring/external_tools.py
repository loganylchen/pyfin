"""Per-candidate scoring pipeline using mappy + f5c eventalign."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from fin.candidates.dataclasses import CandidateSet, TranscriptCandidate

logger = logging.getLogger(__name__)


@dataclass
class ExternalToolPaths:
    """Paths to required external tools."""

    f5c: str = "f5c"
    samtools: str = "samtools"

    def validate(self) -> List[str]:
        """Check that all tools are available on PATH.

        Returns:
            List of missing tool names (empty if all found).
        """
        missing = []
        for name in ("f5c", "samtools"):
            path = getattr(self, name)
            if shutil.which(path) is None:
                missing.append(name)
        return missing


def _signal_format_args(signal_format: str, signal_path: Path) -> List[str]:
    """Return f5c CLI args for the given signal format.

    f5c accepts ``--slow5`` for SLOW5/BLOW5 inputs and ``--pod5`` for POD5
    inputs. Both flags require the path to the signal file.
    """
    fmt = (signal_format or "").lower()
    if fmt in ("slow5", "blow5"):
        return ["--slow5", str(signal_path)]
    if fmt == "pod5":
        return ["--pod5", str(signal_path)]
    raise ValueError(
        f"Unsupported signal_format '{fmt}' (expected 'slow5', 'blow5', or 'pod5')"
    )


def _run(cmd: List[str], description: str, **kwargs) -> subprocess.CompletedProcess:
    """Run a subprocess, raising on failure."""
    logger.info("Running %s: %s", description, " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        logger.error("%s failed (rc=%d): %s", description, result.returncode, result.stderr)
        raise RuntimeError(f"{description} failed: {result.stderr}")
    return result


class ExternalToolRunner:
    """Runs per-candidate mappy alignment + f5c eventalign scoring."""

    def __init__(
        self,
        fastq_path: str,
        signal_path: str,
        signal_format: str,
        work_dir: str,
        tools: Optional[ExternalToolPaths] = None,
    ):
        self.fastq_path = Path(fastq_path)
        self.signal_path = Path(signal_path)
        self.signal_format = signal_format
        self.work_dir = Path(work_dir)
        self.tools = tools or ExternalToolPaths()
        self._f5c_indexed = False

    def build_f5c_index(self):
        """Build f5c index ONCE with all reads. Called at pipeline setup."""
        cmd = [self.tools.f5c, "index", str(self.fastq_path)]
        cmd.extend(_signal_format_args(self.signal_format, self.signal_path))
        _run(cmd, "f5c index")
        self._f5c_indexed = True
        logger.info("f5c index built for %s", self.fastq_path)

    def score_candidates(
        self, candidate_set: CandidateSet, interval_work_dir: Path
    ) -> List[Path]:
        """Score all reads against each candidate separately.

        For each candidate:
        1. Write single-contig FASTA
        2. mappy alignment -> write BAM via pysam
        3. f5c eventalign -> TSV

        Returns list of eventalign TSV paths (one per candidate).
        """
        assert self._f5c_indexed, "Call build_f5c_index() first"

        interval_work_dir.mkdir(parents=True, exist_ok=True)
        tsv_paths = []

        for candidate in candidate_set.candidates:
            if not candidate.sequence:
                logger.warning(
                    "Candidate %s has empty sequence, skipping", candidate.candidate_id
                )
                continue

            cand_dir = interval_work_dir / candidate.candidate_id
            cand_dir.mkdir(parents=True, exist_ok=True)

            # 1. Write single-candidate FASTA
            fasta_path = cand_dir / "candidate.fa"
            write_single_candidate_fasta(candidate, fasta_path)

            # 2. Align ALL reads to this candidate with mappy
            bam_path = cand_dir / "aligned.bam"
            align_reads_with_mappy(candidate, self.fastq_path, bam_path)

            # 3. f5c eventalign
            tsv_path = cand_dir / "eventalign.tsv"
            run_f5c_eventalign(
                reads_fq=self.fastq_path,
                bam=bam_path,
                genome_fa=fasta_path,
                signal=self.signal_path,
                signal_format=self.signal_format,
                output=tsv_path,
                f5c_path=self.tools.f5c,
            )
            tsv_paths.append(tsv_path)

        return tsv_paths


def write_single_candidate_fasta(
    candidate: TranscriptCandidate, fasta_path: Path
) -> None:
    """Write a single candidate as a FASTA file (one contig)."""
    with open(fasta_path, "w") as f:
        seq = candidate.sequence
        f.write(f">{candidate.candidate_id}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def align_reads_with_mappy(
    candidate: TranscriptCandidate, fastq_path: Path, bam_path: Path
) -> None:
    """Align all reads to a single candidate using mappy, write sorted+indexed BAM."""
    import mappy
    import pysam

    # Build mappy index for this candidate's sequence
    aligner = mappy.Aligner(seq=candidate.sequence, preset="map-ont")
    if not aligner:
        raise RuntimeError(
            f"Failed to build mappy index for {candidate.candidate_id}"
        )

    # Parse FASTQ, align each read, collect alignments
    alignments = []
    for name, seq, qual in mappy.fastx_read(str(fastq_path)):
        for hit in aligner.map(seq):
            if hit.is_primary:
                alignments.append((name, seq, qual, hit))

    logger.info(
        "mappy aligned %d reads to candidate %s",
        len(alignments),
        candidate.candidate_id,
    )

    # Write unsorted BAM first, then sort
    unsorted_bam = str(bam_path) + ".unsorted.bam"
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [
                {"SN": candidate.candidate_id, "LN": len(candidate.sequence)}
            ],
        }
    )

    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as outf:
        for name, seq, qual, hit in alignments:
            a = pysam.AlignedSegment(header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = 0 if hit.strand == 1 else 16
            a.reference_id = 0  # single contig
            a.reference_start = hit.r_st
            a.mapping_quality = hit.mapq
            # mappy returns cigar as (length, op); pysam wants (op, length)
            a.cigar = [(op, length) for length, op in hit.cigar]
            if qual:
                a.query_qualities = pysam.qualitystring_to_array(qual)
            outf.write(a)

    # Sort and index
    pysam.sort("-o", str(bam_path), unsorted_bam)
    os.remove(unsorted_bam)
    pysam.index(str(bam_path))


def run_f5c_eventalign(
    reads_fq: Path,
    bam: Path,
    genome_fa: Path,
    signal: Path,
    signal_format: str,
    output: Path,
    f5c_path: str = "f5c",
) -> None:
    """Run f5c eventalign for a single candidate."""
    cmd = [
        f5c_path,
        "eventalign",
        "--rna",
        "--reads",
        str(reads_fq),
        "--bam",
        str(bam),
        "--genome",
        str(genome_fa),
        "--signal-index",
        "--scale-events",
        "--print-read-names",
        "--samples",
    ]
    cmd.extend(_signal_format_args(signal_format, signal))

    result = _run(cmd, f"f5c eventalign ({genome_fa.parent.name})")

    with open(output, "w") as f:
        f.write(result.stdout)

    logger.info("Eventalign output: %s", output)
