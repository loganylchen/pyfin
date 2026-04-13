"""Subprocess wrappers for minimap2, f5c, and samtools."""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from fin.candidates.dataclasses import CandidateSet

logger = logging.getLogger(__name__)


@dataclass
class ExternalToolPaths:
    """Paths to required external tools."""

    minimap2: str = "minimap2"
    f5c: str = "f5c"
    samtools: str = "samtools"

    def validate(self) -> List[str]:
        """Check that all tools are available on PATH.

        Returns:
            List of missing tool names (empty if all found).
        """
        missing = []
        for name in ("minimap2", "f5c", "samtools"):
            path = getattr(self, name)
            if shutil.which(path) is None:
                missing.append(name)
        return missing


def _run(cmd: List[str], description: str, **kwargs) -> subprocess.CompletedProcess:
    """Run a subprocess, raising on failure."""
    logger.info("Running %s: %s", description, " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        logger.error("%s failed (rc=%d): %s", description, result.returncode, result.stderr)
        raise RuntimeError(f"{description} failed: {result.stderr}")
    return result


class ExternalToolRunner:
    """Runs minimap2 + f5c scoring pipeline for a candidate set."""

    def __init__(self, tools: Optional[ExternalToolPaths] = None):
        self.tools = tools or ExternalToolPaths()

    def write_candidate_fasta(
        self, candidate_set: CandidateSet, work_dir: Path
    ) -> Path:
        """Write each candidate as a contig in a single FASTA file.

        Args:
            candidate_set: Candidates to write.
            work_dir: Working directory.

        Returns:
            Path to the candidates FASTA file.
        """
        fasta_path = work_dir / "candidates.fa"
        with open(fasta_path, "w") as f:
            for cand in candidate_set.candidates:
                if not cand.sequence:
                    logger.warning(
                        "Candidate %s has empty sequence, skipping", cand.candidate_id
                    )
                    continue
                f.write(f">{cand.candidate_id}\n")
                # Write sequence in 80-char lines
                seq = cand.sequence
                for i in range(0, len(seq), 80):
                    f.write(seq[i : i + 80] + "\n")
        return fasta_path

    def run_scoring_pipeline(
        self,
        candidate_set: CandidateSet,
        fastq_path: str,
        signal_path: str,
        work_dir: str,
    ) -> Path:
        """Run the full minimap2 + f5c eventalign pipeline.

        Args:
            candidate_set: Candidates (written as reference contigs).
            fastq_path: Path to reads FASTQ file.
            signal_path: Path to signal file (SLOW5/BLOW5).
            work_dir: Working directory for intermediate files.

        Returns:
            Path to the eventalign TSV output.
        """
        work = Path(work_dir)
        work.mkdir(parents=True, exist_ok=True)

        # 1. Write candidate FASTA
        candidates_fa = self.write_candidate_fasta(candidate_set, work)

        # 2. minimap2 align + samtools sort
        aligned_bam = work / "aligned.bam"
        minimap2_cmd = [
            self.tools.minimap2,
            "-ax",
            "map-ont",
            str(candidates_fa),
            str(fastq_path),
        ]
        sort_cmd = [
            self.tools.samtools,
            "sort",
            "-o",
            str(aligned_bam),
        ]
        logger.info("Running minimap2 | samtools sort")
        minimap2_proc = subprocess.Popen(
            minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        sort_proc = subprocess.Popen(
            sort_cmd,
            stdin=minimap2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        minimap2_proc.stdout.close()
        sort_stdout, sort_stderr = sort_proc.communicate()
        minimap2_proc.wait()

        if minimap2_proc.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {minimap2_proc.stderr.read()}")
        if sort_proc.returncode != 0:
            raise RuntimeError(f"samtools sort failed: {sort_stderr.decode()}")

        # 3. samtools index
        _run(
            [self.tools.samtools, "index", str(aligned_bam)],
            "samtools index",
        )

        # 4. f5c index
        _run(
            [
                self.tools.f5c,
                "index",
                "--slow5",
                str(signal_path),
                str(fastq_path),
            ],
            "f5c index",
        )

        # 5. f5c eventalign
        eventalign_tsv = work / "eventalign.tsv"
        eventalign_cmd = [
            self.tools.f5c,
            "eventalign",
            "--rna",
            "--reads",
            str(fastq_path),
            "--bam",
            str(aligned_bam),
            "--genome",
            str(candidates_fa),
            "--slow5",
            str(signal_path),
            "--signal-index",
            "--scale-events",
            "--print-read-names",
            "--samples",
        ]
        result = _run(eventalign_cmd, "f5c eventalign")

        with open(eventalign_tsv, "w") as f:
            f.write(result.stdout)

        logger.info("Eventalign output: %s", eventalign_tsv)
        return eventalign_tsv
