"""Pipeline configuration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class PipelineConfig:
    """Configuration for the full pyfin pipeline."""

    # Input files
    bam_path: str
    gtf_path: Optional[str] = None
    genome_fasta_path: str = ""
    fastq_path: str = ""
    signal_path: str = ""  # SLOW5/BLOW5 file
    signal_format: str = "slow5"  # "pod5" or "slow5"

    # Working directory
    work_dir: str = "./pyfin_work"

    # Candidate discovery
    three_prime_threshold: int = 24
    max_gap: int = 0

    # EM parameters
    em_sigma: float = 1.0
    em_beta: float = 0.5
    em_max_iter: int = 1000
    em_tol: float = 1e-4

    # DTW
    use_gpu: bool = True

    # External tools
    minimap2_path: str = "minimap2"
    f5c_path: str = "f5c"
    samtools_path: str = "samtools"

    # Parallelism (future)
    num_workers: int = 1

    # Limits
    max_reads: Optional[int] = None

    def validate(self):
        """Validate that required paths exist."""
        for attr in ("bam_path", "genome_fasta_path", "fastq_path", "signal_path"):
            val = getattr(self, attr)
            if val and not Path(val).exists():
                raise FileNotFoundError(f"{attr}: {val}")
