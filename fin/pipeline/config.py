"""Pipeline configuration."""

from __future__ import annotations

from dataclasses import dataclass
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
    min_novel_reads: int = 1          # A2: min supporting reads for a novel candidate (after collapsing)
    min_abundance: float = 0.0        # A4: drop quantified transcripts with abundance < this threshold
    min_max_r: float = 0.0            # Phase A T2: drop NOVEL transcripts whose max EM responsibility < this (Cohen's d=1.34)
    min_novel_combined_score: float = 0.0  # Step 3: drop NOVEL transcripts whose combined_score < this (F1-optimal: 0.288 with-GTF, 0.428 no-GTF; Cohen's d=0.70)

    # Canonical-motif alternative expansion: scan ±N bp around each read's
    # CIGAR-derived junction for GT-AG and emit additional chains.
    canonical_search_bp: int = 0           # 0 disables; 2 is the SIRV sweet spot
    max_chains_per_read: int = 16          # cap per-read alternative count (16 is SIRV sweet spot)

    # EM parameters
    em_sigma: float = 1.0
    em_beta: float = 0.5
    em_max_iter: int = 1000
    em_tol: float = 1e-4

    # EM prior / scoring
    score_alpha: float = 0.5          # weight for coherence vs discrimination in combined_score
    prior_weight_cap: float = 10.0    # max multiplicative boost from prior
    use_prior: bool = True            # apply combined_score-derived EM prior (set False for backward compat)

    # DTW
    use_gpu: bool = True
    max_reads_per_interval_for_dtw: int = 2000
    signal_normalize: bool = True   # per-read robust z-score before DTW

    # External tools
    f5c_path: str = "f5c"

    # Parallelism (future)
    num_workers: int = 1

    # Output
    output_gtf: Optional[str] = None
    output_tsv: Optional[str] = None        # scoring TSV output path
    output_bedpe: Optional[str] = None      # fusion BEDPE output path

    # Limits
    max_reads: Optional[int] = None

    # R-matrix persistence (T8: FP-by-EM data infrastructure)
    persist_R_matrix: bool = True  # write R.npy + R_meta.json per interval after EM

    # Fusion detection
    fusion_enabled: bool = False
    fusion_min_support: int = 2
    fusion_max_dist: int = 500     # bp window for breakpoint clustering
    fusion_flank_bp: int = 500     # bp to extract on each side of breakpoint for fusion sequence

    def validate(self):
        """Validate that required paths exist."""
        for attr in ("bam_path", "genome_fasta_path", "fastq_path", "signal_path"):
            val = getattr(self, attr)
            if val and not Path(val).exists():
                raise FileNotFoundError(f"{attr}: {val}")
