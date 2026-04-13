"""Pipeline orchestrator: interval -> candidates -> scoring -> EM -> quantification."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    aggregate_across_intervals,
    quantify_transcripts,
)
from fin.candidates.dataclasses import CandidateSet
from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import GenomicInterval, generate_isolated_intervals
from fin.pipeline.config import PipelineConfig
from fin.scoring.eventalign_parser import (
    ReadCandidateScore,
    build_distance_matrix,
    parse_eventalign_tsv,
)
from fin.scoring.external_tools import ExternalToolPaths, ExternalToolRunner
from fin.scoring.signal_dtw import compute_read_to_read_dtw, extract_signal_segments

logger = logging.getLogger(__name__)


class PipelineRunner:
    """Orchestrates the full pyfin pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._gtf_reader = None
        self._genome_fasta = None
        self._signal_reader = None
        self._tool_runner = None

    def setup(self):
        """Validate tools, open file handles, load references."""
        # Validate external tools
        tool_paths = ExternalToolPaths(
            minimap2=self.config.minimap2_path,
            f5c=self.config.f5c_path,
            samtools=self.config.samtools_path,
        )
        missing = tool_paths.validate()
        if missing:
            raise RuntimeError(f"Missing external tools: {', '.join(missing)}")
        self._tool_runner = ExternalToolRunner(tool_paths)

        # Load GTF
        if self.config.gtf_path:
            from fin.io.io_gtf import GTFReader

            self._gtf_reader = GTFReader(self.config.gtf_path)
            self._gtf_reader.open()
            self._gtf_reader.parse()

        # Load genome FASTA (full sequence per chrom)
        if self.config.genome_fasta_path:
            self._genome_fasta = self._load_genome_fasta(self.config.genome_fasta_path)

        # Open signal reader
        if self.config.signal_path:
            self._signal_reader = self._open_signal_reader()

        # Create work directory
        Path(self.config.work_dir).mkdir(parents=True, exist_ok=True)

        logger.info("Pipeline setup complete")

    def run(self) -> Dict[str, QuantResult]:
        """Run the full pipeline across all intervals.

        Returns:
            Aggregated quantification results.
        """
        # Generate intervals
        result = generate_isolated_intervals(
            self.config.bam_path,
            gtf_path=self.config.gtf_path,
            max_gap=self.config.max_gap,
            max_reads=self.config.max_reads,
        )
        intervals = result["intervals"]
        logger.info("Generated %d intervals", len(intervals))

        # Process each interval
        all_quant_results: List[List[QuantResult]] = []
        for i, interval in enumerate(intervals):
            logger.info("Processing interval %d/%d: %s", i + 1, len(intervals), interval.region_string)
            quant = self.process_interval(interval)
            if quant:
                all_quant_results.append(quant)

        # Aggregate across intervals
        aggregated = aggregate_across_intervals(all_quant_results)

        logger.info(
            "Pipeline complete: %d transcripts quantified", len(aggregated)
        )
        return aggregated

    def process_interval(
        self, interval: GenomicInterval
    ) -> Optional[List[QuantResult]]:
        """Process a single interval through all pipeline phases.

        Phase 1: Candidate discovery
        Phase 2: minimap2 + f5c scoring
        Phase 3: Signal DTW (read-to-read)
        Phase 4: EM assignment
        Phase 5: Quantification
        """
        work_dir = Path(self.config.work_dir) / interval.region_string.replace(":", "_").replace("-", "_")
        work_dir.mkdir(parents=True, exist_ok=True)

        # Get genome sequence for this chromosome
        chrom_seq = ""
        if self._genome_fasta and interval.chrom in self._genome_fasta:
            chrom_seq = self._genome_fasta[interval.chrom]

        # Phase 1: Candidate discovery
        candidate_set = discover_candidates(
            interval=interval,
            bam_path=self.config.bam_path,
            gtf_reader=self._gtf_reader,
            genome_fasta=chrom_seq,
            threshold=self.config.three_prime_threshold,
        )

        if candidate_set.num_candidates == 0:
            logger.info("No candidates for interval %s", interval.region_string)
            return None

        if not candidate_set.read_ids:
            logger.info("No reads for interval %s", interval.region_string)
            return None

        read_ids = sorted(candidate_set.read_ids)
        candidate_ids = candidate_set.candidate_ids()

        # Phase 2: External scoring (minimap2 + f5c)
        eventalign_tsv = self._tool_runner.run_scoring_pipeline(
            candidate_set=candidate_set,
            fastq_path=self.config.fastq_path,
            signal_path=self.config.signal_path,
            work_dir=str(work_dir),
        )

        # Parse eventalign output
        candidate_lengths = {
            c.candidate_id: len(c.sequence) for c in candidate_set.candidates
        }
        scores = parse_eventalign_tsv(str(eventalign_tsv), candidate_lengths)

        # Build distance matrix for EM
        dist_read_to_tx = build_distance_matrix(scores, read_ids, candidate_ids)

        # Phase 3: Signal DTW (read-to-read)
        segments = extract_signal_segments(
            scores, self._signal_reader, signal_format=self.config.signal_format
        )
        dist_read_to_read = compute_read_to_read_dtw(
            segments, read_ids, use_gpu=self.config.use_gpu
        )

        # Phase 4: EM assignment
        R, hard_assignments, _log_likelihoods = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=self.config.em_sigma,
            beta=self.config.em_beta,
            max_iter=self.config.em_max_iter,
            tol=self.config.em_tol,
            verbose=False,
        )

        # Phase 5: Quantification
        quant_results = quantify_transcripts(
            R, hard_assignments, candidate_set.candidates, read_ids
        )

        return quant_results

    def cleanup(self):
        """Close file handles."""
        if self._gtf_reader:
            self._gtf_reader.close()
        if self._signal_reader:
            self._signal_reader.close()

    def _load_genome_fasta(self, path: str) -> Dict[str, str]:
        """Load genome FASTA into a dict of chrom -> sequence."""
        from fin.io.io_fasta import FastaReader

        seqs = {}
        with FastaReader(path) as reader:
            for record in reader.iterate_records():
                seqs[record["header"].split()[0]] = record["sequence"]
        return seqs

    def _open_signal_reader(self):
        """Open the appropriate signal reader."""
        if self.config.signal_format == "pod5":
            from fin.io.io_pod5 import Pod5Reader

            reader = Pod5Reader(self.config.signal_path)
            reader.open()
            return reader
        else:
            from fin.io.io_slow5 import Slow5Reader

            reader = Slow5Reader(self.config.signal_path)
            reader.open()
            return reader
