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
from fin.candidates.discovery import discover_candidates, merge_fusion_candidates
from fin.io.interval_manager import GenomicInterval, generate_isolated_intervals
from fin.pipeline.config import PipelineConfig
from fin.scoring.composite import (
    derive_prior_weights,
    populate_quant_scores,
    score_candidates_composite,
    subsample_reads_for_dtw,
)
from fin.scoring.eventalign_parser import (
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
        tool_paths = ExternalToolPaths(f5c=self.config.f5c_path)
        missing = tool_paths.validate()
        if missing:
            raise RuntimeError(f"Missing external tools: {', '.join(missing)}")

        # Create tool runner and build f5c index ONCE with all reads
        self._tool_runner = ExternalToolRunner(
            fastq_path=self.config.fastq_path,
            signal_path=self.config.signal_path,
            signal_format=self.config.signal_format,
            work_dir=self.config.work_dir,
            tools=tool_paths,
        )
        self._tool_runner.build_f5c_index()

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

        # Resolve gene_ids from GTF annotation
        for cid, qr in aggregated.items():
            if qr.source == "gtf" and self._gtf_reader:
                tx = self._gtf_reader.get_transcript(cid)
                if tx:
                    qr.gene_id = tx.gene_id
            if not qr.gene_id:
                qr.gene_id = qr.candidate_id

        # Write GTF output
        if self.config.output_gtf:
            from fin.io.io_gtf import write_gtf

            write_gtf(aggregated, self.config.output_gtf)
            logger.info("Wrote GTF output: %s", self.config.output_gtf)

        # US-013: Additional output writers
        if self.config.output_tsv:
            from fin.io.io_tsv import write_scoring_tsv

            transcript_lengths = {
                cid: sum(end - start for start, end in qr.exons) if qr.exons else 0
                for cid, qr in aggregated.items()
            }
            write_scoring_tsv(aggregated, transcript_lengths, self.config.output_tsv)
            logger.info("Wrote TSV output: %s", self.config.output_tsv)

        if self.config.fusion_enabled and self.config.output_bedpe:
            from fin.io.io_bedpe import write_fusion_bedpe

            write_fusion_bedpe(aggregated, self.config.output_bedpe)
            logger.info("Wrote BEDPE output: %s", self.config.output_bedpe)

        logger.info(
            "Pipeline complete: %d transcripts quantified", len(aggregated)
        )
        return aggregated

    def process_interval(
        self, interval: GenomicInterval
    ) -> Optional[List[QuantResult]]:
        """Process a single interval through all pipeline phases.

        Phase 1   : Candidate discovery (GTF + novel)
        Phase 1.5 : Fusion candidate augmentation (optional, fusion_enabled)
        Phase 2   : Per-candidate mappy + f5c eventalign
        Phase 3   : Signal DTW -> read-to-read + read-to-tx distance matrices
        Phase 4   : Composite scoring (coherence, discrimination, combined)
        Phase 5   : EM assignment with optional composite-derived prior
        Phase 6   : Probability-weighted quantification (scores populated)
        """
        work_dir = Path(self.config.work_dir) / interval.region_string.replace(":", "_").replace("-", "_")
        work_dir.mkdir(parents=True, exist_ok=True)

        # Get genome sequence for this chromosome
        chrom_seq = ""
        if self._genome_fasta and interval.chrom in self._genome_fasta:
            chrom_seq = self._genome_fasta[interval.chrom]

        # --- Phase 1: Candidate discovery ---
        candidate_set = discover_candidates(
            interval=interval,
            bam_path=self.config.bam_path,
            gtf_reader=self._gtf_reader,
            genome_fasta=chrom_seq,
            threshold=self.config.three_prime_threshold,
        )

        # --- Phase 1.5: Fusion candidate augmentation (optional) ---
        if self.config.fusion_enabled:
            candidate_set = self._augment_with_fusion_candidates(
                candidate_set, interval
            )

        if candidate_set.num_candidates == 0:
            logger.info("No candidates for interval %s", interval.region_string)
            return None

        if not candidate_set.read_ids:
            logger.info("No reads for interval %s", interval.region_string)
            return None

        # Sort once; canonical read axis for every distance matrix.
        read_ids = sorted(candidate_set.read_ids)
        candidate_ids = candidate_set.candidate_ids()

        # --- Phase 2: Per-candidate scoring (mappy + f5c eventalign) ---
        tsv_paths = self._tool_runner.score_candidates(candidate_set, work_dir)

        # Merge per-candidate eventalign results
        candidate_lengths = {
            c.candidate_id: len(c.sequence) for c in candidate_set.candidates
        }
        all_scores = []
        for tsv_path in tsv_paths:
            scores = parse_eventalign_tsv(str(tsv_path), candidate_lengths)
            all_scores.extend(scores)

        # DTW subsampling (m1): uniformly subsample reads when above cap.
        dtw_read_ids = subsample_reads_for_dtw(
            read_ids, self.config.max_reads_per_interval_for_dtw
        )
        if len(dtw_read_ids) != len(read_ids):
            logger.info(
                "Interval %s: subsampling DTW reads %d -> %d",
                interval.region_string,
                len(read_ids),
                len(dtw_read_ids),
            )

        # Read-to-tx distances (aligned with dtw_read_ids for m3 consistency).
        dist_read_to_tx = build_distance_matrix(
            all_scores, dtw_read_ids, candidate_ids
        )

        # --- Phase 3: Signal DTW (read-to-read) ---
        segments = extract_signal_segments(
            all_scores, self._signal_reader, signal_format=self.config.signal_format
        )
        dist_read_to_read = compute_read_to_read_dtw(
            segments, dtw_read_ids, use_gpu=self.config.use_gpu
        )

        # m3 alignment assertion: both distance matrices share the same read axis.
        assert (
            dist_read_to_tx.shape[0]
            == dist_read_to_read.shape[0]
            == len(dtw_read_ids)
        ), (
            f"Read axis mismatch: dist_read_to_tx rows={dist_read_to_tx.shape[0]}, "
            f"dist_read_to_read rows={dist_read_to_read.shape[0]}, "
            f"dtw_read_ids={len(dtw_read_ids)}"
        )
        assert (
            dist_read_to_read.shape[0] == dist_read_to_read.shape[1]
        ), "dist_read_to_read must be square"

        n_reads = len(dtw_read_ids)
        n_tx = len(candidate_ids)

        # --- Phase 4: Composite scoring (uniform R seed) ---
        R_uniform = np.full((n_reads, n_tx), 1.0 / max(n_tx, 1))
        composite_scores = score_candidates_composite(
            candidates=candidate_set.candidates,
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            R=R_uniform,
            alpha=self.config.score_alpha,
            use_gpu=self.config.use_gpu,
        )
        combined_scores_arr = np.array(
            [s.combined for s in composite_scores], dtype=float
        )

        # --- Phase 5: EM assignment with optional composite prior ---
        prior_weights: Optional[np.ndarray] = None
        if self.config.use_prior:
            prior_weights = derive_prior_weights(
                combined_scores_arr, n_tx, self.config.prior_weight_cap
            )

        R, hard_assignments, _log_likelihoods = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=self.config.em_sigma,
            beta=self.config.em_beta,
            max_iter=self.config.em_max_iter,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            prior_weights=prior_weights,
        )

        # --- Phase 6: Quantification + score field population ---
        quant_results = quantify_transcripts(
            R, hard_assignments, candidate_set.candidates, dtw_read_ids
        )
        populate_quant_scores(quant_results, composite_scores)
        return quant_results

    def _augment_with_fusion_candidates(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> CandidateSet:
        """Detect fusion breakpoints for the interval and merge into candidate_set."""
        from fin.fusion import (
            build_fusion_candidates,
            cluster_breakpoints,
            parse_sa_tags,
        )

        region = interval.region_string
        raw_bps = parse_sa_tags(
            self.config.bam_path, region=region, min_mapq=10
        )
        clusters = cluster_breakpoints(
            raw_bps,
            max_dist=self.config.fusion_max_dist,
            min_support=self.config.fusion_min_support,
        )
        if not clusters or build_fusion_candidates is None:
            return candidate_set

        fusion_cands = build_fusion_candidates(
            clusters,
            self._genome_fasta or {},
            flank_bp=self.config.fusion_flank_bp,
        )
        if not fusion_cands:
            return candidate_set

        return merge_fusion_candidates(candidate_set, fusion_cands)

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
