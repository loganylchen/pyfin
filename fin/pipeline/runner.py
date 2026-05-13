"""Pipeline orchestrator: interval -> candidates -> scoring -> EM -> quantification."""

from __future__ import annotations

import json
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
from fin.scoring.diff_region_dtw import compute_diff_region_m4
from fin.scoring.eventalign_parser import (
    build_distance_matrix,
    parse_eventalign_tsv,
)
from fin.scoring.external_tools import ExternalToolPaths, ExternalToolRunner
from fin.scoring.signal_dtw import compute_read_to_read_dtw, extract_signal_segments

logger = logging.getLogger(__name__)


def _project_responsibilities_full(
    R_sub: np.ndarray,
    sub_read_ids: List[str],
    full_read_ids: List[str],
    dist_read_to_tx_full: np.ndarray,
    sigma: float,
    prior_weights: Optional[np.ndarray],
) -> tuple:
    """Project EM responsibilities from subsampled to full read set.

    Subsampled reads keep their EM responsibilities (which include coherence).
    Non-subsampled reads get a coherence-free softmax over the full d_tx using
    the same sigma + (optional) prior. This restores correct quantification
    when DTW subsampling is in effect.
    """
    n_full = len(full_read_ids)
    n_tx = R_sub.shape[1]
    sub_index = {rid: i for i, rid in enumerate(sub_read_ids)}

    R_full = np.zeros((n_full, n_tx), dtype=R_sub.dtype)

    # Coherence-free softmax for non-subsampled reads (numerically stable).
    d_tx = dist_read_to_tx_full
    d_min = d_tx.min(axis=1, keepdims=True)
    R_softmax = np.exp(-(d_tx - d_min) / max(sigma, 1e-6))
    if prior_weights is not None:
        R_softmax = R_softmax * np.asarray(prior_weights).reshape(1, -1)
    row_sums = np.maximum(R_softmax.sum(axis=1, keepdims=True), 1e-10)
    R_softmax = R_softmax / row_sums

    for i, rid in enumerate(full_read_ids):
        si = sub_index.get(rid)
        if si is not None:
            R_full[i, :] = R_sub[si, :]
        else:
            R_full[i, :] = R_softmax[i, :]

    hard_full = np.argmax(R_full, axis=1)
    return R_full, hard_full


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

        # R5 ablation: enable_score_filter=False disables all post-EM filters so
        # the ablation row sees unfiltered EM output (AC7).
        _score_filter_on = getattr(self.config, "enable_score_filter", True)

        # A4: post-EM abundance filter. GTF-sourced and fusion transcripts are
        # exempt so users can still see zero-abundance annotated entries
        # (useful for debugging coverage gaps) and fusion calls are not
        # silently dropped by NOVEL-targeted filters. Only NOVEL candidates
        # are dropped here.
        if _score_filter_on and self.config.min_abundance > 0.0:
            before = len(aggregated)
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if qr.source in ("gtf", "fusion")
                or qr.abundance >= self.config.min_abundance
            }
            dropped = before - len(aggregated)
            if dropped:
                logger.info(
                    "Dropped %d novel transcripts with abundance < %.3f",
                    dropped,
                    self.config.min_abundance,
                )

        # Phase A Tier-2: max_R filter (FP-by-EM). EM responsibility is a
        # strong FP discriminator (Cohen's d=1.34 vs combined_score's 0.70).
        # GTF-sourced and fusion transcripts are exempt to preserve annotation
        # visibility and avoid silently dropping fusion calls.
        if _score_filter_on and self.config.min_max_r > 0.0:
            before = len(aggregated)
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if qr.source in ("gtf", "fusion")
                or qr.max_R >= self.config.min_max_r
            }
            dropped = before - len(aggregated)
            if dropped:
                logger.info(
                    "Dropped %d novel transcripts with max_R < %.3f",
                    dropped,
                    self.config.min_max_r,
                )

        # Step 3: novel combined_score filter (FP-by-scoring). combined_score is
        # the geometric mean of coherence and discrimination (Cohen's d=0.70).
        # The F1-optimal threshold from profile_fp.md is 0.288 for with-GTF and
        # 0.428 for no-GTF runs. GTF-sourced and fusion transcripts are exempt.
        if _score_filter_on and self.config.min_novel_combined_score > 0.0:
            before = len(aggregated)
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if qr.source in ("gtf", "fusion")
                or qr.combined_score >= self.config.min_novel_combined_score
            }
            dropped = before - len(aggregated)
            if dropped:
                logger.info(
                    "Dropped %d novel transcripts with combined_score < %.3f",
                    dropped,
                    self.config.min_novel_combined_score,
                )

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
            min_novel_reads=self.config.min_novel_reads,
            canonical_search_bp=self.config.canonical_search_bp,
            max_chains_per_read=self.config.max_chains_per_read,
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

        # R1 ablation: enable_signal=False → skip signal/EM, use mappy argmax.
        if not getattr(self.config, "enable_signal", True):
            return self._process_interval_mappy_argmax(
                candidate_set, read_ids, interval
            )

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
        # Also build the full-read d_tx so quantification can include
        # non-subsampled reads (P0-1: prevent abundance underestimation when
        # max_reads_per_interval_for_dtw caps the EM input).
        dist_read_to_tx_full = build_distance_matrix(
            all_scores, read_ids, candidate_ids
        )

        # --- Phase 3: Signal DTW (read-to-read) ---
        # m4_source controls which read-to-read matrix is used (AC1).
        m4_src = getattr(self.config, "m4_source", "whole_read")
        if m4_src == "none":
            # AC5: zero matrix (no coherence); not None so EM shape checks pass.
            dist_read_to_read = np.zeros(
                (len(dtw_read_ids), len(dtw_read_ids)), dtype=np.float32
            )
        elif m4_src == "diff_region":
            # Collect per-event records needed by diff-region DTW (AC7-pre).
            scores_with_events = []
            cand_lengths_map = {
                c.candidate_id: len(c.sequence) for c in candidate_set.candidates
            }
            for tsv_path in tsv_paths:
                ev_scores = parse_eventalign_tsv(
                    str(tsv_path), cand_lengths_map, collect_events=True
                )
                scores_with_events.extend(ev_scores)
            scores_by_pair = {
                (s.read_name, s.candidate_id): s for s in scores_with_events
            }
            from fin.ablation.runner import _nan_to_row_mean

            m4_raw = compute_diff_region_m4(
                read_ids=dtw_read_ids,
                candidates=candidate_set.candidates,
                scores_by_pair=scores_by_pair,
                signal_reader=self._signal_reader,
                interval_start=interval.start,
                interval_end=interval.end,
                signal_format=self.config.signal_format,
                use_gpu=self.config.use_gpu,
                normalize=self.config.signal_normalize,
            )
            # AC8: sanitize NaN rows before EM.
            dist_read_to_read = _nan_to_row_mean(m4_raw)
        else:
            # Default: whole-read signal DTW.
            segments = extract_signal_segments(
                all_scores,
                self._signal_reader,
                signal_format=self.config.signal_format,
                normalize=self.config.signal_normalize,
            )
            dist_read_to_read = compute_read_to_read_dtw(
                segments, dtw_read_ids, use_gpu=self.config.use_gpu
            )

        # m3 alignment check: both distance matrices share the same read axis.
        if not (
            dist_read_to_tx.shape[0]
            == dist_read_to_read.shape[0]
            == len(dtw_read_ids)
        ):
            raise ValueError(
                f"Read axis mismatch: dist_read_to_tx rows={dist_read_to_tx.shape[0]}, "
                f"dist_read_to_read rows={dist_read_to_read.shape[0]}, "
                f"dtw_read_ids={len(dtw_read_ids)}"
            )
        if dist_read_to_read.shape[0] != dist_read_to_read.shape[1]:
            raise ValueError("dist_read_to_read must be square")

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

        # Data-adaptive sigma (P0-α): with per-event normalized d_tx the
        # absolute scale is dataset-dependent (~0.5-5 nats). A fixed sigma=1.0
        # often collapses R to one-hot, defeating EM. Use the median per-read
        # range as a robust scale; clip into [min_sigma, em_sigma_max].
        if n_tx >= 2 and n_reads > 0:
            d_max = dist_read_to_tx.max(axis=1)
            d_min = dist_read_to_tx.min(axis=1)
            adaptive = float(np.median(d_max - d_min))
            sigma_use = float(np.clip(
                adaptive if adaptive > 0 else self.config.em_sigma,
                getattr(self.config, "em_sigma_min", 0.05),
                getattr(self.config, "em_sigma_max", 50.0),
            ))
        else:
            sigma_use = self.config.em_sigma

        # R2 ablation: em_max_iter_override=1 forces single-step EM.
        em_max_iter_use = (
            self.config.em_max_iter_override
            if getattr(self.config, "em_max_iter_override", None) is not None
            else self.config.em_max_iter
        )

        R, hard_assignments, _log_likelihoods = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=sigma_use,
            beta=self.config.em_beta,
            max_iter=em_max_iter_use,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            prior_weights=prior_weights,
        )

        # --- Phase 6: Quantification + score field population ---
        # P0-1: when DTW subsampled, project EM responsibilities back to the
        # full read set. Subsampled reads keep their EM responsibilities;
        # non-subsampled reads get a coherence-free softmax over d_tx_full
        # using the same sigma. This avoids dropping ~maxN reads from the
        # quantification, which previously caused systematic abundance
        # under-estimation when an interval had >max_reads_per_interval_for_dtw.
        if len(dtw_read_ids) == len(read_ids):
            R_quant = R
            hard_quant = hard_assignments
            quant_read_ids = dtw_read_ids
        else:
            R_quant, hard_quant = _project_responsibilities_full(
                R_sub=R,
                sub_read_ids=dtw_read_ids,
                full_read_ids=read_ids,
                dist_read_to_tx_full=dist_read_to_tx_full,
                sigma=sigma_use,
                prior_weights=prior_weights,
            )
            quant_read_ids = read_ids

        quant_results = quantify_transcripts(
            R_quant, hard_quant, candidate_set.candidates, quant_read_ids
        )
        populate_quant_scores(quant_results, composite_scores)

        # T8: populate max_R per transcript (max responsibility across all reads)
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R_quant[:, j].max()) if R_quant.shape[0] > 0 else 0.0

        # T8: persist R-matrix to disk for FP-by-EM analysis
        if self.config.persist_R_matrix:
            np.save(str(work_dir / "R.npy"), R_quant.astype(np.float32))
            meta = {
                "read_ids": list(quant_read_ids),
                "candidate_ids": list(candidate_ids),
                "interval_region": interval.region_string,
                "sigma_used": sigma_use,
                "was_subsampled": len(dtw_read_ids) != len(read_ids),
                "n_reads_subsampled": len(dtw_read_ids),
                "n_reads_full": len(read_ids),
            }
            with open(work_dir / "R_meta.json", "w") as _f:
                json.dump(meta, _f)

        return quant_results

    def _process_interval_mappy_argmax(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Optional[List[QuantResult]]:
        """R1 ablation: assign reads via mappy argmax (no signal, no EM).

        Returns a list of QuantResult with abundance = read count, or None if
        no candidates have sequences.
        """
        from fin.ablation.mappy_argmax import (
            mappy_argmax_assignment,
            per_tx_counts_from_argmax,
        )

        candidates = candidate_set.candidates
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [
            (rid, read_sequences.get(rid, "")) for rid in read_ids
        ]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]

        assignment = mappy_argmax_assignment(reads_iter, list(candidates))
        counts = per_tx_counts_from_argmax(assignment, list(candidates))

        quant_results: List[QuantResult] = []
        for cand in candidates:
            cnt = counts.get(cand.candidate_id, 0.0)
            assigned = [
                rid for rid, cid in assignment.items()
                if cid == cand.candidate_id
            ]
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=cnt,
                confidence=1.0 if cnt > 0 else 0.0,
                num_assigned_reads=len(assigned),
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                assigned_read_ids=tuple(assigned),
            )
            quant_results.append(qr)

        logger.info(
            "R1 mappy argmax interval %s: %d reads -> %d candidates",
            interval.region_string,
            len(read_ids),
            len(candidates),
        )
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
        from fin.io.io_fasta import FASTAReader

        seqs = {}
        with FASTAReader(path) as reader:
            for record in reader.iterate_records():
                seqs[record.id] = record.sequence
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
