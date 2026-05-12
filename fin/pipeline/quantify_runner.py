"""Multi-sample quantification runner using GTF-only candidates."""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    aggregate_across_intervals,
    compute_tpm,
    quantify_transcripts,
)
from fin.candidates.dataclasses import TranscriptCandidate
from fin.candidates.discovery import discover_gtf_only
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.scoring.eventalign_parser import (
    build_distance_matrix,
    parse_eventalign_tsv,
)
from fin.scoring.composite import (
    derive_prior_weights,
    populate_quant_scores,
    score_candidates_composite,
    subsample_reads_for_dtw,
)
from fin.scoring.external_tools import ExternalToolPaths, ExternalToolRunner
from fin.scoring.signal_dtw import compute_read_to_read_dtw, extract_signal_segments

logger = logging.getLogger(__name__)


@dataclass
class SampleInput:
    """Input files for a single sample."""

    name: str
    bam_path: str
    fastq_path: str
    signal_path: str


@dataclass
class IntervalAssignment:
    """Per-interval read-to-candidate assignment data for BAM output."""

    read_ids: List[str]
    hard_assignments: np.ndarray
    candidates: List[TranscriptCandidate]
    work_dir: Path  # interval work dir where per-candidate BAMs live


def intervals_from_gtf(gtf_reader) -> List[GenomicInterval]:
    """Generate merged genomic intervals from GTF gene boundaries.

    Groups genes by chromosome, merges overlapping spans, and returns
    one GenomicInterval per merged region.

    Args:
        gtf_reader: Parsed GTFReader instance.

    Returns:
        List of merged GenomicInterval objects.
    """
    # Collect (chrom, start, end, strand) from all genes
    by_chrom: Dict[str, List[Tuple[int, int, str]]] = {}
    for gene in gtf_reader.iterate_genes():
        if gene.start >= gene.end:
            continue
        by_chrom.setdefault(gene.chrom, []).append(
            (gene.start, gene.end, gene.strand)
        )

    intervals = []
    for chrom, spans in sorted(by_chrom.items()):
        # Sort by start position
        spans.sort(key=lambda x: x[0])

        # Merge overlapping spans
        merged_start, merged_end, merged_strand = spans[0]
        for start, end, strand in spans[1:]:
            if start <= merged_end:
                # Overlapping — extend
                merged_end = max(merged_end, end)
            else:
                # Gap — emit previous, start new
                intervals.append(
                    GenomicInterval(
                        chrom=chrom,
                        start=merged_start,
                        end=merged_end,
                        strand=merged_strand,
                    )
                )
                merged_start, merged_end, merged_strand = start, end, strand

        # Emit last interval
        intervals.append(
            GenomicInterval(
                chrom=chrom,
                start=merged_start,
                end=merged_end,
                strand=merged_strand,
            )
        )

    logger.info("Generated %d intervals from GTF genes", len(intervals))
    return intervals


class QuantifyRunner:
    """Orchestrates multi-sample quantification using GTF-only candidates."""

    def __init__(
        self,
        gtf_path: str,
        genome_fasta_path: str,
        samples: List[SampleInput],
        output_dir: str,
        signal_format: str = "slow5",
        use_gpu: bool = True,
        f5c_path: str = "f5c",
        em_sigma: float = 1.0,
        em_beta: float = 0.5,
        em_max_iter: int = 1000,
        em_tol: float = 1e-4,
        signal_normalize: bool = True,
        config: Optional[PipelineConfig] = None,
    ):
        self.gtf_path = gtf_path
        self.genome_fasta_path = genome_fasta_path
        self.samples = samples
        self.output_dir = Path(output_dir)
        self.config = config

        # When config is provided, config wins for all overlapping fields.
        # This keeps a single source of truth for scoring/EM/DTW parameters.
        if config is not None:
            self.signal_format = config.signal_format
            self.use_gpu = config.use_gpu
            self.f5c_path = config.f5c_path
            self.em_sigma = config.em_sigma
            self.em_beta = config.em_beta
            self.em_max_iter = config.em_max_iter
            self.em_tol = config.em_tol
            self.signal_normalize = config.signal_normalize
        else:
            self.signal_format = signal_format
            self.use_gpu = use_gpu
            self.f5c_path = f5c_path
            self.em_sigma = em_sigma
            self.em_beta = em_beta
            self.em_max_iter = em_max_iter
            self.em_tol = em_tol
            self.signal_normalize = signal_normalize

        self._gtf_reader = None
        self._genome_fasta: Dict[str, str] = {}
        self._intervals: List[GenomicInterval] = []

    def setup(self):
        """Load GTF and genome FASTA once, validate tools, generate intervals."""
        # Validate external tools
        tool_paths = ExternalToolPaths(f5c=self.f5c_path)
        missing = tool_paths.validate()
        if missing:
            raise RuntimeError(f"Missing external tools: {', '.join(missing)}")

        # Load GTF
        from fin.io.io_gtf import GTFReader

        self._gtf_reader = GTFReader(self.gtf_path)
        self._gtf_reader.open()
        self._gtf_reader.parse()

        # Load genome FASTA
        from fin.io.io_fasta import FASTAReader

        with FASTAReader(self.genome_fasta_path) as reader:
            for record in reader.iterate_records():
                self._genome_fasta[record.id] = record.sequence

        # Generate intervals from GTF gene boundaries
        self._intervals = intervals_from_gtf(self._gtf_reader)

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(
            "QuantifyRunner setup: %d genes, %d transcripts, %d intervals, %d samples",
            len(self._gtf_reader.genes),
            len(self._gtf_reader.transcripts),
            len(self._intervals),
            len(self.samples),
        )

    def run(self):
        """Run quantification for all samples and write output matrices."""
        # Build transcript length map from GTF (spliced exon length)
        transcript_lengths: Dict[str, int] = {}
        for tx in self._gtf_reader.iterate_transcripts():
            transcript_lengths[tx.transcript_id] = tx.exon_length

        # Run each sample
        all_sample_results: Dict[str, Dict[str, QuantResult]] = {}
        for sample in self.samples:
            logger.info("Processing sample: %s", sample.name)
            results = self._run_sample(sample)
            all_sample_results[sample.name] = results

        # Collect all transcript IDs across all samples
        all_tx_ids = sorted(
            {tid for res in all_sample_results.values() for tid in res}
        )

        # Build count matrix and TPM matrix
        sample_names = [s.name for s in self.samples]

        counts_rows = []
        tpm_rows = []
        for sample_name in sample_names:
            results = all_sample_results[sample_name]
            tpm = compute_tpm(results, transcript_lengths)

            counts_row = []
            tpm_row = []
            for tid in all_tx_ids:
                qr = results.get(tid)
                counts_row.append(qr.abundance if qr else 0.0)
                tpm_row.append(tpm.get(tid, 0.0))
            counts_rows.append(counts_row)
            tpm_rows.append(tpm_row)

        # Write TSV files
        self._write_matrix(
            self.output_dir / "counts.tsv", all_tx_ids, sample_names, counts_rows
        )
        self._write_matrix(
            self.output_dir / "tpm.tsv", all_tx_ids, sample_names, tpm_rows
        )

        logger.info(
            "Quantification complete: %d transcripts x %d samples",
            len(all_tx_ids),
            len(sample_names),
        )

    def _run_sample(self, sample: SampleInput) -> Dict[str, QuantResult]:
        """Run the scoring pipeline for a single sample across all intervals."""
        sample_work_dir = self.output_dir / "work" / sample.name
        sample_work_dir.mkdir(parents=True, exist_ok=True)

        # Create tool runner and build f5c index for this sample
        tool_paths = ExternalToolPaths(f5c=self.f5c_path)
        tool_runner = ExternalToolRunner(
            fastq_path=sample.fastq_path,
            signal_path=sample.signal_path,
            signal_format=self.signal_format,
            work_dir=str(sample_work_dir),
            tools=tool_paths,
        )
        tool_runner.build_f5c_index()

        # Open signal reader for this sample
        signal_reader = self._open_signal_reader(sample.signal_path)

        all_quant_results: List[List[QuantResult]] = []
        all_assignments: List[IntervalAssignment] = []
        for i, interval in enumerate(self._intervals):
            logger.info(
                "Sample %s: interval %d/%d: %s",
                sample.name,
                i + 1,
                len(self._intervals),
                interval.region_string,
            )
            result = self._process_interval(
                interval, sample, tool_runner, signal_reader, sample_work_dir
            )
            if result:
                quant, assignment = result
                all_quant_results.append(quant)
                all_assignments.append(assignment)

        signal_reader.close()

        # Build per-sample assignment BAM
        if all_assignments:
            bam_path = self.output_dir / f"{sample.name}.bam"
            self._build_assignment_bam(sample, all_assignments, bam_path)

        # Aggregate across intervals
        aggregated = aggregate_across_intervals(all_quant_results)

        # Resolve gene_ids from GTF
        for cid, qr in aggregated.items():
            tx = self._gtf_reader.get_transcript(cid)
            if tx:
                qr.gene_id = tx.gene_id
            if not qr.gene_id:
                qr.gene_id = qr.candidate_id

        return aggregated

    def _process_interval(
        self,
        interval: GenomicInterval,
        sample: SampleInput,
        tool_runner: ExternalToolRunner,
        signal_reader,
        sample_work_dir: Path,
    ) -> Optional[Tuple[List[QuantResult], IntervalAssignment]]:
        """Process a single interval for a single sample.

        Returns a tuple of (quantification results, interval assignment data)
        or None if no candidates/reads found.
        """
        work_dir = sample_work_dir / interval.region_string.replace(":", "_").replace("-", "_")
        work_dir.mkdir(parents=True, exist_ok=True)

        chrom_seq = self._genome_fasta.get(interval.chrom, "")

        # GTF-only candidate discovery
        candidate_set = discover_gtf_only(
            interval=interval,
            bam_path=sample.bam_path,
            gtf_reader=self._gtf_reader,
            genome_fasta=chrom_seq,
        )

        if candidate_set.num_candidates == 0 or not candidate_set.read_ids:
            return None

        read_ids = sorted(candidate_set.read_ids)
        if not read_ids:
            return None
        candidate_ids = candidate_set.candidate_ids()

        # Per-candidate scoring (mappy + f5c eventalign)
        tsv_paths = tool_runner.score_candidates(candidate_set, work_dir)

        candidate_lengths = {
            c.candidate_id: len(c.sequence) for c in candidate_set.candidates
        }
        all_scores = []
        for tsv_path in tsv_paths:
            scores = parse_eventalign_tsv(str(tsv_path), candidate_lengths)
            all_scores.extend(scores)

        # DTW subsampling: cap reads for the O(n^2) read-to-read DTW.
        max_dtw = (
            self.config.max_reads_per_interval_for_dtw
            if self.config is not None
            else 2000
        )
        dtw_read_ids = subsample_reads_for_dtw(read_ids, max_dtw)
        if len(dtw_read_ids) != len(read_ids):
            logger.info(
                "Interval %s: subsampling DTW reads %d -> %d",
                interval.region_string,
                len(read_ids),
                len(dtw_read_ids),
            )

        dist_read_to_tx = build_distance_matrix(
            all_scores, dtw_read_ids, candidate_ids
        )

        # Signal DTW (read-to-read)
        segments = extract_signal_segments(
            all_scores,
            signal_reader,
            signal_format=self.signal_format,
            normalize=self.signal_normalize,
        )
        dist_read_to_read = compute_read_to_read_dtw(
            segments, dtw_read_ids, use_gpu=self.use_gpu
        )

        # Shape checks: both distance matrices share the same read axis.
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

        # --- Composite scoring (uniform R seed) ---
        score_alpha = self.config.score_alpha if self.config is not None else 0.5
        R_uniform = np.full((n_reads, n_tx), 1.0 / max(n_tx, 1))
        composite_scores = score_candidates_composite(
            candidates=candidate_set.candidates,
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            R=R_uniform,
            alpha=score_alpha,
            use_gpu=self.use_gpu,
        )
        combined_scores_arr = np.array(
            [s.combined for s in composite_scores], dtype=float
        )

        # --- Determine prior_weights from composite scores ---
        use_prior = self.config.use_prior if self.config is not None else True
        prior_cap = (
            self.config.prior_weight_cap if self.config is not None else 10.0
        )
        prior_weights: Optional[np.ndarray] = None
        if use_prior:
            prior_weights = derive_prior_weights(
                combined_scores_arr, n_tx, prior_cap
            )

        # EM assignment
        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=self.em_sigma,
            beta=self.em_beta,
            max_iter=self.em_max_iter,
            tol=self.em_tol,
            verbose=False,
            use_gpu=self.use_gpu,
            prior_weights=prior_weights,
        )

        # Quantification + score field population
        quant = quantify_transcripts(
            R, hard_assignments, candidate_set.candidates, dtw_read_ids
        )
        populate_quant_scores(quant, composite_scores)

        # T8: populate max_R per transcript
        for j, qr in enumerate(quant):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0

        # T8: persist R-matrix to disk for FP-by-EM analysis
        persist = (
            self.config.persist_R_matrix if self.config is not None else True
        )
        if persist:
            np.save(str(work_dir / "R.npy"), R.astype(np.float32))
            meta = {
                "read_ids": list(dtw_read_ids),
                "candidate_ids": list(candidate_ids),
                "interval_region": interval.region_string,
                "sigma_used": self.em_sigma,
                "was_subsampled": len(dtw_read_ids) != len(read_ids),
                "n_reads_subsampled": len(dtw_read_ids),
                "n_reads_full": len(read_ids),
            }
            with open(work_dir / "R_meta.json", "w") as _f:
                json.dump(meta, _f)

        assignment = IntervalAssignment(
            read_ids=dtw_read_ids,
            hard_assignments=hard_assignments,
            candidates=list(candidate_set.candidates),
            work_dir=work_dir,
        )

        return quant, assignment

    def _build_assignment_bam(
        self,
        sample: SampleInput,
        interval_assignments: List[IntervalAssignment],
        output_path: Path,
    ) -> None:
        """Build a single sorted+indexed BAM with hard-assigned reads per sample.

        Merges reads from per-candidate BAMs, keeping only reads that were
        hard-assigned to each candidate by EM. Reference IDs are remapped
        to a unified header containing all candidate sequences.
        """
        import pysam

        # Collect all unique candidates across intervals → unified header
        seen = {}  # candidate_id → TranscriptCandidate
        for ia in interval_assignments:
            for cand in ia.candidates:
                if cand.candidate_id not in seen:
                    seen[cand.candidate_id] = cand

        # Stable ordering for reference IDs
        all_candidates = sorted(seen.values(), key=lambda c: c.candidate_id)
        cand_to_ref_id = {c.candidate_id: i for i, c in enumerate(all_candidates)}

        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [
                    {"SN": c.candidate_id, "LN": len(c.sequence)}
                    for c in all_candidates
                ],
            }
        )

        unsorted_path = str(output_path) + ".unsorted.bam"
        with pysam.AlignmentFile(unsorted_path, "wb", header=header) as outf:
            for ia in interval_assignments:
                # Group read indices by assigned candidate
                reads_by_candidate: Dict[int, List[str]] = {}
                for read_idx, cand_idx in enumerate(ia.hard_assignments):
                    if cand_idx < 0 or cand_idx >= len(ia.candidates):
                        continue  # unassigned
                    reads_by_candidate.setdefault(int(cand_idx), []).append(
                        ia.read_ids[read_idx]
                    )

                for cand_idx, assigned_reads in reads_by_candidate.items():
                    candidate = ia.candidates[cand_idx]
                    bam_path = ia.work_dir / candidate.candidate_id / "aligned.bam"
                    if not bam_path.exists():
                        logger.warning(
                            "Per-candidate BAM not found: %s", bam_path
                        )
                        continue

                    assigned_set = set(assigned_reads)
                    new_ref_id = cand_to_ref_id[candidate.candidate_id]

                    with pysam.AlignmentFile(str(bam_path), "rb") as inf:
                        for aln in inf:
                            if aln.query_name in assigned_set:
                                a = pysam.AlignedSegment(header)
                                a.query_name = aln.query_name
                                a.query_sequence = aln.query_sequence
                                a.flag = aln.flag
                                a.reference_id = new_ref_id
                                a.reference_start = aln.reference_start
                                a.mapping_quality = aln.mapping_quality
                                a.cigar = aln.cigar
                                if aln.query_qualities is not None:
                                    a.query_qualities = aln.query_qualities
                                outf.write(a)

        # Sort and index
        pysam.sort("-o", str(output_path), unsorted_path)
        os.remove(unsorted_path)
        pysam.index(str(output_path))

        logger.info(
            "Wrote assignment BAM for sample %s: %s", sample.name, output_path
        )

    def _open_signal_reader(self, signal_path: str):
        """Open the appropriate signal reader."""
        if self.signal_format == "pod5":
            from fin.io.io_pod5 import Pod5Reader

            reader = Pod5Reader(signal_path)
            reader.open()
            return reader
        else:
            from fin.io.io_slow5 import Slow5Reader

            reader = Slow5Reader(signal_path)
            reader.open()
            return reader

    def cleanup(self):
        """Close file handles."""
        if self._gtf_reader:
            self._gtf_reader.close()

    @staticmethod
    def _write_matrix(
        path: Path,
        transcript_ids: List[str],
        sample_names: List[str],
        rows: List[List[float]],
    ) -> None:
        """Write a TSV matrix (transcripts as rows, samples as columns)."""
        with open(path, "w") as f:
            f.write("transcript_id\t" + "\t".join(sample_names) + "\n")
            for i, tid in enumerate(transcript_ids):
                values = [f"{rows[s][i]:.4f}" for s in range(len(sample_names))]
                f.write(tid + "\t" + "\t".join(values) + "\n")

        logger.info("Wrote %s (%d transcripts x %d samples)", path, len(transcript_ids), len(sample_names))
