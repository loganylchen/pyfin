#!/usr/bin/env python3
"""
Region-based Transcript Analysis Workflow

This example demonstrates a complete workflow for analyzing nanopore Direct RNA-seq
data at the transcript level within isolated genomic regions.

Workflow:
1. Load input files: BAM (reads mapped to genome), genome reference, 
   transcriptome reference, GTF annotation, and POD5 signal data
2. Separate reads into isolated genomic regions
3. For each region:
   a. Get all candidate transcript sequences from the region
   b. For each read in the region, perform eventalign to all candidate transcripts
   c. Compute DTW pairwise distances among all reads in the region
4. Output analysis results

Usage:
    python region_transcript_analysis_workflow.py \\
        --bam reads.bam \\
        --genome genome.fa \\
        --transcriptome transcripts.fa \\
        --gtf annotation.gtf \\
        --pod5 signals.pod5 \\
        --output results/

Requirements:
    - pysam (for BAM handling)
    - pod5 (for signal reading)
    - CUDA toolkit (optional, for GPU-accelerated DTW)
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict
import json

import numpy as np

# Optional visualization imports
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Optional clustering imports
try:
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import squareform

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Add parent directory to path for development
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

# Import fin modules
from fin.io import (
    generate_isolated_intervals,
    extract_reads_for_interval,
    FASTAReader,
    BamReader,
    GTFReader,
)
from fin.io.io_pod5 import Pod5Reader

# Import eventalign functions
try:
    from fin._f5c import detect_events, eventalign, profile_hmm_eventalign, is_available

    EVENTALIGN_AVAILABLE = is_available()
except ImportError as e:
    EVENTALIGN_AVAILABLE = False
    print(f"Warning: eventalign not available: {e}")

# Import DTW functions
try:
    from fin._dtw import dtw_pairwise, is_available as dtw_is_available

    DTW_AVAILABLE = dtw_is_available()
except ImportError as e:
    DTW_AVAILABLE = False
    print(f"Warning: DTW not available: {e}")

# Import EM assignment function
try:
    from fin.analysis.assignments import em_with_coherence

    ASSIGNMENT_AVAILABLE = True
except ImportError as e:
    ASSIGNMENT_AVAILABLE = False
    print(f"Warning: EM assignment not available: {e}")

from fin.utils.log_config import get_package_logger

logger = get_package_logger(__name__)


class RegionTranscriptAnalyzer:
    """
    Analyzer for performing transcript-level analysis within genomic regions.

    This class coordinates the workflow of:
    1. Separating reads into isolated regions
    2. Getting candidate transcript sequences per region
    3. Performing eventalign for each read to each candidate transcript
    4. Computing DTW distances among reads in each region
    """

    def __init__(
        self,
        bam_path: str,
        genome_path: str,
        transcriptome_path: str,
        gtf_path: str,
        pod5_path: str,
        output_dir: str,
        kmer_size: int = 5,
    ):
        """
        Initialize the analyzer.

        RNA-only: Events are automatically reversed internally to match the
        5'→3' sequence direction. You do NOT need to reverse sequences.

        Args:
            bam_path: Path to BAM file (reads mapped to genome)
            genome_path: Path to genome FASTA reference
            transcriptome_path: Path to transcriptome FASTA reference
            gtf_path: Path to GTF annotation file
            pod5_path: Path to POD5 signal file
            output_dir: Directory for output files
            kmer_size: K-mer size for eventalign (5 or 9, default: 5)
        """
        self.bam_path = Path(bam_path)
        self.genome_path = Path(genome_path)
        self.transcriptome_path = Path(transcriptome_path)
        self.gtf_path = Path(gtf_path)
        self.pod5_path = Path(pod5_path)
        self.output_dir = Path(output_dir)
        self.kmer_size = kmer_size

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Data caches
        self.genome_sequences: Dict[str, str] = {}
        self.transcriptome_sequences: Dict[str, str] = {}
        self.transcript_to_gene: Dict[str, str] = {}
        self.gene_to_transcripts: Dict[str, List[str]] = defaultdict(list)

        # Results storage
        self.regions: List[Any] = []
        self.region_results: Dict[str, Dict] = {}

        logger.info(f"Initialized RegionTranscriptAnalyzer")
        logger.info(f"  BAM: {self.bam_path}")
        logger.info(f"  Genome: {self.genome_path}")
        logger.info(f"  Transcriptome: {self.transcriptome_path}")
        logger.info(f"  GTF: {self.gtf_path}")
        logger.info(f"  POD5: {self.pod5_path}")
        logger.info(f"  Output: {self.output_dir}")

    def load_references(self):
        """Load genome and transcriptome reference sequences."""
        logger.info("Loading reference sequences...")

        # Load genome
        with FASTAReader(str(self.genome_path)) as reader:
            for record in reader.iterate_records():
                self.genome_sequences[record.id] = record.sequence
        logger.info(f"  Loaded {len(self.genome_sequences)} genome sequences")

        # Load transcriptome
        with FASTAReader(str(self.transcriptome_path)) as reader:
            for record in reader.iterate_records():
                self.transcriptome_sequences[record.id] = record.sequence
        logger.info(f"  Loaded {len(self.transcriptome_sequences)} transcript sequences")

    def load_annotation(self):
        """Load GTF annotation and build transcript-gene mappings."""
        logger.info("Loading GTF annotation...")

        with GTFReader(str(self.gtf_path)) as reader:
            reader.parse()

            for tx in reader.iterate_transcripts():
                self.transcript_to_gene[tx.transcript_id] = tx.gene_id
                self.gene_to_transcripts[tx.gene_id].append(tx.transcript_id)

        logger.info(f"  Loaded {len(self.transcript_to_gene)} transcripts")
        logger.info(f"  Across {len(self.gene_to_transcripts)} genes")

    def generate_regions(self) -> List[Any]:
        """
        Generate isolated genomic regions from BAM and GTF.

        Returns:
            List of GenomicInterval objects
        """
        logger.info("Generating isolated regions...")

        result = generate_isolated_intervals(str(self.bam_path), str(self.gtf_path))

        self.regions = result["intervals"]
        fusion_ids = result["fusion_read_ids"]

        logger.info(f"  Generated {len(self.regions)} isolated regions")
        logger.info(f"  Identified {len(fusion_ids)} fusion candidate reads")

        return self.regions

    def get_candidate_transcripts_for_region(
        self, chrom: str, start: int, end: int, strand: Optional[str] = None
    ) -> List[Tuple[str, str]]:
        """
        Get all candidate transcript sequences that overlap a region.

        Args:
            chrom: Chromosome name
            start: Region start (0-based)
            end: Region end
            strand: Strand to filter ('+', '-', or None for both)

        Returns:
            List of (transcript_id, sequence) tuples
        """
        candidates = []

        with GTFReader(str(self.gtf_path)) as reader:
            reader.parse()

            # Get transcripts overlapping this region
            transcripts = reader.get_transcripts_in_region(chrom, start, end)

            for tx in transcripts:
                # Filter by strand if specified
                if strand is not None and tx.strand != strand:
                    continue

                tx_id = tx.transcript_id
                if tx_id in self.transcriptome_sequences:
                    candidates.append((tx_id, self.transcriptome_sequences[tx_id]))

        return candidates

    def get_reads_for_region(
        self, chrom: str, start: int, end: int, strand: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Get all reads mapping to a region, optionally filtered by strand.

        Args:
            chrom: Chromosome name
            start: Region start (0-based)
            end: Region end
            strand: Strand to filter ('+', '-', or None for both)

        Returns:
            List of read dictionaries
        """
        reads = []

        with BamReader(str(self.bam_path)) as reader:
            for read in reader.fetch(region=f"{chrom}:{start}-{end}"):
                # Filter by strand if specified
                if strand is not None:
                    # BAM: is_reverse=False means forward (+), is_reverse=True means reverse (-)
                    read_strand = "-" if (hasattr(read, "is_reverse") and read.is_reverse) else "+"
                    if read_strand != strand:
                        continue
                reads.append(read)

        return reads

    def get_signal_for_read(self, read_id: str, pod5_reader: Pod5Reader) -> Optional[np.ndarray]:
        """
        Get calibrated signal for a read from POD5.

        Args:
            read_id: Read ID
            pod5_reader: Open Pod5Reader instance

        Returns:
            Calibrated signal as numpy array, or None if not found
        """
        result = pod5_reader.get_calibrated_signal(read_id)
        if result is None:
            return None

        signal, metadata = result
        return np.array(signal, dtype=np.float32)

    def eventalign_read_to_transcript(
        self, signal: np.ndarray, transcript_seq: str, use_profile_hmm: bool = True
    ) -> Optional[Dict[str, Any]]:
        """
        Perform eventalign of a read's signal to a transcript sequence.

        For Direct RNA Sequencing (DRS), the RNA strand passes through the pore
        in the 3'→5' direction, so we need to reverse the transcript sequence
        (which is stored 5'→3') before alignment.

        Args:
            signal: Calibrated nanopore signal
            transcript_seq: Transcript sequence (DNA/RNA, stored 5'→3')
            use_profile_hmm: Use Profile HMM (True) or simple ABEA (False)

        Returns:
            Alignment result dictionary, or None if failed
        """
        if not EVENTALIGN_AVAILABLE:
            logger.warning("Eventalign not available")
            return None

        try:
            # RNA-only mode: events are automatically reversed internally
            # Use the original 5'→3' sequence directly - no reversal needed

            if use_profile_hmm:
                result = profile_hmm_eventalign(
                    raw_signal=signal,
                    sequence=transcript_seq,
                    kmer_size=self.kmer_size,
                )
            else:
                result = eventalign(
                    raw_signal=signal,
                    sequence=transcript_seq,
                    kmer_size=self.kmer_size,
                )

            # Store sequence info for reference
            result["sequence"] = transcript_seq

            return result
        except Exception as e:
            logger.warning(f"Eventalign failed: {e}")
            return None

    def compute_dtw_distances(
        self, signals: List[np.ndarray], normalize_length: bool = True, target_length: int = 1000
    ) -> Optional[np.ndarray]:
        """
        Compute pairwise DTW distances among signals.

        Args:
            signals: List of signal arrays
            normalize_length: Whether to normalize signals to same length
            target_length: Target length for normalization

        Returns:
            Distance matrix (n x n), or None if DTW not available
        """
        if not DTW_AVAILABLE:
            logger.warning("DTW not available (requires CUDA)")
            return None

        if len(signals) < 2:
            logger.warning("Need at least 2 signals for pairwise DTW")
            return None

        # Normalize signals to same length if needed
        if normalize_length:
            normalized = []
            for sig in signals:
                if len(sig) == target_length:
                    normalized.append(sig)
                else:
                    # Resample to target length
                    indices = np.linspace(0, len(sig) - 1, target_length)
                    resampled = np.interp(indices, np.arange(len(sig)), sig)
                    normalized.append(resampled.astype(np.float32))
            signals = normalized

        # Stack into 2D array
        signal_matrix = np.stack(signals, axis=0)

        try:
            distance_matrix = dtw_pairwise(signal_matrix)
            return distance_matrix
        except Exception as e:
            logger.warning(f"DTW computation failed: {e}")
            return None

    def cluster_reads_by_dtw(
        self,
        distance_matrix: np.ndarray,
        read_ids: List[str],
        method: str = "ward",
        n_clusters: Optional[int] = None,
        distance_threshold: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Perform hierarchical clustering on reads based on DTW distance matrix.

        Args:
            distance_matrix: Pairwise DTW distance matrix (n x n)
            read_ids: List of read IDs corresponding to matrix rows/columns
            method: Linkage method ('ward', 'complete', 'average', 'single')
            n_clusters: Number of clusters to form (optional)
            distance_threshold: Distance threshold for forming clusters (optional)

        Returns:
            Dictionary with clustering results
        """
        if not SCIPY_AVAILABLE:
            logger.warning("scipy not available for clustering")
            return {"error": "scipy not installed"}

        if len(read_ids) < 2:
            return {"error": "Need at least 2 reads for clustering"}

        try:
            # Convert distance matrix to condensed form for scipy
            # Ensure matrix is symmetric and has zero diagonal
            dist_matrix = np.array(distance_matrix)
            np.fill_diagonal(dist_matrix, 0)
            dist_matrix = (dist_matrix + dist_matrix.T) / 2  # Ensure symmetry

            condensed_dist = squareform(dist_matrix)

            # Perform hierarchical clustering
            linkage_matrix = linkage(condensed_dist, method=method)

            # Determine cluster assignments
            if n_clusters is not None:
                cluster_labels = fcluster(linkage_matrix, n_clusters, criterion="maxclust")
            elif distance_threshold is not None:
                cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion="distance")
            else:
                # Default: use 2 clusters
                cluster_labels = fcluster(linkage_matrix, 2, criterion="maxclust")

            # Build cluster membership
            clusters = defaultdict(list)
            for read_id, label in zip(read_ids, cluster_labels):
                clusters[int(label)].append(read_id)

            return {
                "method": method,
                "n_clusters": len(clusters),
                "cluster_labels": cluster_labels.tolist(),
                "clusters": dict(clusters),
                "linkage_matrix": linkage_matrix.tolist(),
                "read_ids": read_ids,
            }

        except Exception as e:
            logger.error(f"Clustering failed: {e}")
            return {"error": str(e)}

    def plot_dtw_heatmap(
        self,
        distance_matrix: np.ndarray,
        read_ids: List[str],
        output_path: Path,
        title: str = "DTW Distance Matrix",
        cluster_result: Optional[Dict] = None,
        figsize: Tuple[int, int] = (12, 10),
    ) -> Optional[Path]:
        """
        Generate a clustered heatmap visualization of the DTW distance matrix.

        Creates a heatmap with hierarchical clustering dendrograms on both
        rows and columns, reordering the matrix to group similar reads together.

        Args:
            distance_matrix: Pairwise DTW distance matrix (symmetric)
            read_ids: List of read IDs
            output_path: Path to save the heatmap image
            title: Plot title
            cluster_result: Optional clustering result (will compute if not provided)
            figsize: Figure size

        Returns:
            Path to saved image, or None if visualization failed
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for visualization")
            return None

        try:
            dist_matrix = np.array(distance_matrix)
            n_reads = len(read_ids)

            # Create clustered heatmap with dendrograms on both axes
            if SCIPY_AVAILABLE and n_reads >= 2:
                # Get or compute linkage matrix
                if cluster_result and "linkage_matrix" in cluster_result:
                    linkage_matrix = np.array(cluster_result["linkage_matrix"])
                else:
                    # Compute linkage from distance matrix
                    from scipy.spatial.distance import squareform

                    condensed_dist = squareform(dist_matrix)
                    linkage_matrix = linkage(condensed_dist, method="average")

                # Create figure with grid for dendrograms + heatmap + colorbar
                # Layout:
                #   [empty]     [col_dendro]  [empty]
                #   [row_dendro] [heatmap]    [colorbar]
                fig = plt.figure(figsize=(figsize[0] + 3, figsize[1] + 2))

                # Grid ratios: left dendrogram, main heatmap, colorbar
                gs = fig.add_gridspec(
                    2,
                    3,
                    width_ratios=[1.5, 6, 0.3],
                    height_ratios=[1.5, 6],
                    wspace=0.02,
                    hspace=0.02,
                )

                # Top dendrogram (columns)
                ax_col_dendro = fig.add_subplot(gs[0, 1])
                col_dendro = dendrogram(
                    linkage_matrix,
                    orientation="top",
                    no_labels=True,
                    ax=ax_col_dendro,
                    color_threshold=0,
                    above_threshold_color="C0",
                )
                ax_col_dendro.set_xticks([])
                ax_col_dendro.set_yticks([])
                ax_col_dendro.spines["top"].set_visible(False)
                ax_col_dendro.spines["right"].set_visible(False)
                ax_col_dendro.spines["bottom"].set_visible(False)
                ax_col_dendro.spines["left"].set_visible(False)

                # Left dendrogram (rows)
                ax_row_dendro = fig.add_subplot(gs[1, 0])
                row_dendro = dendrogram(
                    linkage_matrix,
                    orientation="left",
                    no_labels=True,
                    ax=ax_row_dendro,
                    color_threshold=0,
                    above_threshold_color="C0",
                )
                ax_row_dendro.set_xticks([])
                ax_row_dendro.set_yticks([])
                ax_row_dendro.spines["top"].set_visible(False)
                ax_row_dendro.spines["right"].set_visible(False)
                ax_row_dendro.spines["bottom"].set_visible(False)
                ax_row_dendro.spines["left"].set_visible(False)

                # Get the ordering from dendrogram (same for both since matrix is symmetric)
                order = row_dendro["leaves"]
                ordered_matrix = dist_matrix[order][:, order]
                ordered_ids = [read_ids[i] for i in order]

                # Main heatmap
                ax_heat = fig.add_subplot(gs[1, 1])
                im = ax_heat.imshow(ordered_matrix, cmap="viridis", aspect="auto")

                # Labels
                if n_reads <= 30:
                    ax_heat.set_xticks(range(n_reads))
                    ax_heat.set_xticklabels(ordered_ids, rotation=90, fontsize=8)
                    ax_heat.set_yticks(range(n_reads))
                    ax_heat.set_yticklabels(ordered_ids, fontsize=8)
                else:
                    ax_heat.set_xticks([])
                    ax_heat.set_yticks([])
                    ax_heat.set_xlabel(f"Reads (n={n_reads})")
                    ax_heat.set_ylabel(f"Reads (n={n_reads})")

                # Colorbar
                ax_cbar = fig.add_subplot(gs[1, 2])
                cbar = plt.colorbar(im, cax=ax_cbar)
                cbar.set_label("DTW Distance")

                # Title at the top
                fig.suptitle(title, fontsize=14, fontweight="bold", y=0.98)

            else:
                # Simple heatmap without dendrogram (scipy not available or < 2 reads)
                fig, ax = plt.subplots(figsize=figsize)
                im = ax.imshow(dist_matrix, cmap="viridis", aspect="auto")

                # Add colorbar
                cbar = plt.colorbar(im, ax=ax, shrink=0.8)
                cbar.set_label("DTW Distance")

                # Labels
                if n_reads <= 30:
                    ax.set_xticks(range(n_reads))
                    ax.set_xticklabels(read_ids, rotation=90, fontsize=8)
                    ax.set_yticks(range(n_reads))
                    ax.set_yticklabels(read_ids, fontsize=8)
                else:
                    ax.set_xlabel(f"Reads (n={n_reads})")
                    ax.set_ylabel(f"Reads (n={n_reads})")

                ax.set_title(title)

            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

            logger.info(f"  Saved clustered DTW heatmap to: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to generate heatmap: {e}")
            import traceback

            traceback.print_exc()
            return None

    def plot_assessment_heatmap(
        self,
        read_alignments: List[Dict[str, Any]],
        output_path: Path,
        metric: str = "overall_score",
        title: str = "Eventalign Quality Assessment",
        figsize: Tuple[int, int] = (14, 10),
    ) -> Optional[Path]:
        """
        Generate a heatmap visualization of eventalign quality assessments.

        Creates a heatmap showing how well each read aligns to each transcript,
        based on the quality assessment metrics.

        Args:
            read_alignments: List of read alignment results from analyze_region()
            output_path: Path to save the heatmap image
            metric: Which metric to display. Options:
                - "overall_score": Overall quality score (0-1)
                - "correlation": Signal-model correlation
                - "event_coverage": Fraction of events aligned
                - "sequence_coverage": Fraction of sequence covered
                - "match_fraction": Fraction of match states
            title: Plot title
            figsize: Figure size

        Returns:
            Path to saved image, or None if visualization failed
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for visualization")
            return None

        if not read_alignments:
            logger.warning("No read alignments to visualize")
            return None

        try:
            # Collect all unique transcripts and reads
            all_transcripts = set()
            all_reads = []
            for ra in read_alignments:
                read_id = ra.get("read_id", "unknown")
                all_reads.append(read_id)
                for aln in ra.get("alignments", []):
                    all_transcripts.add(aln.get("transcript_id", "unknown"))

            transcripts = sorted(list(all_transcripts))
            reads = all_reads

            if not transcripts or not reads:
                logger.warning("No transcripts or reads found")
                return None

            # Build the score matrix
            n_reads = len(reads)
            n_transcripts = len(transcripts)
            score_matrix = np.full((n_reads, n_transcripts), np.nan)

            # Map transcript_id to column index
            tx_to_idx = {tx: i for i, tx in enumerate(transcripts)}

            # Fill in the matrix
            for i, ra in enumerate(read_alignments):
                for aln in ra.get("alignments", []):
                    tx_id = aln.get("transcript_id")
                    if tx_id not in tx_to_idx:
                        continue

                    j = tx_to_idx[tx_id]
                    qa = aln.get("quality_assessment", {})

                    if metric == "overall_score":
                        score_matrix[i, j] = qa.get("overall_score", 0)
                    elif metric == "correlation":
                        fit = qa.get("signal_model_fit", {})
                        score_matrix[i, j] = fit.get("correlation", 0)
                    elif metric == "event_coverage":
                        cov = qa.get("coverage", {})
                        score_matrix[i, j] = cov.get("event_coverage", 0)
                    elif metric == "sequence_coverage":
                        cov = qa.get("coverage", {})
                        score_matrix[i, j] = cov.get("sequence_coverage", 0)
                    elif metric == "match_fraction":
                        hmm = qa.get("hmm_states", {})
                        score_matrix[i, j] = hmm.get("match_fraction", 0)
                    else:
                        score_matrix[i, j] = qa.get("overall_score", 0)

            # Create figure
            fig, ax = plt.subplots(figsize=figsize)

            # Mask NaN values for visualization
            masked_matrix = np.ma.masked_invalid(score_matrix)

            # Choose colormap based on metric
            if metric == "correlation":
                # Correlation can be negative
                cmap = "RdYlGn"
                vmin, vmax = -1, 1
            else:
                # Most metrics are 0-1
                cmap = "YlGnBu"
                vmin, vmax = 0, 1

            im = ax.imshow(masked_matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)

            # Colorbar
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            metric_labels = {
                "overall_score": "Overall Quality Score",
                "correlation": "Signal-Model Correlation",
                "event_coverage": "Event Coverage",
                "sequence_coverage": "Sequence Coverage",
                "match_fraction": "Match State Fraction",
            }
            cbar.set_label(metric_labels.get(metric, metric))

            # Labels
            # Truncate read IDs for display
            read_labels = [r[:15] + "..." if len(r) > 15 else r for r in reads]
            # Truncate transcript IDs for display
            tx_labels = [t[:20] + "..." if len(t) > 20 else t for t in transcripts]

            if n_reads <= 50:
                ax.set_yticks(range(n_reads))
                ax.set_yticklabels(read_labels, fontsize=7)
            else:
                ax.set_ylabel(f"Reads (n={n_reads})")

            if n_transcripts <= 30:
                ax.set_xticks(range(n_transcripts))
                ax.set_xticklabels(tx_labels, rotation=45, ha="right", fontsize=7)
            else:
                ax.set_xlabel(f"Transcripts (n={n_transcripts})")

            ax.set_title(f"{title}\nMetric: {metric_labels.get(metric, metric)}", fontsize=12)
            ax.set_xlabel("Transcript")
            ax.set_ylabel("Read")

            # Add text annotations for small matrices
            if n_reads <= 20 and n_transcripts <= 15:
                for i in range(n_reads):
                    for j in range(n_transcripts):
                        val = score_matrix[i, j]
                        if not np.isnan(val):
                            # Choose text color based on value
                            text_color = "white" if val < 0.5 else "black"
                            ax.text(
                                j,
                                i,
                                f"{val:.2f}",
                                ha="center",
                                va="center",
                                fontsize=6,
                                color=text_color,
                            )

            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

            logger.info(f"  Saved assessment heatmap to: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to generate assessment heatmap: {e}")
            import traceback

            traceback.print_exc()
            return None

    def plot_assessment_multi_metric(
        self,
        read_alignments: List[Dict[str, Any]],
        output_path: Path,
        title: str = "Eventalign Quality Assessment - Multi-Metric",
        figsize: Tuple[int, int] = (20, 12),
    ) -> Optional[Path]:
        """
        Generate a multi-panel heatmap showing multiple assessment metrics.

        Creates a figure with multiple heatmaps side by side, each showing
        a different quality metric for all read-transcript pairs.

        Args:
            read_alignments: List of read alignment results from analyze_region()
            output_path: Path to save the heatmap image
            title: Overall figure title
            figsize: Figure size

        Returns:
            Path to saved image, or None if visualization failed
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for visualization")
            return None

        if not read_alignments:
            logger.warning("No read alignments to visualize")
            return None

        try:
            # Collect all unique transcripts and reads
            all_transcripts = set()
            all_reads = []
            for ra in read_alignments:
                read_id = ra.get("read_id", "unknown")
                all_reads.append(read_id)
                for aln in ra.get("alignments", []):
                    all_transcripts.add(aln.get("transcript_id", "unknown"))

            transcripts = sorted(list(all_transcripts))
            reads = all_reads

            if not transcripts or not reads:
                logger.warning("No transcripts or reads found")
                return None

            n_reads = len(reads)
            n_transcripts = len(transcripts)
            tx_to_idx = {tx: i for i, tx in enumerate(transcripts)}

            # Define metrics to display
            metrics = [
                ("overall_score", "Overall Score", "YlGnBu", 0, 1),
                ("correlation", "Correlation", "RdYlGn", -1, 1),
                ("sequence_coverage", "Seq Coverage", "YlGnBu", 0, 1),
                ("match_fraction", "Match Fraction", "YlGnBu", 0, 1),
            ]

            # Build matrices for each metric
            matrices = {}
            for metric_key, _, _, _, _ in metrics:
                matrices[metric_key] = np.full((n_reads, n_transcripts), np.nan)

            for i, ra in enumerate(read_alignments):
                for aln in ra.get("alignments", []):
                    tx_id = aln.get("transcript_id")
                    if tx_id not in tx_to_idx:
                        continue
                    j = tx_to_idx[tx_id]
                    qa = aln.get("quality_assessment", {})

                    matrices["overall_score"][i, j] = qa.get("overall_score", 0)
                    fit = qa.get("signal_model_fit", {})
                    matrices["correlation"][i, j] = fit.get("correlation", 0)
                    cov = qa.get("coverage", {})
                    matrices["sequence_coverage"][i, j] = cov.get("sequence_coverage", 0)
                    hmm = qa.get("hmm_states", {})
                    matrices["match_fraction"][i, j] = hmm.get("match_fraction", 0)

            # Create multi-panel figure
            fig, axes = plt.subplots(2, 2, figsize=figsize)
            axes = axes.flatten()

            for idx, (metric_key, metric_label, cmap, vmin, vmax) in enumerate(metrics):
                ax = axes[idx]
                matrix = matrices[metric_key]
                masked_matrix = np.ma.masked_invalid(matrix)

                im = ax.imshow(masked_matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
                cbar = plt.colorbar(im, ax=ax, shrink=0.8)
                cbar.set_label(metric_label, fontsize=9)

                # Labels
                read_labels = [r[:12] + ".." if len(r) > 12 else r for r in reads]
                tx_labels = [t[:15] + ".." if len(t) > 15 else t for t in transcripts]

                if n_reads <= 30:
                    ax.set_yticks(range(n_reads))
                    ax.set_yticklabels(read_labels, fontsize=6)
                else:
                    ax.set_ylabel(f"Reads (n={n_reads})", fontsize=9)

                if n_transcripts <= 20:
                    ax.set_xticks(range(n_transcripts))
                    ax.set_xticklabels(tx_labels, rotation=45, ha="right", fontsize=6)
                else:
                    ax.set_xlabel(f"Transcripts (n={n_transcripts})", fontsize=9)

                ax.set_title(metric_label, fontsize=11, fontweight="bold")
                ax.set_xlabel("Transcript", fontsize=9)
                ax.set_ylabel("Read", fontsize=9)

            fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

            logger.info(f"  Saved multi-metric assessment heatmap to: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to generate multi-metric heatmap: {e}")
            import traceback

            traceback.print_exc()
            return None

    def plot_eventalign_signal(
        self,
        signal: np.ndarray,
        eventalign_result: Dict[str, Any],
        output_path: Path,
        read_id: str = "",
        transcript_id: str = "",
        max_events: int = 200,
        figsize: Tuple[int, int] = (20, 8),
    ) -> Optional[Path]:
        """
        Visualize the eventalign result for a single read.

        Creates a single panel showing raw signal with event boundaries,
        event means, and model expected means overlaid. Uses actual signal_start
        and signal_length from eventalign results.

        Args:
            signal: Raw nanopore signal
            eventalign_result: Result from profile_hmm_eventalign()
            output_path: Path to save the figure
            read_id: Read ID for title
            transcript_id: Transcript ID for title
            max_events: Maximum number of events to display (for readability)
            figsize: Figure size

        Returns:
            Path to saved image, or None if visualization failed
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for visualization")
            return None

        try:
            alignment = eventalign_result.get("alignment", [])
            scaling = eventalign_result.get("scaling", {})
            n_events = eventalign_result.get("n_events", 0)
            n_aligned = eventalign_result.get("n_aligned", 0)

            if not alignment:
                logger.warning("No alignment data to visualize")
                return None

            # Sort alignment by event_idx (raw signal order)
            alignment_sorted = sorted(
                [a for a in alignment if a.get("event_idx", -1) >= 0],
                key=lambda x: x.get("event_idx", 0),
            )

            # Limit events for visualization
            if len(alignment_sorted) > max_events:
                alignment_sorted = alignment_sorted[:max_events]

            # Create single-panel figure
            fig, ax = plt.subplots(figsize=figsize)

            # Plot raw signal
            ax.plot(signal, color="steelblue", alpha=0.5, linewidth=0.3, label="Raw signal")

            # Determine signal range to display based on aligned events
            if alignment_sorted:
                min_start = min(a.get("signal_start", 0) for a in alignment_sorted)
                max_end = max(
                    a.get("signal_start", 0) + int(a.get("signal_length", 100))
                    for a in alignment_sorted
                )
                # Add some padding
                view_start = max(0, int(min_start) - 500)
                view_end = min(len(signal), int(max_end) + 500)
            else:
                view_start = 0
                view_end = len(signal)

            # Draw event boundaries, means, and model means using signal_start/signal_length
            for i, aln in enumerate(alignment_sorted):
                signal_start = aln.get("signal_start", 0)
                signal_length = aln.get("signal_length", 0)
                event_mean = aln.get("event_mean", 0)
                model_mean = aln.get("model_mean", 0)
                scaled_model_mean = aln.get("scaled_model_mean", model_mean)
                kmer = aln.get("ref_kmer", "")
                state = aln.get("hmm_state", "M")
                event_idx = aln.get("event_idx", i)
                ref_pos = aln.get("ref_position", 0)

                start_sample = int(signal_start)
                end_sample = int(signal_start + signal_length)

                if end_sample <= view_start or start_sample >= view_end:
                    continue

                # Draw event boundary (vertical line)
                ax.axvline(start_sample, color="gray", alpha=0.2, linewidth=0.5)

                # Draw event mean level (observed)
                ax.hlines(
                    event_mean,
                    start_sample,
                    end_sample,
                    colors="red",
                    linewidth=2.5,
                    alpha=0.9,
                )

                # Draw raw model mean (unscaled)
                if model_mean > 0:
                    ax.hlines(
                        model_mean,
                        start_sample,
                        end_sample,
                        colors="orange",
                        linewidth=1.5,
                        linestyle=":",
                        alpha=0.7,
                    )

                # Draw scaled model expected level
                if scaled_model_mean > 0:
                    ax.hlines(
                        scaled_model_mean,
                        start_sample,
                        end_sample,
                        colors="limegreen",
                        linewidth=2,
                        linestyle="--",
                        alpha=0.8,
                    )

                # Annotate with k-mer and event info (sparse to avoid clutter)
                n_display = len(alignment_sorted)
                if i % max(1, n_display // 25) == 0 and kmer:
                    # Color by HMM state
                    state_colors = {"M": "black", "K": "darkorange", "B": "red"}
                    color = state_colors.get(state, "gray")

                    # Position label above the event mean
                    label_y = max(event_mean, scaled_model_mean, model_mean) + 3
                    ax.text(
                        (start_sample + end_sample) / 2,
                        label_y,
                        f"{kmer}\nidx:{event_idx} pos:{ref_pos}\nevent:{event_mean:.1f} model:{model_mean:.1f} scaled:{scaled_model_mean:.1f}",
                        fontsize=5,
                        ha="center",
                        va="bottom",
                        color=color,
                        rotation=0,
                        bbox=dict(
                            boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"
                        ),
                    )

            # Set axis limits to focus on aligned region
            ax.set_xlim(view_start, view_end)

            # Labels and title
            ax.set_xlabel("Sample index (raw signal)", fontsize=11)
            ax.set_ylabel("Signal level (pA)", fontsize=11)
            ax.set_title(
                f"Eventalign: {read_id[:30]}{'...' if len(read_id) > 30 else ''} → {transcript_id}\n"
                f"Total Events: {n_events}, Aligned: {n_aligned}, "
                f"Scale: {scaling.get('scale', 0):.4f}, Shift: {scaling.get('shift', 0):.2f}, "
                f"Var: {scaling.get('var', 0):.4f}",
                fontsize=12,
            )

            # Add legend
            from matplotlib.lines import Line2D

            legend_elements = [
                Line2D([0], [0], color="steelblue", alpha=0.5, linewidth=1, label="Raw signal"),
                Line2D([0], [0], color="red", linewidth=2.5, label="Event mean (observed)"),
                Line2D(
                    [0],
                    [0],
                    color="orange",
                    linewidth=1.5,
                    linestyle=":",
                    label="Model mean (raw)",
                ),
                Line2D(
                    [0],
                    [0],
                    color="limegreen",
                    linewidth=2,
                    linestyle="--",
                    label="Model mean (scaled)",
                ),
            ]
            ax.legend(handles=legend_elements, loc="upper right", fontsize=9)

            # Add grid for readability
            ax.grid(True, alpha=0.3, linestyle=":")

            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

            logger.info(f"  Saved eventalign visualization to: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to generate eventalign visualization: {e}")
            import traceback

            traceback.print_exc()
            return None

    def assign_reads_to_transcripts(
        self,
        read_alignments: List[Dict[str, Any]],
        dtw_distances: Optional[Dict[str, Any]] = None,
        sigma: float = 1.0,
        beta: float = 0.5,
    ) -> Dict[str, Any]:
        """
        Assign reads to transcripts using EM algorithm with coherence.

        Uses the eventalign quality scores as read-to-transcript distances,
        and optionally DTW distances for read-to-read coherence.

        Args:
            read_alignments: List of read alignment results from analyze_region()
            dtw_distances: Optional DTW distance matrix result
            sigma: Temperature parameter for softmax
            beta: Weight for coherence term (0 = no coherence)

        Returns:
            Dictionary with assignment results
        """
        if not ASSIGNMENT_AVAILABLE:
            logger.warning("EM assignment not available")
            return {"error": "EM assignment module not available"}

        if not read_alignments:
            return {"error": "No read alignments provided"}

        # Build list of reads and transcripts
        reads = []
        transcripts = set()
        for ra in read_alignments:
            reads.append(ra.get("read_id", "unknown"))
            for aln in ra.get("alignments", []):
                transcripts.add(aln.get("transcript_id", "unknown"))

        reads = list(reads)
        transcripts = sorted(list(transcripts))
        n_reads = len(reads)
        n_tx = len(transcripts)

        if n_reads == 0 or n_tx == 0:
            return {"error": "No reads or transcripts found"}

        # Build read-to-transcript distance matrix
        # Use 1 - overall_score as distance (higher quality = lower distance)
        tx_to_idx = {tx: i for i, tx in enumerate(transcripts)}
        dist_read_to_tx = np.ones((n_reads, n_tx)) * 10.0  # Default high distance

        for i, ra in enumerate(read_alignments):
            for aln in ra.get("alignments", []):
                tx_id = aln.get("transcript_id")
                if tx_id in tx_to_idx:
                    j = tx_to_idx[tx_id]
                    qa = aln.get("quality_assessment", {})
                    score = qa.get("overall_score", 0.0)
                    # Convert score to distance (higher score = lower distance)
                    dist_read_to_tx[i, j] = 1.0 - score

        # Build read-to-read distance matrix from DTW or use identity
        if dtw_distances and "matrix" in dtw_distances:
            dtw_matrix = np.array(dtw_distances["matrix"])
            dtw_read_ids = dtw_distances.get("read_ids", [])

            # Map DTW indices to our read indices
            dist_read_to_read = np.zeros((n_reads, n_reads))
            read_to_idx = {r: i for i, r in enumerate(reads)}
            dtw_read_to_idx = {r: i for i, r in enumerate(dtw_read_ids)}

            for i, ri in enumerate(reads):
                for j, rj in enumerate(reads):
                    if ri in dtw_read_to_idx and rj in dtw_read_to_idx:
                        di = dtw_read_to_idx[ri]
                        dj = dtw_read_to_idx[rj]
                        dist_read_to_read[i, j] = dtw_matrix[di, dj]
        else:
            # No DTW, use zeros (no coherence penalty)
            dist_read_to_read = np.zeros((n_reads, n_reads))

        # Normalize distances
        if dist_read_to_tx.max() > 0:
            dist_read_to_tx = dist_read_to_tx / dist_read_to_tx.max()
        if dist_read_to_read.max() > 0:
            dist_read_to_read = dist_read_to_read / dist_read_to_read.max()

        # Run EM assignment
        logger.info(f"Running EM assignment: {n_reads} reads -> {n_tx} transcripts")
        R, hard_assignments, log_likelihoods = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=sigma,
            beta=beta,
            verbose=True,
        )

        # Build result
        assignments = []
        for i, read_id in enumerate(reads):
            tx_idx = hard_assignments[i]
            tx_id = transcripts[tx_idx]
            confidence = R[i, tx_idx]
            assignments.append(
                {
                    "read_id": read_id,
                    "transcript_id": tx_id,
                    "confidence": float(confidence),
                    "probabilities": {tx: float(R[i, j]) for j, tx in enumerate(transcripts)},
                }
            )

        # Group reads by transcript
        reads_per_transcript = defaultdict(list)
        for a in assignments:
            reads_per_transcript[a["transcript_id"]].append(
                {
                    "read_id": a["read_id"],
                    "confidence": a["confidence"],
                }
            )

        return {
            "status": "success",
            "n_reads": n_reads,
            "n_transcripts": n_tx,
            "transcripts": transcripts,
            "reads": reads,
            "assignments": assignments,
            "reads_per_transcript": dict(reads_per_transcript),
            "responsibility_matrix": R.tolist(),
            "log_likelihoods": log_likelihoods,
            "parameters": {"sigma": sigma, "beta": beta},
        }

    def plot_transcript_read_assignments(
        self,
        region: Any,
        assignments: Dict[str, Any],
        read_alignments: List[Dict[str, Any]],
        output_path: Path,
        title: str = "Transcript-Read Assignments",
        figsize_per_panel: Tuple[int, int] = (16, 3),
        max_transcripts: int = 10,
    ) -> Optional[Path]:
        """
        Visualize transcript structures and assigned reads on genomic coordinates.

        Creates a multi-panel figure where each panel shows:
        - Transcript exon structure
        - Reads assigned to that transcript (from BAM mapping)
        - All panels share the same genomic x-axis for comparison

        Args:
            region: GenomicInterval with chrom, start, end
            assignments: Result from assign_reads_to_transcripts()
            read_alignments: List of read alignment results
            output_path: Path to save the figure
            title: Figure title
            figsize_per_panel: Size per transcript panel
            max_transcripts: Maximum number of transcripts to show

        Returns:
            Path to saved image, or None if visualization failed
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for visualization")
            return None

        if assignments.get("status") != "success":
            logger.warning(f"Assignment failed: {assignments.get('error', 'unknown')}")
            return None

        try:
            reads_per_tx = assignments.get("reads_per_transcript", {})
            transcripts = assignments.get("transcripts", [])

            # Limit transcripts to show
            if len(transcripts) > max_transcripts:
                # Sort by number of assigned reads, show top ones
                tx_counts = [(tx, len(reads_per_tx.get(tx, []))) for tx in transcripts]
                tx_counts.sort(key=lambda x: -x[1])
                transcripts = [tx for tx, _ in tx_counts[:max_transcripts]]

            n_panels = len(transcripts)
            if n_panels == 0:
                logger.warning("No transcripts to visualize")
                return None

            # Get genomic region bounds
            region_chrom = region.chrom
            region_start = region.start
            region_end = region.end
            region_strand = region.strand if hasattr(region, "strand") else None

            # Get transcript structures from GTF
            transcript_structures = {}
            with GTFReader(str(self.gtf_path)) as reader:
                reader.parse()
                for tx in reader.iterate_transcripts():
                    if tx.transcript_id in transcripts:
                        tx.sort_features()
                        transcript_structures[tx.transcript_id] = tx

            # Get read mapping info from BAM (filter by region strand)
            read_mapping_info = {}
            with BamReader(str(self.bam_path)) as bam_reader:
                for read in bam_reader.fetch(region=f"{region_chrom}:{region_start}-{region_end}"):
                    read_id = read.query_name
                    if read_id in assignments.get("reads", []):
                        # Filter by strand if region has strand information
                        if region_strand is not None:
                            read_strand = (
                                "-" if (hasattr(read, "is_reverse") and read.is_reverse) else "+"
                            )
                            if read_strand != region_strand:
                                continue
                        # Get read alignment blocks
                        blocks = []
                        if hasattr(read, "get_blocks"):
                            blocks = read.get_blocks()
                        elif hasattr(read, "reference_start") and hasattr(read, "reference_end"):
                            blocks = [(read.reference_start, read.reference_end)]

                        read_mapping_info[read_id] = {
                            "blocks": blocks,
                            "start": (
                                read.reference_start if hasattr(read, "reference_start") else 0
                            ),
                            "end": read.reference_end if hasattr(read, "reference_end") else 0,
                            "is_reverse": read.is_reverse if hasattr(read, "is_reverse") else False,
                        }

            # Create figure with shared x-axis
            fig, axes = plt.subplots(
                n_panels,
                1,
                figsize=(figsize_per_panel[0], figsize_per_panel[1] * n_panels),
                sharex=True,
            )

            if n_panels == 1:
                axes = [axes]

            # Color palette for reads
            colors = plt.cm.tab20.colors

            for panel_idx, tx_id in enumerate(transcripts):
                ax = axes[panel_idx]

                # Get transcript structure
                tx_struct = transcript_structures.get(tx_id)

                # Draw transcript exons
                exon_y = 0.8
                exon_height = 0.15

                if tx_struct:
                    # Draw introns as thin line
                    ax.hlines(
                        exon_y,
                        tx_struct.start,
                        tx_struct.end,
                        colors="darkblue",
                        linewidth=1,
                        alpha=0.5,
                    )

                    # Draw exons as thick boxes
                    for exon_start, exon_end in tx_struct.exons:
                        rect = plt.Rectangle(
                            (exon_start, exon_y - exon_height / 2),
                            exon_end - exon_start,
                            exon_height,
                            facecolor="darkblue",
                            edgecolor="black",
                            linewidth=0.5,
                        )
                        ax.add_patch(rect)

                    # Draw CDS if available (slightly different color)
                    for cds_start, cds_end in tx_struct.cds:
                        rect = plt.Rectangle(
                            (cds_start, exon_y - exon_height / 2),
                            cds_end - cds_start,
                            exon_height,
                            facecolor="navy",
                            edgecolor="black",
                            linewidth=0.5,
                        )
                        ax.add_patch(rect)

                # Draw assigned reads
                assigned_reads = reads_per_tx.get(tx_id, [])
                n_reads_assigned = len(assigned_reads)

                # Stack reads vertically
                read_y_start = 0.55
                read_height = 0.03
                read_spacing = 0.04
                max_stacked = 12  # Maximum reads to stack before condensing

                if n_reads_assigned > max_stacked:
                    # Condense display
                    read_height = 0.02
                    read_spacing = 0.025

                for read_idx, read_info in enumerate(assigned_reads[: min(n_reads_assigned, 30)]):
                    read_id = read_info["read_id"]
                    confidence = read_info["confidence"]

                    mapping = read_mapping_info.get(read_id)
                    if not mapping:
                        continue

                    read_y = read_y_start - (read_idx % max_stacked) * read_spacing
                    color = colors[read_idx % len(colors)]
                    alpha = 0.4 + 0.5 * confidence  # Higher confidence = more opaque

                    # Draw read blocks
                    blocks = mapping.get("blocks", [(mapping["start"], mapping["end"])])
                    for block_start, block_end in blocks:
                        rect = plt.Rectangle(
                            (block_start, read_y - read_height / 2),
                            block_end - block_start,
                            read_height,
                            facecolor=color,
                            edgecolor="gray",
                            linewidth=0.3,
                            alpha=alpha,
                        )
                        ax.add_patch(rect)

                    # Connect blocks with lines (for spliced reads)
                    if len(blocks) > 1:
                        for i in range(len(blocks) - 1):
                            ax.plot(
                                [blocks[i][1], blocks[i + 1][0]],
                                [read_y, read_y],
                                color=color,
                                linewidth=0.5,
                                alpha=alpha * 0.7,
                            )

                # Labels and styling
                ax.set_xlim(region_start, region_end)
                ax.set_ylim(0, 1)

                # Y-axis label with transcript info
                tx_label = tx_id[:25] + "..." if len(tx_id) > 25 else tx_id
                ax.set_ylabel(f"{tx_label}\n({n_reads_assigned} reads)", fontsize=9)
                ax.set_yticks([])

                # Add grid
                ax.axhline(0.65, color="gray", linestyle="--", alpha=0.3, linewidth=0.5)
                ax.grid(True, axis="x", alpha=0.3, linestyle=":")

                # Panel title
                if panel_idx == 0:
                    ax.set_title(
                        f"{title}\nRegion: {region_chrom}:{region_start:,}-{region_end:,}",
                        fontsize=12,
                        fontweight="bold",
                    )

            # X-axis label on bottom panel
            axes[-1].set_xlabel("Genomic Position", fontsize=11)

            # Format x-axis with comma separators
            from matplotlib.ticker import FuncFormatter

            axes[-1].xaxis.set_major_formatter(FuncFormatter(lambda x, p: f"{int(x):,}"))

            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

            logger.info(f"  Saved transcript-read assignment plot to: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to generate transcript-read assignment plot: {e}")
            import traceback

            traceback.print_exc()
            return None

    def plot_assignment_summary(
        self,
        assignments: Dict[str, Any],
        output_path: Path,
        title: str = "Read Assignment Summary",
        figsize: Tuple[int, int] = (14, 10),
    ) -> Optional[Path]:
        """
        Visualize assignment summary statistics.

        Creates a figure with:
        - Bar chart of reads per transcript
        - Confidence distribution histogram
        - Assignment probability heatmap

        Args:
            assignments: Result from assign_reads_to_transcripts()
            output_path: Path to save the figure
            title: Figure title
            figsize: Figure size

        Returns:
            Path to saved image, or None if visualization failed
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for visualization")
            return None

        if assignments.get("status") != "success":
            logger.warning(f"Assignment failed: {assignments.get('error', 'unknown')}")
            return None

        try:
            reads_per_tx = assignments.get("reads_per_transcript", {})
            transcripts = assignments.get("transcripts", [])
            reads = assignments.get("reads", [])
            R = np.array(assignments.get("responsibility_matrix", []))
            all_assignments = assignments.get("assignments", [])

            fig = plt.figure(figsize=figsize)
            gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

            # Panel 1: Reads per transcript (bar chart)
            ax1 = fig.add_subplot(gs[0, 0])
            tx_counts = [(tx, len(reads_per_tx.get(tx, []))) for tx in transcripts]
            tx_counts.sort(key=lambda x: -x[1])

            tx_labels = [t[:20] + ".." if len(t) > 20 else t for t, _ in tx_counts[:15]]
            counts = [c for _, c in tx_counts[:15]]

            bars = ax1.barh(range(len(tx_labels)), counts, color="steelblue", alpha=0.8)
            ax1.set_yticks(range(len(tx_labels)))
            ax1.set_yticklabels(tx_labels, fontsize=8)
            ax1.set_xlabel("Number of Assigned Reads")
            ax1.set_title("Reads per Transcript (Top 15)")
            ax1.invert_yaxis()

            # Add count labels
            for i, (bar, count) in enumerate(zip(bars, counts)):
                ax1.text(
                    bar.get_width() + 0.1,
                    bar.get_y() + bar.get_height() / 2,
                    str(count),
                    va="center",
                    fontsize=8,
                )

            # Panel 2: Confidence distribution
            ax2 = fig.add_subplot(gs[0, 1])
            confidences = [a["confidence"] for a in all_assignments]

            ax2.hist(confidences, bins=20, color="forestgreen", alpha=0.7, edgecolor="black")
            ax2.axvline(
                np.mean(confidences),
                color="red",
                linestyle="--",
                label=f"Mean: {np.mean(confidences):.3f}",
            )
            ax2.axvline(
                np.median(confidences),
                color="orange",
                linestyle="--",
                label=f"Median: {np.median(confidences):.3f}",
            )
            ax2.set_xlabel("Assignment Confidence")
            ax2.set_ylabel("Count")
            ax2.set_title("Confidence Distribution")
            ax2.legend(fontsize=9)

            # Panel 3: Responsibility matrix heatmap
            ax3 = fig.add_subplot(gs[1, :])

            # Sort reads by their primary assignment for better visualization
            read_order = sorted(
                range(len(reads)), key=lambda i: all_assignments[i]["transcript_id"]
            )
            R_sorted = R[read_order, :]

            # Limit display size
            max_reads_display = 50
            max_tx_display = 20

            if R_sorted.shape[0] > max_reads_display:
                R_display = R_sorted[:max_reads_display, :]
                read_labels = [
                    (
                        reads[read_order[i]][:12] + ".."
                        if len(reads[read_order[i]]) > 12
                        else reads[read_order[i]]
                    )
                    for i in range(max_reads_display)
                ]
            else:
                R_display = R_sorted
                read_labels = [
                    (
                        reads[read_order[i]][:12] + ".."
                        if len(reads[read_order[i]]) > 12
                        else reads[read_order[i]]
                    )
                    for i in range(len(reads))
                ]

            if R_display.shape[1] > max_tx_display:
                R_display = R_display[:, :max_tx_display]
                tx_labels = [
                    transcripts[i][:15] + ".." if len(transcripts[i]) > 15 else transcripts[i]
                    for i in range(max_tx_display)
                ]
            else:
                tx_labels = [
                    transcripts[i][:15] + ".." if len(transcripts[i]) > 15 else transcripts[i]
                    for i in range(len(transcripts))
                ]

            im = ax3.imshow(R_display, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
            cbar = plt.colorbar(im, ax=ax3, shrink=0.8)
            cbar.set_label("Assignment Probability")

            if len(read_labels) <= 50:
                ax3.set_yticks(range(len(read_labels)))
                ax3.set_yticklabels(read_labels, fontsize=6)
            else:
                ax3.set_ylabel(f"Reads (n={len(reads)})")

            if len(tx_labels) <= 20:
                ax3.set_xticks(range(len(tx_labels)))
                ax3.set_xticklabels(tx_labels, rotation=45, ha="right", fontsize=7)
            else:
                ax3.set_xlabel(f"Transcripts (n={len(transcripts)})")

            ax3.set_title("Assignment Probability Matrix (sorted by assignment)")
            ax3.set_xlabel("Transcript")
            ax3.set_ylabel("Read")

            fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

            plt.tight_layout()
            plt.savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close(fig)

            logger.info(f"  Saved assignment summary plot to: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to generate assignment summary: {e}")
            import traceback

            traceback.print_exc()
            return None

    def analyze_region(
        self, region, pod5_reader: Pod5Reader, max_reads: int = 10, max_transcripts: int = 50
    ) -> Dict[str, Any]:
        """
        Perform full analysis on a single region.

        Args:
            region: GenomicInterval object
            pod5_reader: Open Pod5Reader instance
            max_reads: Maximum number of reads to process
            max_transcripts: Maximum number of candidate transcripts

        Returns:
            Dictionary with analysis results
        """
        region_id = region.region_string
        logger.info(f"Analyzing region: {region_id}")

        result = {
            "region": {
                "chrom": region.chrom,
                "start": region.start,
                "end": region.end,
                "strand": region.strand,
                "read_count": region.read_count,
            },
            "candidate_transcripts": [],
            "read_alignments": [],
            "dtw_distances": None,
            "summary": {},
        }

        # Step 1: Get candidate transcripts for this region (filter by strand)
        region_strand = region.strand if hasattr(region, "strand") else None
        candidates = self.get_candidate_transcripts_for_region(
            region.chrom, region.start, region.end, strand=region_strand
        )

        if len(candidates) > max_transcripts:
            logger.info(f"  Limiting to {max_transcripts} transcripts (of {len(candidates)})")
            candidates = candidates[:max_transcripts]

        result["candidate_transcripts"] = [
            {"transcript_id": tx_id, "length": len(seq)} for tx_id, seq in candidates
        ]
        logger.info(
            f"  Found {len(candidates)} candidate transcripts on strand {region_strand if region_strand else 'both'}"
        )

        if not candidates:
            logger.warning(f"  No candidate transcripts found for region {region_id}")
            return result

        # Step 2: Get reads for this region (filter by strand)
        region_strand = region.strand if hasattr(region, "strand") else None
        reads = self.get_reads_for_region(
            region.chrom, region.start, region.end, strand=region_strand
        )

        if len(reads) > max_reads:
            logger.info(f"  Limiting to {max_reads} reads (of {len(reads)})")
            reads = reads[:max_reads]

        logger.info(
            f"  Processing {len(reads)} reads on strand {region_strand if region_strand else 'both'}"
        )

        # Step 3: For each read, get signal and align to all transcripts
        read_signals = []

        # Track first read for visualization
        first_read_info = None

        for read in reads:
            read_id = read.query_name

            # Get signal from POD5
            signal = self.get_signal_for_read(read_id, pod5_reader)
            if signal is None:
                logger.debug(f"  No signal found for read {read_id}")
                continue

            read_signals.append((read_id, signal))

            # Align to each candidate transcript
            read_result = {"read_id": read_id, "signal_length": len(signal), "alignments": []}

            for tx_id, tx_seq in candidates:
                if len(tx_seq) < self.kmer_size:
                    continue

                alignment = self.eventalign_read_to_transcript(signal, tx_seq)

                if alignment:
                    # Save first read's first alignment for visualization
                    if first_read_info is None:
                        first_read_info = {
                            "read_id": read_id,
                            "signal": signal,
                            "transcript_id": tx_id,
                            "alignment": alignment,
                        }

                    # Include full alignment details, not just summary
                    alignment_record = {
                        "transcript_id": tx_id,
                        "transcript_length": len(tx_seq),
                        "n_events": alignment.get("n_events", 0),
                        "n_aligned": alignment.get(
                            "n_aligned", alignment.get("n_aligned_pairs", 0)
                        ),
                        "events_per_base": alignment.get("events_per_base", 0),
                        "scaling": alignment.get("scaling", {}),
                        # Quality assessment
                        "quality_assessment": assess_eventalign_quality(
                            alignment, tx_seq, len(signal)
                        ),
                        # Full alignment details per position
                        "alignment_details": [
                            {
                                "ref_position": aln.get("ref_position"),
                                "ref_kmer": aln.get("ref_kmer"),
                                "event_idx": aln.get("event_idx"),
                                "signal_start": aln.get("signal_start"),
                                "signal_length": aln.get("signal_length"),
                                "hmm_state": aln.get("hmm_state"),
                                "event_mean": aln.get("event_mean"),
                                "event_stdv": aln.get("event_stdv"),
                                "event_duration": aln.get("event_duration"),
                                "model_mean": aln.get("model_mean"),
                                "model_stdv": aln.get("model_stdv"),
                                "scaled_model_mean": aln.get("scaled_model_mean"),
                                "scaled_model_stdv": aln.get("scaled_model_stdv"),
                            }
                            for aln in alignment.get("alignment", [])
                        ],
                    }
                    read_result["alignments"].append(alignment_record)

            result["read_alignments"].append(read_result)

        logger.info(f"  Aligned {len(result['read_alignments'])} reads to transcripts")

        # Step 4: Visualize the first read's eventalign result
        if first_read_info is not None:
            logger.info(
                f"  Visualizing eventalign for first read: {first_read_info['read_id'][:20]}..."
            )
            region_safe_id = region_id.replace(":", "_").replace("-", "_")
            eventalign_plot_path = self.output_dir / f"eventalign_{region_safe_id}_first_read.png"

            eventalign_plot = self.plot_eventalign_signal(
                signal=first_read_info["signal"],
                eventalign_result=first_read_info["alignment"],
                output_path=eventalign_plot_path,
                read_id=first_read_info["read_id"],
                transcript_id=first_read_info["transcript_id"],
            )

            if eventalign_plot:
                result["first_read_eventalign_plot"] = str(eventalign_plot)

        # Step 5: Compute DTW distances among read signals
        if len(read_signals) >= 2:
            logger.info(f"  Computing DTW distances for {len(read_signals)} reads...")
            signals_only = [sig for _, sig in read_signals]
            read_ids_list = [rid for rid, _ in read_signals]
            distance_matrix = self.compute_dtw_distances(signals_only)

            if distance_matrix is not None:
                result["dtw_distances"] = {
                    "read_ids": read_ids_list,
                    "matrix": distance_matrix.tolist(),
                }
                logger.info(f"  DTW distance matrix computed: {distance_matrix.shape}")

                # Step 6: Perform hierarchical clustering
                logger.info(f"  Performing hierarchical clustering...")
                cluster_result = self.cluster_reads_by_dtw(
                    distance_matrix, read_ids_list, method="average"
                )

                if "error" not in cluster_result:
                    result["clustering"] = cluster_result
                    logger.info(f"  Identified {cluster_result['n_clusters']} clusters")

                    # Log cluster membership
                    for cluster_id, members in cluster_result["clusters"].items():
                        logger.info(f"    Cluster {cluster_id}: {len(members)} reads")
                else:
                    logger.warning(f"  Clustering failed: {cluster_result['error']}")
                    result["clustering"] = None

                # Step 7: Generate heatmap visualization
                region_safe_id = region_id.replace(":", "_").replace("-", "_")
                heatmap_path = self.output_dir / f"dtw_heatmap_{region_safe_id}.png"

                heatmap_file = self.plot_dtw_heatmap(
                    distance_matrix,
                    read_ids_list,
                    heatmap_path,
                    title=f"DTW Distance Matrix - {region_id}",
                    cluster_result=result.get("clustering"),
                )

                if heatmap_file:
                    result["heatmap_path"] = str(heatmap_file)

        # Step 8: Generate assessment heatmaps
        if result["read_alignments"]:
            region_safe_id = region_id.replace(":", "_").replace("-", "_")

            # Single metric heatmap (overall score)
            assessment_heatmap_path = self.output_dir / f"assessment_heatmap_{region_safe_id}.png"
            assessment_heatmap = self.plot_assessment_heatmap(
                result["read_alignments"],
                assessment_heatmap_path,
                metric="overall_score",
                title=f"Eventalign Quality - {region_id}",
            )
            if assessment_heatmap:
                result["assessment_heatmap_path"] = str(assessment_heatmap)

            # Multi-metric heatmap
            multi_metric_path = self.output_dir / f"assessment_multi_{region_safe_id}.png"
            multi_metric_heatmap = self.plot_assessment_multi_metric(
                result["read_alignments"],
                multi_metric_path,
                title=f"Eventalign Quality Assessment - {region_id}",
            )
            if multi_metric_heatmap:
                result["assessment_multi_metric_path"] = str(multi_metric_heatmap)

        # Step 9: Assign reads to transcripts using EM algorithm
        if result["read_alignments"] and ASSIGNMENT_AVAILABLE:
            region_safe_id = region_id.replace(":", "_").replace("-", "_")
            logger.info(f"  Running EM read-to-transcript assignment...")

            assignments = self.assign_reads_to_transcripts(
                result["read_alignments"],
                dtw_distances=result.get("dtw_distances"),
                sigma=1.0,
                beta=0.5 if result.get("dtw_distances") else 0.0,
            )

            if assignments.get("status") == "success":
                result["assignments"] = assignments
                logger.info(
                    f"  Assigned {assignments['n_reads']} reads to {assignments['n_transcripts']} transcripts"
                )

                # Step 10: Visualize transcript-read assignments on genome
                tx_read_plot_path = self.output_dir / f"transcript_reads_{region_safe_id}.png"
                tx_read_plot = self.plot_transcript_read_assignments(
                    region,
                    assignments,
                    result["read_alignments"],
                    tx_read_plot_path,
                    title=f"Read Assignments - {region_id}",
                )
                if tx_read_plot:
                    result["transcript_reads_plot_path"] = str(tx_read_plot)

                # Step 11: Visualize assignment summary
                assignment_summary_path = (
                    self.output_dir / f"assignment_summary_{region_safe_id}.png"
                )
                assignment_summary = self.plot_assignment_summary(
                    assignments,
                    assignment_summary_path,
                    title=f"Assignment Summary - {region_id}",
                )
                if assignment_summary:
                    result["assignment_summary_plot_path"] = str(assignment_summary)
            else:
                logger.warning(f"  Assignment failed: {assignments.get('error', 'unknown')}")

        # Summary statistics
        result["summary"] = {
            "num_transcripts": len(candidates),
            "num_reads_processed": len(result["read_alignments"]),
            "num_reads_with_signal": len(read_signals),
            "avg_alignments_per_read": (
                sum(len(r["alignments"]) for r in result["read_alignments"])
                / max(len(result["read_alignments"]), 1)
            ),
            "dtw_computed": result["dtw_distances"] is not None,
            "clustering_computed": result.get("clustering") is not None,
            "n_clusters": (
                result.get("clustering", {}).get("n_clusters", 0) if result.get("clustering") else 0
            ),
            "assignment_computed": result.get("assignments") is not None,
        }

        return result

    def run_analysis(
        self,
        max_regions: Optional[int] = None,
        max_reads_per_region: int = 10,
        max_transcripts_per_region: int = 50,
    ) -> Dict[str, Any]:
        """
        Run the complete analysis workflow.

        Args:
            max_regions: Maximum number of regions to analyze (None = all)
            max_reads_per_region: Maximum reads per region
            max_transcripts_per_region: Maximum transcript candidates per region

        Returns:
            Complete analysis results dictionary
        """
        logger.info("=" * 70)
        logger.info("Starting Region-based Transcript Analysis Workflow")
        logger.info("=" * 70)

        # Step 1: Load references
        self.load_references()

        # Step 2: Load annotation
        self.load_annotation()

        # Step 3: Generate isolated regions
        regions = self.generate_regions()

        if max_regions and len(regions) > max_regions:
            logger.info(f"Limiting analysis to {max_regions} regions")
            regions = regions[:max_regions]

        # Step 4: Analyze each region
        logger.info(f"Analyzing {len(regions)} regions...")

        with Pod5Reader(str(self.pod5_path)) as pod5_reader:
            for i, region in enumerate(regions):
                region_id = region.region_string
                logger.info(f"\n[{i+1}/{len(regions)}] Processing {region_id}")

                region_result = self.analyze_region(
                    region,
                    pod5_reader,
                    max_reads=max_reads_per_region,
                    max_transcripts=max_transcripts_per_region,
                )

                self.region_results[region_id] = region_result
                break

        # Compile final results
        final_results = {
            "workflow": "region_transcript_analysis",
            "parameters": {
                "bam": str(self.bam_path),
                "genome": str(self.genome_path),
                "transcriptome": str(self.transcriptome_path),
                "gtf": str(self.gtf_path),
                "pod5": str(self.pod5_path),
                "mode": "RNA-only",
                "kmer_size": self.kmer_size,
            },
            "summary": {
                "num_regions_analyzed": len(self.region_results),
                "total_transcripts": sum(
                    r["summary"]["num_transcripts"] for r in self.region_results.values()
                ),
                "total_reads_processed": sum(
                    r["summary"]["num_reads_processed"] for r in self.region_results.values()
                ),
                "regions_with_dtw": sum(
                    1 for r in self.region_results.values() if r["summary"]["dtw_computed"]
                ),
                "regions_with_clustering": sum(
                    1
                    for r in self.region_results.values()
                    if r["summary"].get("clustering_computed", False)
                ),
                "total_clusters": sum(
                    r["summary"].get("n_clusters", 0) for r in self.region_results.values()
                ),
                "regions_with_assignments": sum(
                    1
                    for r in self.region_results.values()
                    if r["summary"].get("assignment_computed", False)
                ),
            },
            "regions": self.region_results,
        }

        # Save results
        output_file = self.output_dir / "analysis_results.json"
        with open(output_file, "w") as f:
            json.dump(final_results, f, indent=2)
        logger.info(f"\nResults saved to: {output_file}")

        return final_results


def assess_eventalign_quality(
    eventalign_result: Dict[str, Any],
    sequence: str,
    signal_length: int,
) -> Dict[str, Any]:
    """
    Comprehensive quality assessment of eventalign results.

    Evaluates the alignment quality based on multiple metrics including
    coverage, signal-model correlation, HMM state distribution, event
    continuity, and scaling parameters.

    Args:
        eventalign_result: Result dictionary from profile_hmm_eventalign()
        sequence: The reference sequence that was aligned to
        signal_length: Length of the raw signal array

    Returns:
        Dictionary containing quality metrics and overall assessment
    """
    alignment = eventalign_result.get("alignment", [])
    scaling = eventalign_result.get("scaling", {})
    n_events = eventalign_result.get("n_events", 0)
    n_aligned = eventalign_result.get("n_aligned", len(alignment))

    if not alignment or n_events == 0:
        return {
            "status": "failed",
            "reason": "No alignment data",
            "overall_score": 0.0,
        }

    # Filter valid alignments (with positive event_idx)
    valid_alignments = [a for a in alignment if a.get("event_idx", -1) >= 0]

    if not valid_alignments:
        return {
            "status": "failed",
            "reason": "No valid aligned events",
            "overall_score": 0.0,
        }

    # ========== 1. Coverage Metrics ==========
    # Event coverage: what fraction of detected events were aligned
    event_coverage = n_aligned / n_events if n_events > 0 else 0.0

    # Sequence coverage: what fraction of sequence positions have alignments
    aligned_positions = set(a.get("ref_position", -1) for a in valid_alignments)
    aligned_positions.discard(-1)
    seq_len = len(sequence)
    sequence_coverage = len(aligned_positions) / seq_len if seq_len > 0 else 0.0

    # Signal coverage: what fraction of raw signal is covered by aligned events
    signal_starts = [a.get("signal_start", 0) for a in valid_alignments]
    signal_lengths = [a.get("signal_length", 0) for a in valid_alignments]
    if signal_starts and signal_lengths:
        min_signal = min(signal_starts)
        max_signal = max(s + l for s, l in zip(signal_starts, signal_lengths))
        aligned_signal_span = max_signal - min_signal
        signal_coverage = aligned_signal_span / signal_length if signal_length > 0 else 0.0
    else:
        signal_coverage = 0.0
        min_signal = 0
        max_signal = 0

    # ========== 2. HMM State Distribution ==========
    hmm_states = [a.get("hmm_state", "?") for a in valid_alignments]
    state_counts = {}
    for s in hmm_states:
        state_counts[s] = state_counts.get(s, 0) + 1

    total_states = len(hmm_states)
    match_fraction = state_counts.get("M", 0) / total_states if total_states > 0 else 0.0
    bad_event_fraction = state_counts.get("B", 0) / total_states if total_states > 0 else 0.0
    kmer_skip_fraction = state_counts.get("K", 0) / total_states if total_states > 0 else 0.0

    # ========== 3. Signal-Model Correlation ==========
    event_means = []
    model_means = []
    for a in valid_alignments:
        em = a.get("event_mean")
        mm = a.get("scaled_model_mean", a.get("model_mean"))
        if em is not None and mm is not None and mm > 0:
            event_means.append(em)
            model_means.append(mm)

    if len(event_means) >= 2:
        event_means_arr = np.array(event_means)
        model_means_arr = np.array(model_means)

        # Pearson correlation
        correlation = np.corrcoef(event_means_arr, model_means_arr)[0, 1]
        if np.isnan(correlation):
            correlation = 0.0

        # Mean absolute error
        mae = np.mean(np.abs(event_means_arr - model_means_arr))

        # Root mean squared error
        rmse = np.sqrt(np.mean((event_means_arr - model_means_arr) ** 2))

        # Normalized RMSE (as fraction of signal range)
        signal_range = np.max(event_means_arr) - np.min(event_means_arr)
        nrmse = rmse / signal_range if signal_range > 0 else 1.0

        # R-squared
        ss_res = np.sum((event_means_arr - model_means_arr) ** 2)
        ss_tot = np.sum((event_means_arr - np.mean(event_means_arr)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    else:
        correlation = 0.0
        mae = float("inf")
        rmse = float("inf")
        nrmse = 1.0
        r_squared = 0.0

    # ========== 4. Event Continuity (gaps in event indices) ==========
    event_indices = sorted([a.get("event_idx", 0) for a in valid_alignments])
    if len(event_indices) >= 2:
        index_diffs = np.diff(event_indices)
        # Ideally consecutive events have diff of 1 or small values
        continuity_score = np.mean(index_diffs <= 3) if len(index_diffs) > 0 else 0.0
        max_gap = int(np.max(index_diffs)) if len(index_diffs) > 0 else 0
        mean_gap = float(np.mean(index_diffs)) if len(index_diffs) > 0 else 0.0
    else:
        continuity_score = 0.0
        max_gap = 0
        mean_gap = 0.0

    # ========== 5. Sequence Continuity (gaps in ref positions) ==========
    ref_positions = sorted([a.get("ref_position", 0) for a in valid_alignments])
    if len(ref_positions) >= 2:
        pos_diffs = np.diff(ref_positions)
        seq_continuity_score = np.mean(pos_diffs <= 2) if len(pos_diffs) > 0 else 0.0
        max_pos_gap = int(np.max(pos_diffs)) if len(pos_diffs) > 0 else 0
    else:
        seq_continuity_score = 0.0
        max_pos_gap = 0

    # ========== 6. Events per Base ==========
    events_per_base = n_aligned / seq_len if seq_len > 0 else 0.0
    # Typical range is 1-5 events per base for good alignments
    epb_score = (
        min(1.0, events_per_base / 2.0)
        if events_per_base <= 5
        else max(0.0, 1.0 - (events_per_base - 5) / 10)
    )

    # ========== 7. Event Duration Statistics ==========
    durations = [
        a.get("event_duration", 0) for a in valid_alignments if a.get("event_duration", 0) > 0
    ]
    if durations:
        mean_duration = float(np.mean(durations))
        std_duration = float(np.std(durations))
        min_duration = float(np.min(durations))
        max_duration = float(np.max(durations))
    else:
        mean_duration = 0.0
        std_duration = 0.0
        min_duration = 0.0
        max_duration = 0.0

    # ========== 8. Scaling Parameter Quality ==========
    scale = scaling.get("scale", 1.0)
    shift = scaling.get("shift", 0.0)
    var = scaling.get("var", 1.0)

    # Reasonable ranges for RNA (based on typical nanopore models)
    # scale: ~0.8-1.2, shift: ~-20 to +20, var: ~0.8-1.2
    scale_ok = 0.5 <= scale <= 2.0
    shift_ok = -50 <= shift <= 50
    var_ok = 0.5 <= var <= 2.0
    scaling_quality = sum([scale_ok, shift_ok, var_ok]) / 3.0

    # ========== 9. Soft-clipping Assessment ==========
    # Check how much of the signal was soft-clipped (before/after alignment)
    if signal_starts:
        pre_clip_samples = int(min(signal_starts))
        post_clip_samples = signal_length - int(
            max(s + l for s, l in zip(signal_starts, signal_lengths))
        )
        pre_clip_fraction = pre_clip_samples / signal_length if signal_length > 0 else 0.0
        post_clip_fraction = post_clip_samples / signal_length if signal_length > 0 else 0.0
        total_clip_fraction = pre_clip_fraction + post_clip_fraction
    else:
        pre_clip_samples = 0
        post_clip_samples = 0
        pre_clip_fraction = 0.0
        post_clip_fraction = 0.0
        total_clip_fraction = 0.0

    # ========== 10. Overall Quality Score ==========
    # Weighted combination of key metrics
    weights = {
        "correlation": 0.25,  # Signal-model correlation is very important
        "match_fraction": 0.15,  # High match state fraction is good
        "event_coverage": 0.10,  # Good event utilization
        "sequence_coverage": 0.15,  # Good sequence coverage
        "continuity": 0.10,  # Continuous alignment without big gaps
        "scaling": 0.10,  # Reasonable scaling parameters
        "epb": 0.05,  # Reasonable events per base
        "low_bad_events": 0.10,  # Low bad event fraction
    }

    # Normalize correlation to 0-1 (correlation can be negative)
    norm_correlation = (correlation + 1) / 2

    overall_score = (
        weights["correlation"] * norm_correlation
        + weights["match_fraction"] * match_fraction
        + weights["event_coverage"] * event_coverage
        + weights["sequence_coverage"] * sequence_coverage
        + weights["continuity"] * continuity_score
        + weights["scaling"] * scaling_quality
        + weights["epb"] * epb_score
        + weights["low_bad_events"] * (1 - bad_event_fraction)
    )

    # Quality grade
    if overall_score >= 0.8:
        grade = "Excellent"
    elif overall_score >= 0.6:
        grade = "Good"
    elif overall_score >= 0.4:
        grade = "Fair"
    elif overall_score >= 0.2:
        grade = "Poor"
    else:
        grade = "Failed"

    # Build assessment result
    assessment = {
        "status": "success",
        "overall_score": round(overall_score, 4),
        "grade": grade,
        # Coverage metrics
        "coverage": {
            "event_coverage": round(event_coverage, 4),
            "sequence_coverage": round(sequence_coverage, 4),
            "signal_coverage": round(signal_coverage, 4),
            "aligned_positions": len(aligned_positions),
            "sequence_length": seq_len,
            "n_events": n_events,
            "n_aligned": n_aligned,
            "events_per_base": round(events_per_base, 3),
        },
        # HMM state distribution
        "hmm_states": {
            "match_fraction": round(match_fraction, 4),
            "bad_event_fraction": round(bad_event_fraction, 4),
            "kmer_skip_fraction": round(kmer_skip_fraction, 4),
            "state_counts": state_counts,
        },
        # Signal-model fit
        "signal_model_fit": {
            "correlation": round(correlation, 4),
            "r_squared": round(r_squared, 4),
            "mae": round(mae, 3),
            "rmse": round(rmse, 3),
            "nrmse": round(nrmse, 4),
        },
        # Continuity metrics
        "continuity": {
            "event_continuity_score": round(continuity_score, 4),
            "sequence_continuity_score": round(seq_continuity_score, 4),
            "max_event_gap": max_gap,
            "mean_event_gap": round(mean_gap, 2),
            "max_position_gap": max_pos_gap,
        },
        # Event duration statistics
        "event_duration": {
            "mean": round(mean_duration, 3),
            "std": round(std_duration, 3),
            "min": round(min_duration, 3),
            "max": round(max_duration, 3),
        },
        # Scaling parameters
        "scaling": {
            "scale": round(scale, 4),
            "shift": round(shift, 3),
            "var": round(var, 4),
            "quality": round(scaling_quality, 2),
            "scale_in_range": scale_ok,
            "shift_in_range": shift_ok,
            "var_in_range": var_ok,
        },
        # Soft-clipping info
        "soft_clipping": {
            "pre_clip_samples": pre_clip_samples,
            "post_clip_samples": max(0, post_clip_samples),
            "pre_clip_fraction": round(pre_clip_fraction, 4),
            "post_clip_fraction": round(max(0, post_clip_fraction), 4),
            "total_clip_fraction": round(max(0, total_clip_fraction), 4),
            "aligned_signal_range": [int(min_signal), int(max_signal)],
        },
        # Diagnostic messages
        "diagnostics": [],
    }

    # Add diagnostic messages
    diagnostics = assessment["diagnostics"]

    if correlation < 0.5:
        diagnostics.append(
            f"Low signal-model correlation ({correlation:.3f}). Possible sequence mismatch or poor scaling."
        )
    if event_coverage < 0.3:
        diagnostics.append(
            f"Low event coverage ({event_coverage:.1%}). Many events were not aligned."
        )
    if sequence_coverage < 0.5:
        diagnostics.append(
            f"Low sequence coverage ({sequence_coverage:.1%}). Large portions of sequence not aligned."
        )
    if bad_event_fraction > 0.3:
        diagnostics.append(
            f"High bad event fraction ({bad_event_fraction:.1%}). Signal quality may be poor."
        )
    if max_gap > 50:
        diagnostics.append(
            f"Large event index gap ({max_gap}). Possible signal dropout or alignment issue."
        )
    if not scale_ok:
        diagnostics.append(f"Unusual scale parameter ({scale:.4f}). Expected 0.5-2.0.")
    if not shift_ok:
        diagnostics.append(f"Unusual shift parameter ({shift:.1f}). Expected -50 to +50.")
    if events_per_base < 0.5:
        diagnostics.append(
            f"Very low events per base ({events_per_base:.2f}). Possible sequence length mismatch."
        )
    if events_per_base > 10:
        diagnostics.append(
            f"Very high events per base ({events_per_base:.2f}). Possible over-segmentation."
        )
    if total_clip_fraction > 0.5:
        diagnostics.append(
            f"High soft-clipping ({total_clip_fraction:.1%}). Much of signal not aligned."
        )

    if not diagnostics:
        diagnostics.append("No issues detected. Alignment quality is good.")

    return assessment


def print_assessment_report(assessment: Dict[str, Any], verbose: bool = True) -> None:
    """
    Print a formatted assessment report to stdout.

    Args:
        assessment: Result from assess_eventalign_quality()
        verbose: If True, print detailed metrics; otherwise print summary
    """
    print("\n" + "=" * 70)
    print("EVENTALIGN QUALITY ASSESSMENT REPORT")
    print("=" * 70)

    if assessment.get("status") == "failed":
        print(f"\n❌ Assessment FAILED: {assessment.get('reason', 'Unknown error')}")
        return

    score = assessment["overall_score"]
    grade = assessment["grade"]

    # Grade emoji
    grade_emoji = {"Excellent": "🌟", "Good": "✅", "Fair": "⚠️", "Poor": "❌", "Failed": "💀"}
    emoji = grade_emoji.get(grade, "❓")

    print(f"\n{emoji} Overall Quality: {grade} (Score: {score:.2%})")
    print("-" * 70)

    # Coverage summary
    cov = assessment["coverage"]
    print(f"\n📊 COVERAGE:")
    print(
        f"   Events aligned:    {cov['n_aligned']:,} / {cov['n_events']:,} ({cov['event_coverage']:.1%})"
    )
    print(
        f"   Sequence covered:  {cov['aligned_positions']:,} / {cov['sequence_length']:,} positions ({cov['sequence_coverage']:.1%})"
    )
    print(f"   Signal covered:    {cov['signal_coverage']:.1%}")
    print(f"   Events per base:   {cov['events_per_base']:.2f}")

    # HMM states
    hmm = assessment["hmm_states"]
    print(f"\n🔄 HMM STATES:")
    print(f"   Match (M):         {hmm['match_fraction']:.1%}")
    print(f"   Bad Event (B):     {hmm['bad_event_fraction']:.1%}")
    print(f"   K-mer Skip (K):    {hmm['kmer_skip_fraction']:.1%}")

    # Signal-model fit
    fit = assessment["signal_model_fit"]
    print(f"\n📈 SIGNAL-MODEL FIT:")
    print(f"   Correlation (r):   {fit['correlation']:.4f}")
    print(f"   R-squared:         {fit['r_squared']:.4f}")
    print(f"   MAE:               {fit['mae']:.2f} pA")
    print(f"   RMSE:              {fit['rmse']:.2f} pA")

    if verbose:
        # Continuity
        cont = assessment["continuity"]
        print(f"\n🔗 CONTINUITY:")
        print(f"   Event continuity:  {cont['event_continuity_score']:.1%}")
        print(f"   Sequence continuity: {cont['sequence_continuity_score']:.1%}")
        print(f"   Max event gap:     {cont['max_event_gap']}")
        print(f"   Max position gap:  {cont['max_position_gap']}")

        # Event duration
        dur = assessment["event_duration"]
        print(f"\n⏱️ EVENT DURATION:")
        print(f"   Mean:              {dur['mean']:.1f} samples")
        print(f"   Std:               {dur['std']:.1f} samples")
        print(f"   Range:             [{dur['min']:.1f}, {dur['max']:.1f}]")

        # Scaling
        scl = assessment["scaling"]
        print(f"\n⚖️ SCALING PARAMETERS:")
        print(f"   Scale:             {scl['scale']:.4f} {'✓' if scl['scale_in_range'] else '⚠️'}")
        print(f"   Shift:             {scl['shift']:.2f} {'✓' if scl['shift_in_range'] else '⚠️'}")
        print(f"   Var:               {scl['var']:.4f} {'✓' if scl['var_in_range'] else '⚠️'}")

        # Soft-clipping
        clip = assessment["soft_clipping"]
        print(f"\n✂️ SOFT-CLIPPING:")
        print(
            f"   Pre-clip:          {clip['pre_clip_samples']:,} samples ({clip['pre_clip_fraction']:.1%})"
        )
        print(
            f"   Post-clip:         {clip['post_clip_samples']:,} samples ({clip['post_clip_fraction']:.1%})"
        )
        print(
            f"   Aligned range:     [{clip['aligned_signal_range'][0]:,}, {clip['aligned_signal_range'][1]:,}]"
        )

    # Diagnostics
    print(f"\n💡 DIAGNOSTICS:")
    for diag in assessment["diagnostics"]:
        print(f"   • {diag}")

    print("\n" + "=" * 70)


def run_demo():
    """Run a demonstration with synthetic data."""
    print("=" * 70)
    print("Region-based Transcript Analysis Workflow - DEMO MODE")
    print("=" * 70)
    print()
    print("This workflow requires the following input files:")
    print("  1. BAM file - reads mapped to the genome")
    print("  2. Genome FASTA - reference genome sequences")
    print("  3. Transcriptome FASTA - reference transcript sequences")
    print("  4. GTF file - gene/transcript annotations")
    print("  5. POD5 file - nanopore raw signal data")
    print()
    print("The workflow performs:")
    print("  1. Separate reads into isolated genomic regions")
    print("  2. For each region, identify candidate transcripts")
    print("  3. Align each read's signal to each candidate transcript")
    print("  4. Compute DTW distances among reads in the region")
    print()

    # Demonstrate with synthetic data if eventalign is available
    if EVENTALIGN_AVAILABLE:
        print("Demonstrating eventalign with synthetic data...")
        print("-" * 70)

        # Generate synthetic signal
        np.random.seed(42)
        signal = np.random.randn(1000).astype(np.float32) * 10 + 120
        sequence = "ACGUACGUACGUACGUACGUACGUACGU"

        # RNA-only mode: events are reversed internally to match 5'→3' sequence
        # No need to reverse the sequence - use it directly
        print(f"Sequence (5'→3'): {sequence}")
        print()

        result = profile_hmm_eventalign(raw_signal=signal, sequence=sequence, kmer_size=5)

        print(f"Signal length: {len(signal)} samples")
        print(f"Sequence length: {len(sequence)} bases")
        print(f"Events detected: {result['n_events']}")
        print(f"Alignment records: {result['n_aligned']}")
        print(
            f"Scaling: scale={result['scaling']['scale']:.3f}, "
            f"shift={result['scaling']['shift']:.1f}"
        )
        print()
        print("First 5 alignment records:")
        for aln in result["alignment"][:5]:
            print(
                f"  Position {aln['ref_position']}: {aln['ref_kmer']} "
                f"(state={aln['hmm_state']}, event_idx={aln.get('event_idx', -1)}, "
                f"signal_start={aln.get('signal_start', 0)}, "
                f"event_mean={aln['event_mean']:.1f}, scaled_model_mean={aln.get('scaled_model_mean', 0):.1f})"
            )

        # Run quality assessment
        print()
        print("-" * 70)
        print("Running eventalign quality assessment...")
        assessment = assess_eventalign_quality(result, sequence, len(signal))
        print_assessment_report(assessment, verbose=True)

    else:
        print("Note: eventalign extension not available.")
        print("Build with: pip install -e .")

    print()
    print("To run the full workflow:")
    print("  python region_transcript_analysis_workflow.py \\")
    print("      --bam reads.bam \\")
    print("      --genome genome.fa \\")
    print("      --transcriptome transcripts.fa \\")
    print("      --gtf annotation.gtf \\")
    print("      --pod5 signals.pod5 \\")
    print("      --output results/")


def main():
    parser = argparse.ArgumentParser(
        description="Region-based Transcript Analysis Workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("--bam", "-b", type=str, help="Path to BAM file (reads mapped to genome)")
    parser.add_argument("--genome", "-g", type=str, help="Path to genome FASTA reference")
    parser.add_argument(
        "--transcriptome", "-t", type=str, help="Path to transcriptome FASTA reference"
    )
    parser.add_argument("--gtf", type=str, help="Path to GTF annotation file")
    parser.add_argument("--pod5", "-p", type=str, help="Path to POD5 signal file")
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="./analysis_results",
        help="Output directory (default: ./analysis_results)",
    )
    parser.add_argument(
        "--max-regions", type=int, default=None, help="Maximum number of regions to analyze"
    )
    parser.add_argument(
        "--max-reads", type=int, default=100, help="Maximum reads per region (default: 100)"
    )
    parser.add_argument(
        "--max-transcripts",
        type=int,
        default=50,
        help="Maximum transcript candidates per region (default: 50)",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=5,
        choices=[5, 9],
        help="K-mer size for eventalign (default: 5)",
    )
    # RNA-only mode: no --dna option needed
    parser.add_argument("--demo", action="store_true", help="Run demonstration with synthetic data")

    args = parser.parse_args()

    # Run demo if requested or if no input files provided
    if args.demo or not any([args.bam, args.genome, args.transcriptome, args.gtf, args.pod5]):
        run_demo()
        return

    # Validate required arguments
    required = ["bam", "genome", "transcriptome", "gtf", "pod5"]
    missing = [arg for arg in required if not getattr(args, arg)]

    if missing:
        print(f"Error: Missing required arguments: {', '.join('--' + m for m in missing)}")
        print("Run with --demo to see example usage.")
        sys.exit(1)

    # Validate files exist
    for arg in required:
        path = getattr(args, arg)
        if not Path(path).exists():
            print(f"Error: File not found: {path}")
            sys.exit(1)

    # Run analysis
    analyzer = RegionTranscriptAnalyzer(
        bam_path=args.bam,
        genome_path=args.genome,
        transcriptome_path=args.transcriptome,
        gtf_path=args.gtf,
        pod5_path=args.pod5,
        output_dir=args.output,
        kmer_size=args.kmer_size,
    )

    results = analyzer.run_analysis(
        max_regions=args.max_regions,
        max_reads_per_region=args.max_reads,
        max_transcripts_per_region=args.max_transcripts,
    )

    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print(f"Regions analyzed: {results['summary']['num_regions_analyzed']}")
    print(f"Total transcripts: {results['summary']['total_transcripts']}")
    print(f"Total reads processed: {results['summary']['total_reads_processed']}")
    print(f"Regions with DTW: {results['summary']['regions_with_dtw']}")
    print(f"Regions with clustering: {results['summary']['regions_with_clustering']}")
    print(f"Total clusters identified: {results['summary']['total_clusters']}")
    print(f"\nOutput files saved to: {args.output}/")
    print(f"  - analysis_results.json          : Complete analysis data")
    print(f"  - eventalign_*_first_read.png    : Signal-sequence alignment visualization")
    print(f"  - dtw_heatmap_*.png              : Read-read DTW distance heatmaps")
    print(f"  - assessment_heatmap_*.png       : Eventalign quality heatmaps")
    print(f"  - assessment_multi_*.png         : Multi-metric quality assessment")
    print(f"  - transcript_reads_*.png         : Transcript structure with assigned reads")
    print(f"  - assignment_summary_*.png       : Read assignment summary statistics")


if __name__ == "__main__":
    main()
