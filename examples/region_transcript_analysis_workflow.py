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
        self, chrom: str, start: int, end: int
    ) -> List[Tuple[str, str]]:
        """
        Get all candidate transcript sequences that overlap a region.

        Args:
            chrom: Chromosome name
            start: Region start (0-based)
            end: Region end

        Returns:
            List of (transcript_id, sequence) tuples
        """
        candidates = []

        with GTFReader(str(self.gtf_path)) as reader:
            reader.parse()

            # Get transcripts overlapping this region
            transcripts = reader.get_transcripts_in_region(chrom, start, end)

            for tx in transcripts:
                tx_id = tx.transcript_id
                if tx_id in self.transcriptome_sequences:
                    candidates.append((tx_id, self.transcriptome_sequences[tx_id]))

        return candidates

    def get_reads_for_region(self, chrom: str, start: int, end: int) -> List[Dict[str, Any]]:
        """
        Get all reads mapping to a region.

        Args:
            chrom: Chromosome name
            start: Region start (0-based)
            end: Region end

        Returns:
            List of read dictionaries
        """
        reads = []

        with BamReader(str(self.bam_path)) as reader:
            for read in reader.fetch(region=f"{chrom}:{start}-{end}"):
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
                scaled_model_mean = aln.get("scaled_model_mean", aln.get("model_mean", 0))
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

                # Draw model expected level
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
                    label_y = max(event_mean, scaled_model_mean) + 3
                    ax.text(
                        (start_sample + end_sample) / 2,
                        label_y,
                        f"{kmer}\nidx:{event_idx}\npos:{ref_pos}",
                        fontsize=6,
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
                    color="limegreen",
                    linewidth=2,
                    linestyle="--",
                    label="Model mean (expected)",
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

        # Step 1: Get candidate transcripts for this region
        candidates = self.get_candidate_transcripts_for_region(
            region.chrom, region.start, region.end
        )

        if len(candidates) > max_transcripts:
            logger.info(f"  Limiting to {max_transcripts} transcripts (of {len(candidates)})")
            candidates = candidates[:max_transcripts]

        result["candidate_transcripts"] = [
            {"transcript_id": tx_id, "length": len(seq)} for tx_id, seq in candidates
        ]
        logger.info(f"  Found {len(candidates)} candidate transcripts")

        if not candidates:
            logger.warning(f"  No candidate transcripts found for region {region_id}")
            return result

        # Step 2: Get reads for this region
        reads = self.get_reads_for_region(region.chrom, region.start, region.end)

        if len(reads) > max_reads:
            logger.info(f"  Limiting to {max_reads} reads (of {len(reads)})")
            reads = reads[:max_reads]

        logger.info(f"  Processing {len(reads)} reads")

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
            },
            "regions": self.region_results,
        }

        # Save results
        output_file = self.output_dir / "analysis_results.json"
        with open(output_file, "w") as f:
            json.dump(final_results, f, indent=2)
        logger.info(f"\nResults saved to: {output_file}")

        return final_results


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
                f"(state={aln['hmm_state']}, event_mean={aln['event_mean']:.1f})"
            )
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
    print(f"\nResults saved to: {args.output}/analysis_results.json")
    print(f"Eventalign plots: {args.output}/eventalign_*_first_read.png")
    print(f"DTW heatmaps: {args.output}/dtw_heatmap_*.png")


if __name__ == "__main__":
    main()
