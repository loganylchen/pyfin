"""
Integration matrix module for combining DTW and completeness metrics.

Provides functionality to:
1. Build completeness matrices (reads x isoforms)
2. Cluster reads by signal similarity using DTW
3. Integrate completeness scores with clusters to validate isoforms
4. Generate validation confidence scores
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional, Set
from dataclasses import dataclass
import logging
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform

from fin.core.isoform import IsoformSequence, ValidatedIsoform, IsoformEvidence
from fin.io.io_manager import ReadData
from fin.core.completeness import CompletenessResult
from fin.core.dtw_gpu import dtw_pairwise_matrix

logger = logging.getLogger(__name__)


@dataclass
class IntegrationMetrics:
    """
    Metrics from integrating completeness and DTW clustering.

    Attributes:
        isoform_id: ID of the isoform being evaluated
        read_support: Number of reads supporting this isoform
        avg_completeness: Mean completeness across supporting reads
        dtw_cluster_consistency: Proportion of reads in same cluster
        validation_score: Combined score (0-1) indicating validity
        confidence: Overall confidence (0-1) - higher is more confident
        is_validated: Whether isoform passes validation thresholds
    """
    isoform_id: str
    read_support: int = 0
    avg_completeness: float = 0.0
    dtw_cluster_consistency: float = 0.0
    validation_score: float = 0.0
    confidence: float = 0.0
    is_validated: bool = False


def build_completeness_matrix(
    reads: List[ReadData],
    isoforms: List[IsoformSequence],
    completeness_scores: Dict[Tuple[str, str], CompletenessResult]
) -> np.ndarray:
    """
    Build matrix of completeness scores (reads x isoforms).

    Args:
        reads: List of ReadData objects
        isoforms: List of IsoformSequence objects
        completeness_scores: Dictionary mapping (read_id, isoform_id) -> CompletenessResult

    Returns:
        N x M numpy array where N = num reads, M = num isoforms
        Each entry is the completeness score (0-1)
    """
    num_reads = len(reads)
    num_isoforms = len(isoforms)

    # Create mapping from IDs to indices
    read_id_to_idx = {read.read_id: i for i, read in enumerate(reads)}
    isoform_id_to_idx = {isoform.isoform_id: i for i, isoform in enumerate(isoforms)}

    # Initialize matrix
    matrix = np.zeros((num_reads, num_isoforms), dtype=np.float32)

    # Fill in scores
    for (read_id, isoform_id), comp_result in completeness_scores.items():
        if read_id in read_id_to_idx and isoform_id in isoform_id_to_idx:
            read_idx = read_id_to_idx[read_id]
            isoform_idx = isoform_id_to_idx[isoform_id]
            matrix[read_idx, isoform_idx] = comp_result.completeness_score

    return matrix


def cluster_reads_by_signal(
    read_signals: List[np.ndarray],
    dtw_threshold: float,
    use_cuda: bool = True,
    min_cluster_size: int = 2
) -> np.ndarray:
    """
    Cluster reads based on signal similarity using DTW distances.

    Uses connected components on a graph where edges connect reads with
    DTW distance below threshold.

    Args:
        read_signals: List of signal arrays for each read
        dtw_threshold: Maximum DTW distance for two reads to be connected
        use_cuda: Use GPU acceleration for DTW calculation
        min_cluster_size: Minimum reads per cluster to be considered valid

    Returns:
        Array of cluster labels (int) for each read
    """
    if not read_signals:
        logger.warning("No read signals provided for clustering")
        return np.array([], dtype=int)

    if len(read_signals) == 1:
        # Single read -> own cluster
        return np.array([0], dtype=int)

    # Compute pairwise DTW distance matrix
    logger.info(f"Computing pairwise DTW distances for {len(read_signals)} reads...")

    try:
        distance_matrix = dtw_pairwise_matrix(
            read_signals,
            open_start=False,
            open_end=False,
            use_cuda=use_cuda
        )

        if distance_matrix is None or np.isnan(distance_matrix).all():
            logger.warning("DTW distance matrix computation failed")
            return np.arange(len(read_signals))

        # Normalize distances by median (robust to outliers)
        nonzero_distances = distance_matrix[distance_matrix > 0]
        if len(nonzero_distances) > 0:
            median_dist = np.median(nonzero_distances)
            if median_dist > 0:
                normalized_distances = distance_matrix / median_dist
            else:
                normalized_distances = distance_matrix
        else:
            normalized_distances = distance_matrix

        # Create adjacency matrix: connect if distance < threshold
        # Use 1 - normalized_distance as similarity
        adjacency_matrix = csr_matrix(normalized_distances < dtw_threshold, dtype=bool)

        # Remove self-loops
        adjacency_matrix.setdiag(False)

        # Identify connected components
        n_components, labels = connected_components(
            adjacency_matrix,
            directed=False,
            return_labels=True
        )

        logger.info(f"Found {n_components} clusters from {len(read_signals)} reads")

        # Filter clusters by size (optional)
        if min_cluster_size > 1:
            cluster_counts = np.bincount(labels)
            small_clusters = np.where(cluster_counts < min_cluster_size)[0]

            if len(small_clusters) > 0:
                # Reassign small clusters to noise (-1)
                for i, label in enumerate(labels):
                    if label in small_clusters:
                        labels[i] = -1

                # Renumber remaining clusters
                unique_labels = sorted(set(labels) - {-1})
                label_mapping = {old: new for new, old in enumerate(unique_labels)}
                for i, label in enumerate(labels):
                    if label != -1:
                        labels[i] = label_mapping[label]

        return labels

    except Exception as e:
        logger.error(f"DTW clustering failed: {e}")
        # Fallback: each read is its own cluster
        return np.arange(len(read_signals))


def validate_isoforms(
    completeness_matrix: np.ndarray,
    dtw_cluster_labels: np.ndarray,
    min_read_support: int,
    min_completeness: float,
    min_cluster_consistency: float = 0.6
) -> Dict[int, List[int]]:
    """
    Validate which isoforms are real based on completeness and cluster consistency.

    Validation criteria:
    1. Minimum number of reads with high completeness
    2. Reads must cluster together (same or adjacent DTW clusters)
    3. Cluster consistency above threshold

    Args:
        completeness_matrix: N x M matrix (reads x isoforms)
        dtw_cluster_labels: Cluster assignment for each read
        min_read_support: Minimum reads needed to validate
        min_completeness: Minimum completeness score per read
        min_cluster_consistency: Minimum proportion in same cluster

    Returns:
        Dict mapping isoform_index -> list of read_indices that support it
    """
    if completeness_matrix.shape[0] == 0:
        return {}

    num_reads, num_isoforms = completeness_matrix.shape
    validated_isoforms = {}

    for isoform_idx in range(num_isoforms):
        # Find reads with high completeness for this isoform
        qualifying_reads = []
        completeness_scores = []

        for read_idx in range(num_reads):
            score = completeness_matrix[read_idx, isoform_idx]
            if score >= min_completeness:
                qualifying_reads.append(read_idx)
                completeness_scores.append(score)

        # Check minimum support
        if len(qualifying_reads) < min_read_support:
            continue

        # Check if reads cluster together
        if len(dtw_cluster_labels) > 0:
            # Get cluster labels for qualifying reads
            read_clusters = [dtw_cluster_labels[i] for i in qualifying_reads]

            # Count most common cluster
            cluster_counts = {}
            for cluster in read_clusters:
                cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1

            if not cluster_counts:
                continue

            most_common = max(cluster_counts, key=cluster_counts.get)
            cluster_consistency = cluster_counts[most_common] / len(qualifying_reads)

            # Check consistency threshold
            if cluster_consistency < min_cluster_consistency:
                continue

            # Also check that reads don't have high completeness to multiple isoforms
            # (indicates ambiguity)
            ambiguity_check = True
            for read_idx in qualifying_reads:
                max_completeness = np.max(completeness_matrix[read_idx])
                if max_completeness > 0.9:
                    # Check if this isoform is the best
                    best_isoform = np.argmax(completeness_matrix[read_idx])
                    if best_isoform != isoform_idx:
                        ambiguity_check = False
                        break

            if not ambiguity_check:
                continue

        # If all checks pass, add to validated
        validated_isoforms[isoform_idx] = qualifying_reads

        logger.debug(
            f"Isoform {isoform_idx}: {len(qualifying_reads)} reads, "
            f"cluster consistency: {cluster_consistency:.2f}"
        )

    return validated_isoforms


def calculate_integration_metrics(
    isoform_idx: int,
    supporting_reads: List[int],
    completeness_matrix: np.ndarray,
    dtw_cluster_labels: np.ndarray
) -> IntegrationMetrics:
    """
    Calculate integration metrics for a validated isoform.

    Args:
        isoform_idx: Index of the isoform
        supporting_reads: List of read indices that support this isoform
        completeness_matrix: N x M completeness matrix
        dtw_cluster_labels: Cluster assignments for each read

    Returns:
        IntegrationMetrics object
    """
    if not supporting_reads:
        return IntegrationMetrics(isoform_id=str(isoform_idx))

    metrics = IntegrationMetrics(isoform_id=str(isoform_idx))
    metrics.read_support = len(supporting_reads)

    # Calculate average completeness
    completeness_scores = [
        completeness_matrix[read_idx, isoform_idx]
        for read_idx in supporting_reads
    ]
    metrics.avg_completeness = np.mean(completeness_scores)

    # Calculate cluster consistency
    if len(dtw_cluster_labels) > 0:
        read_clusters = [dtw_cluster_labels[i] for i in supporting_reads]
        most_common_cluster = max(set(read_clusters), key=read_clusters.count)
        metrics.dtw_cluster_consistency = read_clusters.count(most_common_cluster) / len(supporting_reads)

    # Calculate validation score (weighted combination)
    # High completeness + high cluster consistency = high validation
    metrics.validation_score = (
        0.6 * metrics.avg_completeness +
        0.4 * metrics.dtw_cluster_consistency
    )

    # Calculate confidence score (read support weighted)
    metrics.confidence = metrics.validation_score * (1 - np.exp(-metrics.read_support / 5.0))

    # Determine if validated (all criteria met)
    metrics.is_validated = (
        metrics.read_support >= 5 and
        metrics.avg_completeness >= 0.8 and
        metrics.dtw_cluster_consistency >= 0.6
    )

    return metrics


def select_best_isoform_per_read(
    completeness_matrix: np.ndarray,
    validated_isoforms: Dict[int, List[int]]
) -> Dict[int, int]:
    """
    For each read, select the validated isoform with highest completeness.

    Args:
        completeness_matrix: N x M completeness matrix
        validated_isoforms: Dict mapping isoform_idx -> list of read_indices

    Returns:
        Dict mapping read_idx -> best_isoform_idx
    """
    read_assignments = {}

    # For each read, find all validated isoforms it supports
    read_to_isoforms = {}
    for isoform_idx, read_indices in validated_isoforms.items():
        for read_idx in read_indices:
            if read_idx not in read_to_isoforms:
                read_to_isoforms[read_idx] = []
            read_to_isoforms[read_idx].append(isoform_idx)

    # For reads with multiple isoforms, choose the best
    for read_idx, isoform_indices in read_to_isoforms.items():
        if len(isoform_indices) == 1:
            read_assignments[read_idx] = isoform_indices[0]
        else:
            # Choose isoform with highest completeness
            best_isoform = max(
                isoform_indices,
                key=lambda i: completeness_matrix[read_idx, i]
            )
            read_assignments[read_idx] = best_isoform

    return read_assignments


def build_validation_summary(
    validated_isoforms: Dict[int, List[int]],
    completeness_matrix: np.ndarray,
    dtw_cluster_labels: np.ndarray,
    isoform_sequences: List[IsoformSequence]
) -> List[Dict[str, Any]]:
    """
    Build summary statistics for validated isoforms.

    Args:
        validated_isoforms: ISO index -> read indices dict
        completeness_matrix: N x M completeness matrix
        dtw_cluster_labels: DTW cluster assignments
        isoform_sequences: List of isoform sequences

    Returns:
        List of dictionaries with summary statistics
    """
    summaries = []

    for isoform_idx, read_indices in validated_isoforms.items():
        if isoform_idx >= len(isoform_sequences):
            continue

        isoform = isoform_sequences[isoform_idx]
        metrics = calculate_integration_metrics(
            isoform_idx, read_indices, completeness_matrix, dtw_cluster_labels
        )

        summary = {
            "isoform_id": isoform.isoform_id,
            "gene_id": isoform.gene_id,
            "read_support": metrics.read_support,
            "avg_completeness": f"{metrics.avg_completeness:.3f}",
            "cluster_consistency": f"{metrics.dtw_cluster_consistency:.3f}",
            "validation_score": f"{metrics.validation_score:.3f}",
            "confidence": f"{metrics.confidence:.3f}",
            "is_validated": metrics.is_validated,
            "num_exons": isoform.num_exons,
            "length": isoform.length
        }
        summaries.append(summary)

    # Sort by confidence
    summaries.sort(key=lambda x: float(x["confidence"]), reverse=True)

    return summaries
