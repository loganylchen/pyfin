"""EM algorithm with coherence regularization for read-to-transcript assignment.

Supports CuPy GPU acceleration when available, with automatic numpy fallback.
"""

import logging
import numpy as np
from typing import Optional

logger = logging.getLogger(__name__)

# Try CuPy for GPU-accelerated matrix ops
try:
    import cupy as cp

    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False


def _get_array_module(use_gpu=True):
    """Return (xp, to_device, to_host) tuple for numpy or cupy."""
    if use_gpu and CUPY_AVAILABLE:
        return cp, cp.asarray, cp.asnumpy
    return np, lambda x: x, lambda x: x


def em_with_coherence(
    dist_read_to_tx,
    dist_read_to_read,
    sigma=1.0,
    beta=0.5,
    max_iter=1000,
    tol=1e-4,
    min_sigma=0.01,
    verbose=True,
    use_gpu=True,
    prior_weights: Optional[np.ndarray] = None,
):
    """
    EM algorithm with coherence regularization for clustering reads to transcripts.

    This algorithm assigns reads to transcripts while encouraging reads assigned
    to the same transcript to be similar to each other (coherence). Combines:
    1. Read-to-transcript distance (alignment quality)
    2. Read-to-read coherence within each transcript cluster

    Args:
        dist_read_to_tx: np.ndarray, shape (n_reads, n_tx)
            Distance matrix from reads to transcripts (lower = better alignment)
        dist_read_to_read: np.ndarray, shape (n_reads, n_reads)
            Pairwise distance matrix between reads (e.g., DTW, edit distance)
        sigma: float, default=1.0
            Temperature parameter for softmax. Lower = harder assignments.
        beta: float, default=0.5
            Weight for coherence term. Higher = stronger clustering pressure.
            - beta=0: ignore read-read similarity, pure read-tx assignment
            - beta>0: encourage similar reads to cluster together
        max_iter: int, default=1000
            Maximum number of EM iterations
        tol: float, default=1e-4
            Convergence tolerance (max absolute change in responsibilities)
        min_sigma: float, default=0.01
            Minimum allowed sigma to prevent numerical issues
        verbose: bool, default=True
            Print convergence information
        use_gpu: bool, default=True
            Use CuPy GPU acceleration when available
        prior_weights: np.ndarray or None, shape (n_tx,), default=None
            Optional prior weights over transcripts. When provided, the initial
            responsibility matrix is scaled by the normalized prior before row
            normalization:  R_init[i,j] = exp(-d_tx[i,j]/sigma) * pi[j], where
            pi = prior_weights / prior_weights.sum(). Must be non-negative and
            finite. If all weights are zero, falls back to uniform prior with a
            warning. When None (default), behavior is identical to the original
            (uniform) initialization.

    Returns:
        tuple: (R, hard_assignments, log_likelihoods)
            - R: np.ndarray, shape (n_reads, n_tx)
                Soft assignment probabilities (responsibility matrix)
                R[i, j] = P(read_i belongs to transcript_j | data)
            - hard_assignments: np.ndarray, shape (n_reads,)
                Hard cluster assignments (argmax of R)
            - log_likelihoods: list
                Log-likelihood at each iteration (for convergence diagnostics)

    Algorithm:
        E-step: Compute coherence penalty for each (read, transcript) pair
                coherence[i,j] = expected distance from read i to other reads in cluster j
        M-step: Update responsibilities using:
                energy = dist_read_to_tx + beta * coherence
                R = softmax(-energy / sigma)

    Example:
        >>> # Assign reads to transcripts with coherence
        >>> R, assignments, lls = em_with_coherence(
        ...     dist_read_to_tx=dtw_distances,
        ...     dist_read_to_read=read_similarity,
        ...     sigma=1.0,
        ...     beta=0.5
        ... )
        >>> # Get transcript assignment for each read
        >>> for read_idx, tx_idx in enumerate(assignments):
        ...     confidence = R[read_idx, tx_idx]
        ...     print(f"Read {read_idx} -> Transcript {tx_idx} (p={confidence:.3f})")
    """
    # Input validation
    n_reads, n_tx = dist_read_to_tx.shape
    if dist_read_to_read.shape != (n_reads, n_reads):
        raise ValueError(
            f"dist_read_to_read must be ({n_reads}, {n_reads}), "
            f"got {dist_read_to_read.shape}"
        )
    if sigma < min_sigma:
        raise ValueError(f"sigma must be >= {min_sigma}, got {sigma}")
    if beta < 0:
        raise ValueError(f"beta must be non-negative, got {beta}")
    if max_iter < 0:
        raise ValueError(f"max_iter must be non-negative, got {max_iter}")

    # Ensure distance matrices are non-negative
    if dist_read_to_tx.min() < 0 or dist_read_to_read.min() < 0:
        raise ValueError("Distance matrices must be non-negative")

    # Validate and normalize prior_weights
    pi = None
    if prior_weights is not None:
        prior_weights = np.asarray(prior_weights, dtype=float)
        if prior_weights.shape != (n_tx,):
            raise ValueError(
                f"prior_weights must be ({n_tx},), got {prior_weights.shape}"
            )
        if np.any(prior_weights < 0) or not np.all(np.isfinite(prior_weights)):
            raise ValueError("prior_weights must be non-negative and finite")
        s = prior_weights.sum()
        if s <= 0:
            if verbose:
                logger.warning("prior_weights sum to 0; falling back to uniform prior")
            pi = None
        else:
            pi = prior_weights / s

    # Select array backend (cupy GPU or numpy CPU)
    xp, to_device, to_host = _get_array_module(use_gpu)
    if xp is not np and verbose:
        logger.info("GPU (CuPy) acceleration enabled for EM")

    # Transfer to device
    d_tx = to_device(dist_read_to_tx)
    d_rr = to_device(dist_read_to_read)

    # Initialize responsibility matrix: R[i, j] = P(read_i -> transcript_j)
    # Numerically stable softmax: subtract row-min before exp so the best
    # distance maps to exp(0)=1. The constant factor cancels in row
    # normalization, so this is mathematically identical to exp(-d/σ).
    # Without this shift, eventalign distances of ~10^3 plus 1e6 missing-pair
    # fallbacks underflow to 0 → R is all zeros → abundance is all zeros.
    d_tx_min = d_tx.min(axis=1, keepdims=True)
    R = xp.exp(-(d_tx - d_tx_min) / sigma)
    if pi is not None:
        pi_dev = to_device(pi.reshape(1, -1))  # shape (1, n_tx) for broadcasting
        R = R * pi_dev
    row_sums = R.sum(axis=1, keepdims=True)
    row_sums = xp.maximum(row_sums, 1e-10)
    R = R / row_sums

    log_likelihoods = []

    for it in range(max_iter):
        # E-step: coherence[i,j] = (dist_read_to_read @ R)[:,j] / R[:,j].sum()
        coherence_numerator = d_rr @ R  # (n_reads, n_tx) — main GPU bottleneck
        cluster_weights = R.sum(axis=0, keepdims=True)  # (1, n_tx)
        cluster_weights = xp.maximum(cluster_weights, 1e-10)
        coherence_penalty = coherence_numerator / cluster_weights

        # M-step: energy = dist + beta * coherence, then softmax
        # Row-shift for numerical stability (see init step above for rationale).
        energy = d_tx + beta * coherence_penalty
        energy_min = energy.min(axis=1, keepdims=True)
        R_new = xp.exp(-(energy - energy_min) / sigma)
        row_sums = R_new.sum(axis=1, keepdims=True)
        row_sums = xp.maximum(row_sums, 1e-10)
        R_new = R_new / row_sums

        # Convergence check
        diff = float(xp.abs(R - R_new).max())

        # Free energy / ELBO surrogate: -<E>_R + sigma * H(R).
        # Without the entropy term the "log-likelihood" trace is just the
        # expected energy, so it never reflects soft-assignment uncertainty
        # and convergence diagnostics under-report when EM is doing useful
        # uncertainty smoothing. Sigma scales the entropy back into the same
        # units as the energy term (since E/sigma is what is exponentiated).
        eps = 1e-12
        entropy = float(-xp.sum(R_new * xp.log(R_new + eps)))
        ll = float(-xp.sum(R_new * energy)) + sigma * entropy
        log_likelihoods.append(ll)

        R = R_new

        if verbose and (it % 100 == 0 or diff < tol):
            logger.debug(
                "Iter %4d: log-likelihood = %12.4f, max_diff = %.6f", it, ll, diff
            )

        if diff < tol:
            if verbose:
                logger.info("EM converged at iteration %d", it)
            break
    else:
        if verbose:
            logger.info("EM reached max iterations (%d) without converging", max_iter)

    # Transfer back to host (numpy)
    R_host = to_host(R)
    hard_assignments = np.argmax(R_host, axis=1)

    if verbose:
        unique_assignments, counts = np.unique(hard_assignments, return_counts=True)
        avg_confidence = R_host[np.arange(n_reads), hard_assignments].mean()
        logger.info(
            "EM cluster stats: active_tx=%d/%d, reads_per_cluster min=%d max=%d "
            "mean=%.1f median=%.1f, avg_confidence=%.3f",
            len(unique_assignments), n_tx,
            counts.min(), counts.max(),
            counts.mean(), float(np.median(counts)),
            avg_confidence,
        )

    return R_host, hard_assignments, log_likelihoods
