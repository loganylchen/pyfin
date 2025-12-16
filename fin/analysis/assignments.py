import numpy as np


def em_with_coherence(
    dist_read_to_tx,
    dist_read_to_read,
    sigma=1.0,
    beta=0.5,
    max_iter=1000,
    tol=1e-4,
    min_sigma=0.01,
    verbose=True
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
    assert dist_read_to_read.shape == (n_reads, n_reads), \
        f"dist_read_to_read must be ({n_reads}, {n_reads}), got {dist_read_to_read.shape}"
    assert sigma >= min_sigma, f"sigma must be >= {min_sigma}"
    assert beta >= 0, "beta must be non-negative"
    assert max_iter > 0, "max_iter must be positive"
    
    # Ensure distance matrices are non-negative
    if dist_read_to_tx.min() < 0 or dist_read_to_read.min() < 0:
        raise ValueError("Distance matrices must be non-negative")
    
    # Initialize responsibility matrix: R[i, j] = P(read_i -> transcript_j)
    R = np.exp(-dist_read_to_tx / sigma)
    row_sums = R.sum(axis=1, keepdims=True)
    # Handle rows with all zeros (reads that don't match any transcript)
    row_sums = np.maximum(row_sums, 1e-10)
    R = R / row_sums
    
    log_likelihoods = []
    
    for it in range(max_iter):
        # E-step: Compute coherence penalty (vectorized for efficiency)
        # For each transcript j, coherence[i,j] = weighted avg distance from read i to cluster j
        # coherence[i,j] = sum_k R[k,j] * dist[i,k] / sum_k R[k,j]
        #                = (dist_read_to_read[i,:] @ R[:,j]) / R[:,j].sum()
        
        # Shape: (n_reads, n_tx)
        # Matrix multiplication: dist_read_to_read @ R gives weighted distances
        coherence_numerator = dist_read_to_read @ R  # (n_reads, n_tx)
        
        # Sum of responsibilities per transcript (cluster sizes)
        cluster_weights = R.sum(axis=0, keepdims=True)  # (1, n_tx)
        cluster_weights = np.maximum(cluster_weights, 1e-10)  # Avoid division by zero
        
        # Compute coherence penalty
        coherence_penalty = coherence_numerator / cluster_weights  # (n_reads, n_tx)
        
        # M-step: Update responsibilities with coherence-aware energy
        energy = dist_read_to_tx + beta * coherence_penalty
        
        # Compute new responsibilities (softmax with temperature sigma)
        R_new = np.exp(-energy / sigma)
        row_sums = R_new.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-10)
        R_new = R_new / row_sums
        
        # Check convergence
        diff = np.abs(R - R_new).max()
        
        # Compute log-likelihood (negative energy weighted by responsibilities)
        # This is -E[energy] where expectation is over R
        ll = -np.sum(R_new * energy)
        log_likelihoods.append(ll)
        
        # Update responsibilities
        R = R_new
        
        if verbose and (it % 100 == 0 or diff < tol):
            print(f"Iter {it:4d}: log-likelihood = {ll:12.4f}, max_diff = {diff:.6f}")
        
        if diff < tol:
            if verbose:
                print(f"✓ EM converged at iteration {it}")
            break
    else:
        if verbose:
            print(f"⚠ EM reached max iterations ({max_iter}) without converging")
    
    # Compute hard assignments
    hard_assignments = np.argmax(R, axis=1)
    
    # Print cluster statistics
    if verbose:
        unique_assignments, counts = np.unique(hard_assignments, return_counts=True)
        print(f"\nCluster statistics:")
        print(f"  Active transcripts: {len(unique_assignments)} / {n_tx}")
        print(f"  Reads per cluster: min={counts.min()}, max={counts.max()}, "
              f"mean={counts.mean():.1f}, median={np.median(counts):.1f}")
        
        # Compute average confidence
        avg_confidence = np.mean([R[i, hard_assignments[i]] for i in range(n_reads)])
        print(f"  Average assignment confidence: {avg_confidence:.3f}")
    
    return R, hard_assignments, log_likelihoods
