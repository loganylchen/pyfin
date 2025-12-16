def em_with_coherence(
	dist_read_to_tx,  		 # (n_reads, n_tx)
    dist_read_to_read,       # (n_reads, n_reads)
    sigma=1.0,               # scale for read-tx distance
    beta=0.5,                # weight for coherence term (higher = stronger clustering)
    max_iter=1000,
    tol=1e-4):
    n_reads, n_tx = dist_read_to_tx.shape
    # Initialize responsibility matrix: R[i, j] = P(z_i = j | data)
    R = np.exp(-dist_read_to_tx / sigma)
    R /= R.sum(axis=1, keepdims=True)  # normalize
    log_likelihoods = []
    for it in range(max_iter):
        # E-step: already in R
        # Compute current "cluster coherence": for each transcript j,
        # compute expected average pairwise distance among reads assigned to j
        coherence_penalty = np.zeros((n_reads, n_tx))
        for j in range(n_tx):
            if R[:, j].sum() < 1e-6:
                continue
            # Weighted average distance from read i to all others in cluster j
            for i in range(n_reads):
                # E[ D(r_i, r_k) | both in j ] ≈ sum_k R[k,j] * D(i,k) / sum_k R[k,j]
                weighted_dist = np.dot(R[:, j], dist_read_to_read[i, :])
                total_weight = R[:, j].sum()
                if total_weight > 0:
                    coherence_penalty[i, j] = weighted_dist / total_weight
        # M-step: update responsibilities with coherence-aware energy
        energy = dist_read_to_tx + beta * coherence_penalty
        R_new = np.exp(-energy / sigma)
        R_new /= R_new.sum(axis=1, keepdims=True)
        # Check convergence
        diff = np.abs(R - R_new).max()
        R = R_new
        # Log-likelihood (optional)
        ll = -np.sum(R * energy)
        log_likelihoods.append(ll)
        if diff < tol:
            print(f"EM converged at iter {it}")
            break
    # Return soft assignments (R) and hard assignments
    hard_assignments = np.argmax(R, axis=1)
    return R, hard_assignments, log_likelihoods