#!/usr/bin/env python3
"""
Example: Batch pairwise DTW distance computation

This demonstrates the GPU-accelerated batch pairwise DTW function,
which is ideal for computing distance matrices for clustering,
classification, or similarity search tasks.
"""

import numpy as np
import time
from fin._dtw import dtw_pairwise, dtw_distance, is_available, cleanup

def main():
    if not is_available():
        print("ERROR: CUDA DTW is not available")
        print("Please rebuild the package with CUDA support")
        return
    
    print("="*80)
    print("Batch Pairwise DTW Example")
    print("="*80)
    
    # Example 1: Small batch
    print("\n1. Computing pairwise distances for 5 sequences")
    print("-"*80)
    
    num_sequences = 5
    seq_length = 100
    
    # Generate random sequences
    np.random.seed(42)
    sequences = np.random.randn(num_sequences, seq_length).astype(np.float32)
    
    print(f"Sequences shape: {sequences.shape}")
    
    # Compute pairwise distances
    start = time.perf_counter()
    distance_matrix = dtw_pairwise(sequences)
    elapsed_ms = (time.perf_counter() - start) * 1000
    
    print(f"Distance matrix shape: {distance_matrix.shape}")
    print(f"Computation time: {elapsed_ms:.2f} ms")
    print("\nDistance matrix:")
    print(distance_matrix)
    
    # Verify: diagonal should be all zeros
    print(f"\nDiagonal (should be all zeros): {np.diag(distance_matrix)}")
    
    # Verify: matrix should be symmetric
    is_symmetric = np.allclose(distance_matrix, distance_matrix.T)
    print(f"Matrix is symmetric: {is_symmetric}")
    
    # Example 2: Compare batch vs individual computation
    print("\n\n2. Comparing batch vs individual computation")
    print("-"*80)
    
    num_sequences = 10
    seq_length = 200
    sequences = np.random.randn(num_sequences, seq_length).astype(np.float32)
    
    # Batch computation
    print(f"\nBatch computation ({num_sequences} sequences, length {seq_length}):")
    start = time.perf_counter()
    batch_distances = dtw_pairwise(sequences)
    batch_time_ms = (time.perf_counter() - start) * 1000
    num_pairs = (num_sequences * (num_sequences - 1)) // 2
    print(f"  Total time: {batch_time_ms:.2f} ms")
    print(f"  Time per pair: {batch_time_ms / num_pairs:.2f} ms ({num_pairs} pairs)")
    
    # Individual computation
    print("\nIndividual computation (for comparison):")
    individual_distances = np.zeros((num_sequences, num_sequences), dtype=np.float32)
    
    start = time.perf_counter()
    for i in range(num_sequences):
        for j in range(i+1, num_sequences):
            dist = dtw_distance(sequences[i], sequences[j])
            individual_distances[i, j] = dist
            individual_distances[j, i] = dist
    individual_time_ms = (time.perf_counter() - start) * 1000
    
    print(f"  Total time: {individual_time_ms:.2f} ms")
    print(f"  Time per pair: {individual_time_ms / num_pairs:.2f} ms")
    
    speedup = individual_time_ms / batch_time_ms
    print(f"\n  Speedup: {speedup:.2f}x faster with batch API")
    
    # Verify results match
    max_diff = np.max(np.abs(batch_distances - individual_distances))
    print(f"  Max difference: {max_diff:.6f}")
    if max_diff < 0.01:
        print("  ✓ Results match!")
    else:
        print("  ⚠ Results differ (may need investigation)")
    
    # Example 3: Larger batch to show GPU advantage
    print("\n\n3. Large batch computation (GPU advantage)")
    print("-"*80)
    
    num_sequences = 50
    seq_length = 500
    sequences = np.random.randn(num_sequences, seq_length).astype(np.float32)
    num_pairs = (num_sequences * (num_sequences - 1)) // 2
    
    print(f"Computing {num_pairs} pairwise distances")
    print(f"Sequences: {num_sequences} × {seq_length} samples")
    
    start = time.perf_counter()
    distance_matrix = dtw_pairwise(sequences)
    elapsed_ms = (time.perf_counter() - start) * 1000
    
    print(f"\nBatch computation time: {elapsed_ms:.2f} ms")
    print(f"Average time per pair: {elapsed_ms / num_pairs:.2f} ms")
    print(f"Throughput: {num_pairs / (elapsed_ms / 1000):.1f} pairs/second")
    
    # Example 4: Using distance matrix for clustering
    print("\n\n4. Using distance matrix for clustering")
    print("-"*80)
    
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
        from scipy.spatial.distance import squareform
        
        num_sequences = 20
        seq_length = 100
        
        # Create sequences with some structure (3 clusters)
        np.random.seed(42)
        cluster1 = np.random.randn(7, seq_length) + 0.0
        cluster2 = np.random.randn(7, seq_length) + 3.0
        cluster3 = np.random.randn(6, seq_length) - 3.0
        sequences = np.vstack([cluster1, cluster2, cluster3]).astype(np.float32)
        
        # Compute pairwise distances
        distance_matrix = dtw_pairwise(sequences)
        
        # Convert to condensed distance matrix for scipy
        condensed_distances = squareform(distance_matrix)
        
        # Hierarchical clustering
        linkage_matrix = linkage(condensed_distances, method='average')
        
        # Get cluster assignments
        clusters = fcluster(linkage_matrix, t=3, criterion='maxclust')
        
        print("Cluster assignments:")
        print(f"  Cluster 1 (expected 0-6):   {np.where(clusters == 1)[0].tolist()}")
        print(f"  Cluster 2 (expected 7-13):  {np.where(clusters == 2)[0].tolist()}")
        print(f"  Cluster 3 (expected 14-19): {np.where(clusters == 3)[0].tolist()}")
        
    except ImportError:
        print("scipy not available - skipping clustering example")
        print("Install with: pip install scipy")
    
    # Example 5: Memory efficiency
    print("\n\n5. Memory efficiency for large batches")
    print("-"*80)
    
    num_sequences = 100
    seq_length = 1000
    sequences = np.random.randn(num_sequences, seq_length).astype(np.float32)
    num_pairs = (num_sequences * (num_sequences - 1)) // 2
    
    print(f"Computing {num_pairs} pairwise distances")
    print(f"Input size: {sequences.nbytes / 1024 / 1024:.2f} MB")
    print(f"Output size: {num_sequences * num_sequences * 4 / 1024 / 1024:.2f} MB")
    
    start = time.perf_counter()
    distance_matrix = dtw_pairwise(sequences)
    elapsed_ms = (time.perf_counter() - start) * 1000
    
    print(f"\nComputation time: {elapsed_ms:.2f} ms")
    print(f"Throughput: {num_pairs / (elapsed_ms / 1000):.1f} pairs/second")
    
    # Cleanup
    cleanup()
    print("\n✓ CUDA resources cleaned up")
    print("\n" + "="*80)
    print("Examples complete!")
    print("="*80)

if __name__ == "__main__":
    main()
