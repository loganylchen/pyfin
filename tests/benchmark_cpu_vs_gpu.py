#!/usr/bin/env python3
"""
Benchmark: CPU vs GPU Performance for Eventalign and DTW

This script benchmarks the performance of:
1. Event alignment (eventalign) - CPU vs GPU (CUDA)
2. DTW pairwise distance - CPU fallback vs CUDA

Requirements:
    pip install numpy matplotlib tabulate

Build Options:
    - CPU only: pip install -e .
    - GPU (CUDA): CUDA_HOME=/usr/local/cuda pip install -e .

Usage:
    python benchmark_cpu_vs_gpu.py [--signal-sizes 1000,5000,10000]
                                   [--n-sequences 5,10,20]
                                   [--iterations 3]
"""

import argparse
import time
import sys
import numpy as np

# Check what's available - CPU eventalign
try:
    from fin._f5c import (
        profile_hmm_eventalign,
        eventalign_is_available,
    )

    EVENTALIGN_CPU_AVAILABLE = eventalign_is_available()
except ImportError:
    EVENTALIGN_CPU_AVAILABLE = False

# Check what's available - CUDA eventalign
try:
    from fin._f5c import (
        profile_hmm_eventalign_cuda,
        eventalign_cuda_is_available,
    )

    EVENTALIGN_CUDA_AVAILABLE = eventalign_cuda_is_available()
except ImportError:
    EVENTALIGN_CUDA_AVAILABLE = False

# Check DTW CUDA
try:
    from fin._dtw import (
        dtw_pairwise,
        is_available as dtw_is_available,
    )

    DTW_CUDA_AVAILABLE = dtw_is_available()
except ImportError:
    DTW_CUDA_AVAILABLE = False

# Try to import tabulate for nice tables
try:
    from tabulate import tabulate

    TABULATE_AVAILABLE = True
except ImportError:
    TABULATE_AVAILABLE = False


def print_table(headers, rows, title=None):
    """Print a formatted table."""
    if title:
        print(f"\n{'='*60}")
        print(f" {title}")
        print(f"{'='*60}")

    if TABULATE_AVAILABLE:
        print(tabulate(rows, headers=headers, tablefmt="grid"))
    else:
        # Simple fallback
        col_widths = [
            max(len(str(h)), max(len(str(row[i])) for row in rows)) for i, h in enumerate(headers)
        ]
        header_line = " | ".join(f"{h:<{col_widths[i]}}" for i, h in enumerate(headers))
        print(header_line)
        print("-" * len(header_line))
        for row in rows:
            print(" | ".join(f"{str(v):<{col_widths[i]}}" for i, v in enumerate(row)))


def generate_random_sequence(length: int) -> str:
    """Generate a random RNA sequence."""
    bases = "ACGU"  # RNA-only mode
    return "".join(np.random.choice(list(bases), size=length))


def generate_synthetic_signal(
    sequence: str, samples_per_base: int = 15, noise_level: float = 0.5
) -> np.ndarray:
    """
    Generate synthetic nanopore-like signal from a sequence.

    This creates a signal with:
    - Different current levels for each base
    - Gaussian noise
    - Some level variation per event
    """
    # Approximate current levels for each base (pA-like values)
    base_levels = {
        "A": 65.0,
        "C": 85.0,
        "G": 75.0,
        "T": 70.0,
        "U": 68.0,  # Similar to T for RNA
    }

    signal_parts = []
    for base in sequence:
        level = base_levels.get(base, 70.0)
        # Add variation per event
        n_samples = samples_per_base + np.random.randint(-3, 4)
        n_samples = max(5, n_samples)  # At least 5 samples
        event_level = level + np.random.randn() * 2.0  # Base-level variation
        event_signal = event_level + np.random.randn(n_samples) * noise_level
        signal_parts.append(event_signal)

    return np.concatenate(signal_parts).astype(np.float32)


def benchmark_eventalign_cpu(signal_sizes: list, iterations: int = 3) -> list:
    """
    Benchmark CPU eventalign performance across different signal sizes.

    Returns:
        List of (signal_size, seq_length, n_events, time_ms, events_per_ms) tuples
    """
    if not EVENTALIGN_CPU_AVAILABLE:
        print("WARNING: CPU Eventalign is not available - skipping benchmark")
        return []

    results = []

    for signal_size in signal_sizes:
        # Create sequence that will generate approximately this signal size
        seq_length = signal_size // 15  # ~15 samples per base
        seq_length = max(20, seq_length)  # At least 20 bases

        sequence = generate_random_sequence(seq_length)
        signal = generate_synthetic_signal(sequence)

        times = []
        n_events = 0

        for i in range(iterations):
            start = time.perf_counter()
            result = profile_hmm_eventalign(
                raw_signal=signal, sequence=sequence, kmer_size=5, events_per_base=3.0
            )
            elapsed = time.perf_counter() - start
            times.append(elapsed)
            n_events = result.get("n_events", 0)

        avg_time_ms = np.mean(times) * 1000
        events_per_ms = n_events / avg_time_ms if avg_time_ms > 0 else 0

        results.append(
            (len(signal), seq_length, n_events, f"{avg_time_ms:.2f}", f"{events_per_ms:.1f}")
        )

    return results


def benchmark_eventalign_cuda(signal_sizes: list, iterations: int = 3) -> list:
    """
    Benchmark CUDA (GPU) eventalign performance across different signal sizes.

    Returns:
        List of (signal_size, seq_length, n_events, time_ms, events_per_ms) tuples
    """
    if not EVENTALIGN_CUDA_AVAILABLE:
        print("WARNING: CUDA Eventalign is not available - skipping benchmark")
        return []

    results = []

    for signal_size in signal_sizes:
        # Create sequence that will generate approximately this signal size
        seq_length = signal_size // 15  # ~15 samples per base
        seq_length = max(20, seq_length)  # At least 20 bases

        sequence = generate_random_sequence(seq_length)
        signal = generate_synthetic_signal(sequence)

        times = []
        n_events = 0

        for i in range(iterations):
            start = time.perf_counter()
            result = profile_hmm_eventalign_cuda(
                raw_signal=signal, sequence=sequence, kmer_size=5, events_per_base=3.0
            )
            elapsed = time.perf_counter() - start
            times.append(elapsed)
            n_events = result.get("n_events", 0)

        avg_time_ms = np.mean(times) * 1000
        events_per_ms = n_events / avg_time_ms if avg_time_ms > 0 else 0

        results.append(
            (len(signal), seq_length, n_events, f"{avg_time_ms:.2f}", f"{events_per_ms:.1f}")
        )

    return results


def benchmark_dtw_pairwise_cuda(
    n_sequences_list: list, seq_length: int = 500, iterations: int = 3
) -> list:
    """
    Benchmark DTW pairwise distance computation with CUDA.

    Returns:
        List of (n_sequences, n_pairs, time_ms, pairs_per_ms) tuples
    """
    if not DTW_CUDA_AVAILABLE:
        print("WARNING: DTW CUDA is not available - skipping benchmark")
        return []

    results = []

    for n_seq in n_sequences_list:
        # Generate random sequences
        sequences = [np.random.randn(seq_length).astype(np.float32) for _ in range(n_seq)]
        n_pairs = n_seq * (n_seq - 1) // 2  # Number of unique pairs

        times = []

        for i in range(iterations):
            start = time.perf_counter()
            distance_matrix = dtw_pairwise(sequences)
            elapsed = time.perf_counter() - start
            times.append(elapsed)

        avg_time_ms = np.mean(times) * 1000
        pairs_per_ms = n_pairs / avg_time_ms if avg_time_ms > 0 else 0

        results.append((n_seq, n_pairs, f"{avg_time_ms:.2f}", f"{pairs_per_ms:.1f}"))

    return results


def benchmark_dtw_pairwise_cpu(
    n_sequences_list: list, seq_length: int = 500, iterations: int = 3
) -> list:
    """
    Benchmark DTW pairwise distance computation with pure NumPy (CPU).

    This is a reference CPU implementation for comparison.
    """

    def dtw_distance_numpy(s1: np.ndarray, s2: np.ndarray) -> float:
        """Simple DTW implementation using NumPy."""
        n, m = len(s1), len(s2)

        # Create cost matrix
        dtw = np.full((n + 1, m + 1), np.inf)
        dtw[0, 0] = 0

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                cost = abs(s1[i - 1] - s2[j - 1])
                dtw[i, j] = cost + min(dtw[i - 1, j], dtw[i, j - 1], dtw[i - 1, j - 1])

        return dtw[n, m]

    def dtw_pairwise_cpu(sequences: list) -> np.ndarray:
        """Compute pairwise DTW distances using CPU."""
        n = len(sequences)
        distances = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = dtw_distance_numpy(sequences[i], sequences[j])
                distances[i, j] = d
                distances[j, i] = d
        return distances

    results = []

    for n_seq in n_sequences_list:
        # Generate random sequences
        sequences = [np.random.randn(seq_length).astype(np.float32) for _ in range(n_seq)]
        n_pairs = n_seq * (n_seq - 1) // 2

        times = []

        for i in range(iterations):
            start = time.perf_counter()
            distance_matrix = dtw_pairwise_cpu(sequences)
            elapsed = time.perf_counter() - start
            times.append(elapsed)

        avg_time_ms = np.mean(times) * 1000
        pairs_per_ms = n_pairs / avg_time_ms if avg_time_ms > 0 else 0

        results.append((n_seq, n_pairs, f"{avg_time_ms:.2f}", f"{pairs_per_ms:.1f}"))

    return results


def print_system_info():
    """Print system and availability information."""
    print("\n" + "=" * 60)
    print(" SYSTEM INFORMATION")
    print("=" * 60)

    # Check CUDA availability
    cuda_info = "Not checked"
    try:
        import subprocess

        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            cuda_info = result.stdout.strip()
        else:
            cuda_info = "nvidia-smi failed"
    except Exception as e:
        cuda_info = f"Not available ({e})"

    info_rows = [
        ["Python Version", sys.version.split()[0]],
        ["NumPy Version", np.__version__],
        ["GPU (nvidia-smi)", cuda_info],
        ["", ""],
        ["Eventalign CPU", "Yes ✓" if EVENTALIGN_CPU_AVAILABLE else "No ✗"],
        ["Eventalign CUDA", "Yes ✓" if EVENTALIGN_CUDA_AVAILABLE else "No ✗"],
        ["DTW CUDA", "Yes ✓" if DTW_CUDA_AVAILABLE else "No ✗"],
    ]

    print_table(["Component", "Status"], info_rows)

    # Build instructions if something is missing
    if not EVENTALIGN_CPU_AVAILABLE or not EVENTALIGN_CUDA_AVAILABLE or not DTW_CUDA_AVAILABLE:
        print("\n" + "-" * 60)
        print(" BUILD INSTRUCTIONS")
        print("-" * 60)
        if not EVENTALIGN_CPU_AVAILABLE:
            print("\nTo enable CPU eventalign:")
            print("  pip install -e /path/to/pyfin")
        if not EVENTALIGN_CUDA_AVAILABLE or not DTW_CUDA_AVAILABLE:
            print("\nTo enable CUDA extensions (GPU eventalign + DTW):")
            print("  1. Install CUDA Toolkit (nvcc must be in PATH)")
            print("  2. CUDA_HOME=/usr/local/cuda pip install -e /path/to/pyfin")


def main():
    parser = argparse.ArgumentParser(description="Benchmark CPU vs GPU performance")
    parser.add_argument(
        "--signal-sizes",
        type=str,
        default="1000,5000,10000,20000",
        help="Comma-separated signal sizes to test (default: 1000,5000,10000,20000)",
    )
    parser.add_argument(
        "--n-sequences",
        type=str,
        default="5,10,15,20",
        help="Comma-separated number of sequences for DTW pairwise (default: 5,10,15,20)",
    )
    parser.add_argument(
        "--seq-length",
        type=int,
        default=500,
        help="Sequence length for DTW benchmark (default: 500)",
    )
    parser.add_argument(
        "--iterations", type=int, default=3, help="Number of iterations per benchmark (default: 3)"
    )
    parser.add_argument(
        "--skip-cpu-dtw", action="store_true", help="Skip CPU DTW benchmark (can be slow)"
    )

    args = parser.parse_args()

    signal_sizes = [int(x) for x in args.signal_sizes.split(",")]
    n_sequences_list = [int(x) for x in args.n_sequences.split(",")]

    # Print system info
    print_system_info()

    # Benchmark Eventalign CPU
    print("\n")
    eventalign_cpu_results = benchmark_eventalign_cpu(signal_sizes, args.iterations)
    if eventalign_cpu_results:
        print_table(
            ["Signal Size", "Seq Length", "Events", "Time (ms)", "Events/ms"],
            eventalign_cpu_results,
            title="EVENTALIGN BENCHMARK (CPU - Profile HMM)",
        )

    # Benchmark Eventalign CUDA
    print("\n")
    eventalign_cuda_results = benchmark_eventalign_cuda(signal_sizes, args.iterations)
    if eventalign_cuda_results:
        print_table(
            ["Signal Size", "Seq Length", "Events", "Time (ms)", "Events/ms"],
            eventalign_cuda_results,
            title="EVENTALIGN BENCHMARK (CUDA - Profile HMM)",
        )

    # Show eventalign speedup comparison if both available
    if eventalign_cpu_results and eventalign_cuda_results:
        print("\n")
        comparison_rows = []
        for cpu_row, cuda_row in zip(eventalign_cpu_results, eventalign_cuda_results):
            cpu_time = float(cpu_row[3])
            cuda_time = float(cuda_row[3])
            speedup = f"{cpu_time/cuda_time:.1f}x" if cuda_time > 0 else "N/A"
            comparison_rows.append((cpu_row[0], cpu_row[2], cpu_row[3], cuda_row[3], speedup))

        print_table(
            ["Signal Size", "Events", "CPU Time (ms)", "CUDA Time (ms)", "Speedup"],
            comparison_rows,
            title="EVENTALIGN SPEEDUP COMPARISON (CPU vs CUDA)",
        )

    # Benchmark DTW CUDA
    print("\n")
    dtw_cuda_results = benchmark_dtw_pairwise_cuda(
        n_sequences_list, args.seq_length, args.iterations
    )
    if dtw_cuda_results:
        print_table(
            ["N Sequences", "N Pairs", "Time (ms)", "Pairs/ms"],
            dtw_cuda_results,
            title=f"DTW PAIRWISE BENCHMARK (CUDA, seq_len={args.seq_length})",
        )

    # Benchmark DTW CPU (optional - can be slow)
    if not args.skip_cpu_dtw:
        print("\n")
        # Use smaller sizes for CPU to avoid very long runtimes
        cpu_n_sequences = [n for n in n_sequences_list if n <= 15]
        if cpu_n_sequences:
            dtw_cpu_results = benchmark_dtw_pairwise_cpu(
                cpu_n_sequences, args.seq_length, max(1, args.iterations // 2)
            )
            if dtw_cpu_results:
                print_table(
                    ["N Sequences", "N Pairs", "Time (ms)", "Pairs/ms"],
                    dtw_cpu_results,
                    title=f"DTW PAIRWISE BENCHMARK (CPU NumPy, seq_len={args.seq_length})",
                )

        # Show speedup comparison if both available
        if dtw_cuda_results and dtw_cpu_results:
            print("\n")
            print_table(
                ["N Sequences", "CPU Time (ms)", "CUDA Time (ms)", "Speedup"],
                [
                    (
                        cpu_row[0],
                        cpu_row[2],
                        cuda_row[2],
                        (
                            f"{float(cpu_row[2])/float(cuda_row[2]):.1f}x"
                            if float(cuda_row[2]) > 0
                            else "N/A"
                        ),
                    )
                    for cpu_row, cuda_row in zip(
                        dtw_cpu_results, dtw_cuda_results[: len(dtw_cpu_results)]
                    )
                ],
                title="DTW SPEEDUP COMPARISON (CPU vs CUDA)",
            )
    else:
        print("\nNote: CPU DTW benchmark skipped (use --no-skip-cpu-dtw to enable)")

    print("\n" + "=" * 60)
    print(" BENCHMARK COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
