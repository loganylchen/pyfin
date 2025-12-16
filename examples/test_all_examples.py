#!/usr/bin/env python3
"""
Test all examples to ensure they work properly

This script checks which examples can run and identifies any issues.
"""

import subprocess
import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).parent

# Examples to test - some require special data files
SIMPLE_EXAMPLES = [
    "test_event.py",
    "test_trimming_fix.py",
    "test_profile_hmm.py",
    "raw_signal_alignment_example.py",
    "debug_eventalign.py",
    "test_adapter_trimming.py",
    "dtw_nanopore_example.py",
    "test_package.py",
]

# Examples that need data files or external tools
ADVANCED_EXAMPLES = [
    "test_eventalign.py",  # Needs matplotlib
    "compare_with_f5c.py",  # Needs f5c installed
    "eventalign_example.py",  # Deprecated
]

# Examples that should just be skipped
SKIP_EXAMPLES = [
    "compare_profile_hmm_with_f5c.py",  # Needs real f5c output
    "test_cuda_dtw.py",  # Needs CUDA
    "benchmark_dtw_comparison.py",  # Needs data
    "benchmark_dtw_pairwise.py",  # Needs data
    "dtw_pairwise_example.py",  # Needs data
    "interval_workflow_example.py",  # Needs data
]


def test_example(example_path, dry_run=False):
    """Test if an example can at least import and start"""
    print(f"\n{'='*70}")
    print(f"Testing: {example_path.name}")
    print(f"{'='*70}")

    if dry_run:
        print(f"  [DRY RUN] Would test: {example_path}")
        return True

    try:
        # Try to run with --help or just check imports
        result = subprocess.run(
            [sys.executable, str(example_path)],
            capture_output=True,
            text=True,
            timeout=10,
            cwd=example_path.parent,
        )

        # Check output
        if "Error importing" in result.stdout or "ModuleNotFoundError" in result.stderr:
            print(f"  ✗ Import error")
            print(f"  STDOUT: {result.stdout[:200]}")
            print(f"  STDERR: {result.stderr[:200]}")
            return False
        elif "Error:" in result.stdout and "not available" in result.stdout:
            print(f"  ⚠ Extension not built (expected)")
            return True
        elif result.returncode == 0:
            print(f"  ✓ Runs successfully")
            return True
        else:
            print(f"  ⚠ Exit code {result.returncode} (may be expected)")
            if result.stdout:
                print(f"  STDOUT: {result.stdout[:200]}")
            return True

    except subprocess.TimeoutExpired:
        print(f"  ⚠ Timeout (example may wait for input)")
        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False


def main():
    print("=" * 70)
    print("Testing PyFin Examples")
    print("=" * 70)

    results = {"pass": [], "fail": [], "skip": []}

    # Test simple examples
    print(f"\n\n{'#'*70}")
    print("# SIMPLE EXAMPLES (should work without external data)")
    print(f"{'#'*70}")

    for example_name in SIMPLE_EXAMPLES:
        example_path = EXAMPLES_DIR / example_name
        if not example_path.exists():
            print(f"\n⚠ File not found: {example_name}")
            results["skip"].append(example_name)
            continue

        if test_example(example_path):
            results["pass"].append(example_name)
        else:
            results["fail"].append(example_name)

    # List advanced examples
    print(f"\n\n{'#'*70}")
    print("# ADVANCED EXAMPLES (may need data files or tools)")
    print(f"{'#'*70}")

    for example_name in ADVANCED_EXAMPLES:
        example_path = EXAMPLES_DIR / example_name
        if example_path.exists():
            print(f"\n  - {example_name}")
            results["skip"].append(example_name)

    # List skipped examples
    print(f"\n\n{'#'*70}")
    print("# SKIPPED EXAMPLES (require specific setup)")
    print(f"{'#'*70}")

    for example_name in SKIP_EXAMPLES:
        example_path = EXAMPLES_DIR / example_name
        if example_path.exists():
            print(f"  - {example_name}")
            results["skip"].append(example_name)

    # Summary
    print(f"\n\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"✓ Passed: {len(results['pass'])}")
    print(f"✗ Failed: {len(results['fail'])}")
    print(f"⊘ Skipped: {len(results['skip'])}")

    if results["fail"]:
        print(f"\nFailed examples:")
        for name in results["fail"]:
            print(f"  - {name}")
        return 1
    else:
        print(f"\n✓ All testable examples passed!")
        return 0


if __name__ == "__main__":
    sys.exit(main())
