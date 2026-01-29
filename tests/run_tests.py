#!/usr/bin/env python3
"""
Run all PyFIN eventalign tests and generate a comprehensive report.

Usage:
    python run_tests.py              # Run all tests
    python run_tests.py --unit       # Run only unit tests
    python run_tests.py --integration # Run only integration tests
    python run_tests.py --report     # Generate detailed HTML report
"""

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).parent.parent
TESTS_DIR = PROJECT_ROOT / "tests"
RESULTS_DIR = PROJECT_ROOT / "test_results"


def run_tests(unit_only: bool = False, integration_only: bool = False, 
              report: bool = False, verbose: bool = True) -> int:
    """Run tests with specified options."""
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Build pytest command
    cmd = [sys.executable, "-m", "pytest"]
    
    # Add test directory
    cmd.append(str(TESTS_DIR))
    
    # Add markers filter
    if unit_only:
        cmd.extend(["-m", "unit"])
    elif integration_only:
        cmd.extend(["-m", "integration"])
    
    # Add verbosity
    if verbose:
        cmd.extend(["-v", "-s"])
    
    # Add report generation
    if report:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        html_report = RESULTS_DIR / f"test_report_{timestamp}.html"
        junit_report = RESULTS_DIR / f"test_results_{timestamp}.xml"
        
        # Check if pytest-html is available
        try:
            import pytest_html
            cmd.extend(["--html", str(html_report), "--self-contained-html"])
        except ImportError:
            print("Note: pytest-html not installed, skipping HTML report")
        
        # Always generate JUnit XML
        cmd.extend(["--junitxml", str(junit_report)])
    
    # Add color output
    cmd.append("--color=yes")
    
    print("=" * 60)
    print("Running PyFIN Eventalign Tests")
    print("=" * 60)
    print(f"\nCommand: {' '.join(cmd)}\n")
    
    # Run tests
    result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
    
    if report:
        print(f"\nTest results saved to: {RESULTS_DIR}")
    
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description="Run PyFIN eventalign tests")
    parser.add_argument("--unit", action="store_true", 
                        help="Run only unit tests")
    parser.add_argument("--integration", action="store_true",
                        help="Run only integration tests")
    parser.add_argument("--report", action="store_true",
                        help="Generate HTML and XML reports")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Less verbose output")
    
    args = parser.parse_args()
    
    returncode = run_tests(
        unit_only=args.unit,
        integration_only=args.integration,
        report=args.report,
        verbose=not args.quiet
    )
    
    sys.exit(returncode)


if __name__ == "__main__":
    main()
