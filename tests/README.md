# PyFIN Eventalign Test Suite

## Overview

This test suite validates PyFIN's eventalign implementation against f5c's reference output. The goal is to achieve result parity between the two implementations.

## Test Structure

```
tests/
├── conftest.py                    # pytest configuration
├── run_tests.py                   # Test runner script
├── generate_f5c_reference.sh      # Generate f5c reference outputs
├── verify_test_data.py            # Quick verification of test data
├── test_eventalign_vs_f5c.py     # Main integration tests (f5c comparison)
├── test_event_detection.py        # Unit tests for event detection
├── test_scaling.py               # Unit tests for scaling estimation
├── README.md                      # This file
└── testdata/                      # Test data directory
    ├── RNA004.test.pod5           # Raw signal data (POD5 format)
    ├── RNA004.test.fq.gz          # Basecalled sequences (FASTQ)
    ├── RNA004.test.bam            # Alignments to reference
    ├── RNA004.test.tsv.gz         # f5c eventalign reference output
    ├── RNA004.test.blow5          # Raw signal data (BLOW5 format)
    └── test.fa                    # Reference sequences (69 SIRV transcripts)
```

## Test Data

Test data is located in `tests/testdata/`:

### Current Files

| File | Description | Size |
|------|-------------|------|
| `RNA004.test.pod5` | Raw signal data (POD5 format) | 116.5 MB |
| `RNA004.test.fq.gz` | Basecalled sequences | 2.3 MB |
| `RNA004.test.bam` | Alignments to SIRV reference | 2.8 MB |
| `RNA004.test.tsv.gz` | **f5c eventalign output** (ground truth) | 353.6 MB |
| `RNA004.test.blow5` | Raw signal data (BLOW5 format) | 117.8 MB |
| `test.fa` | Reference sequences (69 SIRV transcripts) | 0.1 MB |

### f5c TSV Format

The `RNA004.test.tsv.gz` file contains f5c eventalign output with the following columns:

```
contig, position, reference_kmer, read_name, strand, event_index,
event_level_mean, event_stdv, event_length, model_kmer, model_mean,
model_stdv, standardized_level, start_idx, end_idx, samples
```

### Verifying Test Data

Run the verification script to check test data integrity:

```bash
python tests/verify_test_data.py
```

## Running Tests

### Prerequisites

```bash
# Install test dependencies
pip install pytest pytest-html

# Ensure PyFIN is installed
pip install -e .
```

### Run All Tests

```bash
# From project root
pytest tests/ -v -s

# Or use the test runner
python tests/run_tests.py
```

### Run Specific Test Categories

```bash
# Unit tests only
pytest tests/ -m unit -v -s
python tests/run_tests.py --unit

# Integration tests only
pytest tests/ -m integration -v -s
python tests/run_tests.py --integration
```

### Generate Reports

```bash
# Generate HTML and XML reports
python tests/run_tests.py --report

# Reports saved to test_results/
```

## Test Categories

### Unit Tests (`-m unit`)

1. **Event Detection Tests** (`test_event_detection.py`)
   - Validates event detection produces reasonable events
   - Checks event coverage of signal
   - Verifies RNA-specific parameters

2. **Scaling Tests** (`test_scaling.py`)
   - Validates Model of the Mean (MoM) scaling estimation
   - Compares scaling parameters with f5c estimates
   - Tests pore model loading

### Integration Tests (`-m integration`)

1. **f5c Comparison Tests** (`test_eventalign_vs_f5c.py`)
   - Full pipeline comparison with f5c reference output
   - Position-level accuracy metrics
   - Event-level correlation and RMSE
   - Coverage analysis

## Test Metrics

| Metric | Current Baseline | Target |
|--------|-----------------|--------|
| Position Match Rate | ~50% | >95% |
| Event Index Correlation | ~0.8 | >0.99 |
| Event Mean RMSE | ~5 pA | <0.5 pA |
| Coverage Overlap | ~50% | >95% |

## Interpreting Results

### Position Match Rate
Percentage of reference positions where both f5c and PyFIN have event alignments.
- Low value indicates missing Profile HMM implementation or alignment QC issues.

### Event Index Correlation
Pearson correlation between event indices at matched positions.
- High correlation (>0.9) indicates events are generally in the same order.
- Low correlation indicates fundamental alignment differences.

### Event Mean RMSE
Root mean square error between event mean values at matched positions.
- High RMSE indicates scaling parameter differences.
- Should be <1 pA for proper parity.

### Coverage Overlap (Jaccard Index)
Intersection over union of covered reference positions.
- Measures how well the two implementations agree on which positions to align.

## Debugging Test Failures

1. **Event detection issues**: Check `test_event_detection.py` output
2. **Scaling issues**: Check `test_scaling.py` comparison with f5c
3. **Alignment issues**: Run `python tests/test_eventalign_vs_f5c.py --report` for detailed position comparison

## Adding New Tests

1. Add test functions to appropriate file
2. Use pytest fixtures from `conftest.py`
3. Mark tests with `@pytest.mark.unit` or `@pytest.mark.integration`
4. Update this README with new test descriptions
