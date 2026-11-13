# Tests

```bash
pip install -e ".[dev]"
pip install "krill==0.1.6" --no-deps \
    --index-url https://loganylchen.github.io/krill-dist/simple/
pytest
```

`pytest` with no arguments runs `tests/unit` and `tests/integration`
(`pyproject.toml`). The whole suite takes one to two minutes on one core and
needs no GPU and no data beyond what is in this directory.

## Layout

| Path | What it holds |
|---|---|
| `unit/` | One file per module or behaviour. |
| `integration/` | Whole-pipeline runs: the command line end to end, the fusion analysis, determinism across worker counts. |
| `data/sirv1_mini/` | 91 reads of one SIRV locus **with their raw signal** (RNA002). The only fixture that reaches the signal scoring; `scripts/smoke_test.sh` runs on it too. |
| `data/fusion_k562_bcr_abl1/` | K562 reads over BCR and ABL1, for the fusion analysis. No signal. |
| `conftest.py` | Puts this checkout first on `sys.path`; gives every test with `integration` in its node id the `integration` marker. |

Each fixture directory has a README saying where the data came from and how it
was built.

## The signal backend is required

The end-to-end tests on `data/sirv1_mini` score raw signal with krill. Without
krill they would have to skip, and a suite that skips its only signal tests
proves nothing about the signal path, so CI first checks that krill imports
and builds an aligner (`.github/workflows/python-package.yml`). Install it
from the index above; the package called `krill` on PyPI is something else.

## What is skipped, and why

With krill installed and no GPU, these skip:

- Three tests that need CuPy on a CUDA device: the two in
  `integration/test_gpu_cpu_equiv.py` and one case of
  `unit/test_em_abundance_feedback.py`.
- `integration/test_cpu_fallback.py::test_cpu_fallback_end_to_end` and
  `integration/test_discovery_gtf_vs_no_gtf.py`: written against an RNA004
  dataset, `tests/testdata/RNA004.test.*`, that is not in the repository.
- One test of `integration/test_threads_determinism.py` and one case of
  `unit/test_io_fasta.py`: need data files that are not in the repository.

That is seven skips in a clean Python 3.11 environment, and nothing else
should skip. A new skip in CI is a regression to look at, not noise.

## Selecting tests

```bash
pytest tests/unit                      # unit tests only
pytest -m "not integration"            # leaves out everything under integration/
pytest -m "not slow"
pytest tests/integration/test_sirv1_end_to_end.py    # the signal path, ~10 s
```

Markers are `slow`, `gpu` and `integration`; `--strict-markers` is on, so an
unregistered marker is an error.

## What the suite does not cover

- RNA004 signal. The pore is detected from file headers in unit tests, on
  synthetic headers; no test scores RNA004 signal. See
  [`docs/RNA004_STATUS.md`](../docs/RNA004_STATUS.md).
- A CUDA build of the signal backend.
- Whole-genome behaviour: memory on extreme-depth intervals and accuracy are
  checked by benchmark runs, not by this suite.
