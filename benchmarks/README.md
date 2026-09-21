# Benchmark utilities

> **This directory is analysis and plotting utilities, not a reproduction
> harness.** The `run_benchmark.sh` skeleton that used to be advertised here
> wrote `"status": "ran"` markers without invoking any tool; it was removed on
> 2026-09-21. The runs that produce the published comparison live in
> `experiments/prod_validation/` -- see the Benchmarking section of the root
> README.

Intended comparison targets for **pyfin**, long-read isoform quantification tools:
[Bambu](https://github.com/GoekeLab/bambu), [IsoQuant](https://github.com/ablab/IsoQuant),
[JAFFAL](https://github.com/Oshlack/JAFFA), and [LongGF](https://github.com/WGLab/LongGF).

---

## Supported Tools

| Tool      | Required binary | Notes                                      |
|-----------|-----------------|--------------------------------------------|
| pyfin     | `python`        | Always available (this project)            |
| bambu     | `Rscript`       | Requires R + Bioconductor `bambu` package  |
| isoquant  | `isoquant.py`   | Python-based; must be on PATH              |
| jaffal    | `bpipe`         | Requires Bpipe pipeline runner             |
| longgf    | `LongGF`        | Compiled binary must be on PATH            |

## The utilities that remain

| Script | What it does |
|---|---|
| `compare_results.py` | counts transcript features in per-tool GTFs and writes a comparison TSV |
| `evaluate_ablation.py`, `ablation_5rows.py` | ablation-table helpers |
| `viz_*.py`, `replot_class_partitioned.py` | plotting used during development |
| `map_drna.sh` | direct-RNA mapping helper |

These take explicit paths; none of them orchestrates a comparison run. The
runs that produced the published comparison are in
`experiments/prod_validation/` -- see the Benchmarking section of the root
README for which script does what, and for the reproduction gap that is still
open.

## Comparison TSV Schema

`comparison.tsv` is a tab-separated file with a header row:

```
tool	metric	value
pyfin	status	ran
pyfin	dataset	/path/to/data.bam
pyfin	num_transcripts	1234
bambu	status	ran
...
```

| Column  | Description                              |
|---------|------------------------------------------|
| tool    | Tool name (directory name under out/)    |
| metric  | Metric key (from result.json or GTF count) |
| value   | String value                             |

---

## Running compare_results.py Standalone

```bash
python benchmarks/compare_results.py \
    --input-dir out/ \
    --output out/comparison.tsv
```

---

## Environment Setup

### R / Bambu

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("bambu")
```

### IsoQuant

```bash
pip install isoquant  # or follow https://github.com/ablab/IsoQuant
```

### JAFFAL

Follow the [JAFFA install guide](https://github.com/Oshlack/JAFFA/wiki).
Ensure `bpipe` is on PATH.

### LongGF

Download a pre-built binary from the
[LongGF releases page](https://github.com/WGLab/LongGF/releases) and add to PATH.

---

## Running Smoke Tests

```bash
python -m pytest tests/integration/test_benchmark_smoke.py -v
```

These tests exercise `compare_results.py` only and do not depend on any external tool being installed; they skip if that script itself is missing.
