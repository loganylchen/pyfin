# Benchmark Harness

Compare **pyfin** against established long-read isoform quantification tools:
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

Missing tools are **skipped** with a `[SKIP]` message — they do not cause the harness to fail.

---

## Quick Start

```bash
# Make executable once
chmod +x benchmarks/run_benchmark.sh

# Run pyfin only (no external tools needed)
bash benchmarks/run_benchmark.sh \
    --dataset /path/to/data.bam \
    --tools pyfin \
    --output-dir out/

# Run pyfin + bambu (bambu skipped if Rscript absent)
bash benchmarks/run_benchmark.sh \
    --dataset /path/to/data.bam \
    --tools pyfin,bambu \
    --output-dir out/

# Run the full suite
bash benchmarks/run_benchmark.sh \
    --dataset /path/to/data.bam \
    --tools pyfin,bambu,isoquant,jaffal,longgf \
    --output-dir out/

# Show help
bash benchmarks/run_benchmark.sh --help
```

---

## Output Layout

```
out/
├── pyfin/
│   └── result.json          # status marker written by harness
├── bambu/
│   └── result.json
├── isoquant/
│   └── result.json
└── comparison.tsv           # aggregated metrics (auto-generated)
```

Per-tool directories may also contain `*.gtf` outputs; `compare_results.py`
will count transcript features automatically.

---

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

All tests pass even when external tools are absent — missing tools trigger `pytest.skip`.
