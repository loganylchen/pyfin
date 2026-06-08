"""Determinism / zero-regression test for interval-level process parallelism.

Running the SIRV assembly serially (--threads 1) and in parallel (--threads 4),
both --no-gpu, must yield the SAME set of assembled transcript structures with the
SAME abundances. Novel transcript ids are random uuid4 (discovery.py), so we compare
on a structural signature (chrom, strand, sorted exon intervals) rather than
byte-identical text.

Guarded: skips when krill (the m2_em default signal backend) or the SIRV testdata
is unavailable. Slow (runs the full pipeline twice).
"""
from __future__ import annotations

import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent.parent
TD = REPO / "testdata"

_HAVE_KRILL = True
try:  # m2_em default needs krill
    import krill  # noqa: F401
except Exception:  # pragma: no cover - environment dependent
    _HAVE_KRILL = False

_HAVE_DATA = all(
    (TD / f).exists()
    for f in (
        "mapped.bam",
        "SIRV.genome.gtf",
        "SIRV.genome.fa",
        "mapped.fq.gz",
        "mapped.blow5",
    )
)

_TX = re.compile(r'transcript_id "([^"]+)"')
_AB = re.compile(r'abundance "([0-9.eE+-]+)"')


def _run(out_dir: Path, threads: int) -> None:
    args = [
        sys.executable, "-m", "fin.cli",
        "--bam", str(TD / "mapped.bam"),
        "--gtf", str(TD / "SIRV.genome.gtf"),
        "--genome", str(TD / "SIRV.genome.fa"),
        "--fastq", str(TD / "mapped.fq.gz"),
        "--signal", str(TD / "mapped.blow5"),
        "--output-dir", str(out_dir),
        "--no-gpu",
        "--threads", str(threads),
    ]
    r = subprocess.run(args, capture_output=True, text=True, cwd=str(REPO))
    assert r.returncode == 0, f"threads={threads} failed:\n{r.stderr[-3000:]}"


def _structures(gtf: Path) -> dict:
    """Map structural signature -> abundance for each assembled transcript.

    Signature = (chrom, strand, sorted tuple of exon (start, end)).
    """
    exons: dict = defaultdict(list)
    meta: dict = {}
    ab: dict = {}
    for ln in gtf.read_text().splitlines():
        if ln.startswith("#") or "\t" not in ln:
            continue
        f = ln.split("\t")
        if len(f) < 9:
            continue
        m = _TX.search(f[8])
        if not m:
            continue
        tid = m.group(1)
        if f[2] == "exon":
            exons[tid].append((int(f[3]), int(f[4])))
            meta[tid] = (f[0], f[6])
        elif f[2] == "transcript":
            am = _AB.search(f[8])
            ab[tid] = float(am.group(1)) if am else 0.0
    out: dict = {}
    for tid, ex in exons.items():
        chrom, strand = meta[tid]
        sig = (chrom, strand, tuple(sorted(ex)))
        out[sig] = ab.get(tid, 0.0)
    return out


@pytest.mark.slow
@pytest.mark.integration
@pytest.mark.skipif(not _HAVE_KRILL, reason="krill not importable")
@pytest.mark.skipif(not _HAVE_DATA, reason="SIRV testdata not present")
def test_threads_parallel_matches_serial(tmp_path: Path) -> None:
    serial = tmp_path / "serial"
    parallel = tmp_path / "parallel"
    serial.mkdir()
    parallel.mkdir()

    _run(serial, threads=1)
    _run(parallel, threads=4)

    s = _structures(serial / "assembly.gtf")
    p = _structures(parallel / "assembly.gtf")

    assert set(s) == set(p), (
        f"structure sets differ: serial-only={len(set(s)-set(p))}, "
        f"parallel-only={len(set(p)-set(s))}"
    )
    for sig in s:
        assert s[sig] == pytest.approx(p[sig], rel=1e-6, abs=1e-9), (
            f"abundance mismatch for {sig}: serial={s[sig]} parallel={p[sig]}"
        )
