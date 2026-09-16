"""The scoring cache must key on what was scored, not on a file existing.

`expression_rescore.score()` reused any `cmp.stats` already present in a cell
directory. The two GENCODE datasets share one work root and differ only in
scope -- `gencode_sub` is chr1/19/21, `gencode` is the whole genome -- so a
genome-wide campaign would land on the same cell paths as the subset runs that
preceded it and report the SUBSET's numbers under genome-wide labels. Nothing
would look wrong: the tables fill in, the Sn/Pr are plausible, and the whole
comparison is silently answering the previous question.

Found by Codex before the campaign ran.
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "experiments" / "prod_validation"))
from expression_rescore import score  # noqa: E402

STATS = """#= Summary for dataset: query.gtf
#-----------------| Sensitivity | Precision  |
    Transcript level:    42.0     |    77.0    |
 Query mRNAs :    1234 in 1000 loci
Matching transcripts:     500
"""


def _cell(tmp_path, query, truth):
    cell = tmp_path / "cells" / "full" / "pyfin"
    cell.mkdir(parents=True)
    (cell / "cmp.stats").write_text(STATS)
    import hashlib
    (cell / "cache_key.json").write_text(json.dumps({
        "query": [str(query.resolve()), query.stat().st_size, query.stat().st_mtime_ns],
        "truth": hashlib.sha256(truth.read_bytes()).hexdigest(),
        "image": __import__("expression_rescore").IMAGE,
    }))
    return cell


def _job(query, truth, cell):
    return ("S", "full", "pyfin", "fam", str(query), str(truth), 7398, cell, False)


@pytest.fixture
def no_docker(monkeypatch):
    """Observe the recompute path without gffcompare. A call means the cache was
    correctly judged stale; the returned row must then carry the recompute's
    outcome, never the numbers that were sitting in the directory."""
    import expression_rescore

    calls = []

    class Result:
        returncode = 0

    monkeypatch.setattr(expression_rescore.subprocess, "run",
                        lambda *a, **k: calls.append(a) or Result())
    return calls


@pytest.fixture
def inputs(tmp_path):
    query = tmp_path / "assembly.gtf"
    query.write_text("chr1\tt\ttranscript\t1\t9\t.\t+\t.\tid\n")
    truth = tmp_path / "truth.gtf"
    truth.write_text("chr1\tt\ttranscript\t1\t9\t.\t+\t.\tid\n")
    return query, truth


def test_an_unchanged_cell_is_reused(inputs, tmp_path):
    query, truth = inputs
    cell = _cell(tmp_path, query, truth)
    row = score(_job(query, truth, cell))
    assert row["tx_sn"] == 42.0 and row["tx_pr"] == 77.0
    assert row["n_query"] == 1234


def test_a_different_truth_must_not_return_the_cached_numbers(inputs, tmp_path, no_docker):
    """The exact failure: same cell path, different evaluation universe. The
    subset truth's 42.0/77.0 must NOT come back labelled genome-wide."""
    query, truth = inputs
    cell = _cell(tmp_path, query, truth)
    other = tmp_path / "truth_genomewide.gtf"
    other.write_text("chr1\tt\ttranscript\t1\t9\t.\t+\t.\tid\n"
                     "chr7\tt\ttranscript\t1\t9\t.\t+\t.\tid\n")
    row = score(_job(query, other, cell))
    assert len(no_docker) == 1, "a changed truth must force a recompute"
    assert "tx_sn" not in row and row["status"] == "no_stats"


def test_a_rewritten_query_also_invalidates(inputs, tmp_path, no_docker):
    query, truth = inputs
    cell = _cell(tmp_path, query, truth)
    query.write_text("chr1\tt\ttranscript\t1\t9\t.\t+\t.\tid\n"
                     "chr1\tt\ttranscript\t20\t99\t.\t+\t.\tid2\n")
    row = score(_job(query, truth, cell))
    assert len(no_docker) == 1
    assert "tx_sn" not in row


def test_a_cell_with_stats_but_no_key_is_treated_as_stale(inputs, tmp_path, no_docker):
    """Cells written before this fix carry no key. They must be recomputed
    rather than trusted -- their provenance is exactly what is in doubt."""
    query, truth = inputs
    cell = _cell(tmp_path, query, truth)
    (cell / "cache_key.json").unlink()
    row = score(_job(query, truth, cell))
    assert len(no_docker) == 1
    assert "tx_sn" not in row


def test_a_missing_query_is_reported_not_cached(tmp_path):
    truth = tmp_path / "truth.gtf"
    truth.write_text("chr1\tt\ttranscript\t1\t9\t.\t+\t.\tid\n")
    row = score(_job(tmp_path / "absent.gtf", truth, tmp_path / "cell"))
    assert row["status"] == "missing_query"


def test_regenerating_an_identical_truth_keeps_the_cache(inputs, tmp_path, no_docker):
    """build_truth() rewrites the truth file on every invocation with mode "w",
    so its mtime changes even when nothing did. Keying on that would make every
    cell stale on every run and push all 390+ competitor GTFs -- some 230 MB --
    back through gffcompare. Found by Codex."""
    query, truth = inputs
    cell = _cell(tmp_path, query, truth)
    original = truth.read_text()
    truth.write_text(original)  # same bytes, new mtime
    assert truth.stat().st_mtime_ns != json.loads((cell / "cache_key.json").read_text()).get("mtime", -1)
    row = score(_job(query, truth, cell))
    assert no_docker == [], "an unchanged truth must not force a recompute"
    assert row["tx_sn"] == 42.0
