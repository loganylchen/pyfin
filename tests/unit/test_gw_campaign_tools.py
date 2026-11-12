"""The campaign's bookkeeping: which cell a marker names, and the report's statistics."""
from __future__ import annotations

import hashlib
import importlib.util
import math
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parents[2] / "experiments" / "prod_validation"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cells = _load("gw_cells")
report = _load("gw_report")


def _load_caller(name, monkeypatch):
    """The scorer and the compactor import their neighbours by name."""
    monkeypatch.syspath_prepend(str(_DIR))
    return _load(name)


def _lock_is_held(camp):
    try:
        with cells.scoring_lock(camp, wait=False):
            return False
    except BlockingIOError:
        return True


def _campaign(tmp_path, part_a, part_b, done_a=(), done_b=()):
    """Two manifests with their marker directories; returns (camp, parts)."""
    parts = []
    for name, done_dir, offset, rows, done in (
        ("manifest.tsv", "_done", 0, part_a, done_a),
        ("manifest_b.tsv", "_done_b", 1000, part_b, done_b),
    ):
        text = "".join(f"{s}\t{c}\n" for s, c in rows)
        (tmp_path / name).write_text(text)
        (tmp_path / done_dir).mkdir(exist_ok=True)
        for n in done:
            (tmp_path / done_dir / f"{n}.ok").touch()
        parts.append((name, done_dir, offset, hashlib.md5(text.encode()).hexdigest()))
    return tmp_path, tuple(parts)


A = [("lib1", "full"), ("lib1", "p00"), ("lib2", "full")]
B = [("lib2", "p00"), ("lib3", "full")]


# ---------------------------------------------------------------- markers -> cells
def test_a_marker_names_the_physical_line_of_its_own_manifest(tmp_path):
    camp, parts = _campaign(tmp_path, A, B, done_a=(1, 3), done_b=(2,))
    assert cells.finished_cells(camp, parts) == [
        (1, "lib1", "full"), (3, "lib2", "full"), (1002, "lib3", "full")]


def test_no_markers_is_no_cells_not_an_error(tmp_path):
    camp, parts = _campaign(tmp_path, A, B)
    assert cells.finished_cells(camp, parts) == []


def test_a_missing_manifest_is_refused_rather_than_skipped(tmp_path):
    # skipping part B would publish a score table holding part A only
    camp, parts = _campaign(tmp_path, A, B, done_a=(1,), done_b=(1,))
    (camp / "manifest_b.tsv").unlink()
    with pytest.raises(cells.ManifestError, match="manifest_b.tsv.*cannot be read"):
        cells.finished_cells(camp, parts)


def test_an_edited_manifest_is_refused(tmp_path):
    # the same markers would now name other cells
    camp, parts = _campaign(tmp_path, A, B, done_a=(1,))
    (camp / "manifest.tsv").write_text("".join(f"{s}\t{c}\n" for s, c in reversed(A)))
    with pytest.raises(cells.ManifestError, match="md5.*pinned"):
        cells.finished_cells(camp, parts)


@pytest.mark.parametrize("marker", ["0", "4", "x", "1.5"])
def test_a_marker_that_is_not_a_line_of_the_manifest_is_refused(tmp_path, marker):
    camp, parts = _campaign(tmp_path, A, B)
    (camp / "_done" / f"{marker}.ok").touch()
    with pytest.raises(cells.ManifestError, match="not a line number"):
        cells.finished_cells(camp, parts)


def test_one_cell_in_both_parts_is_refused(tmp_path):
    camp, parts = _campaign(tmp_path, A, [("lib1", "p00")], done_a=(2,), done_b=(1,))
    with pytest.raises(cells.ManifestError, match=r"\(lib1, p00\) is also"):
        cells.finished_cells(camp, parts)


def test_a_blank_line_is_a_malformed_manifest_not_a_shifted_one(tmp_path):
    # the sbatch scripts index physical lines; dropping the blank would move every later cell
    text = "lib1\tfull\n\nlib2\tfull\n"
    (tmp_path / "manifest.tsv").write_text(text)
    (tmp_path / "_done").mkdir()
    (tmp_path / "_done" / "3.ok").touch()
    parts = (("manifest.tsv", "_done", 0, hashlib.md5(text.encode()).hexdigest()),)
    with pytest.raises(cells.ManifestError, match="manifest.tsv:2"):
        cells.finished_cells(tmp_path, parts)


def test_the_pinned_hashes_are_the_ones_the_sbatch_scripts_pin():
    pinned = {name: md5 for name, _, _, md5 in cells.PARTS}
    assert f"MANIFEST_MD5={pinned['manifest.tsv']}" in (_DIR / "gw_campaign.sbatch").read_text()
    assert f"MANIFEST_MD5={pinned['manifest_b.tsv']}" in (_DIR / "gw_campaign_b.sbatch").read_text()


# ---------------------------------------------------------------- scorer / compactor lock
def test_the_compactor_cannot_take_the_lock_a_scoring_pass_holds(tmp_path):
    with cells.scoring_lock(tmp_path, wait=True):
        with pytest.raises(BlockingIOError):
            with cells.scoring_lock(tmp_path, wait=False):
                pass
    with cells.scoring_lock(tmp_path, wait=False):   # released: now it can
        pass


def _scorer(tmp_path, monkeypatch, score):
    """The scorer on a two-cell campaign in tmp_path, with gffcompare replaced by ``score``."""
    camp, parts = _campaign(tmp_path, A, B, done_a=(1,), done_b=(2,))
    scorer = _load_caller("gw_headtohead", monkeypatch)
    seen = {}

    def finished(root):
        seen["markers_read_under_lock"] = _lock_is_held(root)
        return cells.finished_cells(root, parts)

    monkeypatch.setattr(scorer, "CAMP", camp)
    monkeypatch.setattr(scorer, "BENCH", tmp_path / "bench")      # no competitor files
    monkeypatch.setattr(scorer, "finished_cells", finished)
    monkeypatch.setattr(scorer, "expressed_ids", lambda *a: set())
    monkeypatch.setattr(scorer, "build_truth", lambda master, ids, out: Path(out).touch() or 7)
    monkeypatch.setattr(scorer, "score", score)
    monkeypatch.setattr("sys.argv", ["gw_headtohead.py"])
    return scorer, camp, seen


def _fake_score(job):
    sample, cond, tool, family, _query, _truth, n_truth, cell, _force = job
    Path(cell).mkdir(parents=True, exist_ok=True)
    for name in ("query.gtf", "cmp.annotated.gtf", "cmp.stats"):
        (Path(cell) / name).write_text("x")
    return dict(sample=sample, cond=cond, tool=tool, family=family, status="ok",
                truth_transcripts=n_truth, tx_sn=50.0, tx_pr=50.0, tx_f1=50.0, n_query=1)


def test_a_scoring_pass_reads_the_markers_and_cleans_up_under_the_lock(tmp_path, monkeypatch, capsys):
    scorer, camp, seen = _scorer(tmp_path, monkeypatch, _fake_score)
    table_path = camp / "h2h" / "scores.tsv"
    table_path.parent.mkdir()
    table_path.write_text("the previous pass\n")
    replace, drop = scorer.os.replace, scorer._drop_scratch

    def publishing(src, dst):
        # the new table is complete beside the old one, which is still the old one
        seen["published_from"] = Path(src).name
        seen["new_rows_before_publish"] = len(Path(src).read_text().splitlines())
        seen["old_table_until_publish"] = Path(dst).read_text()
        seen["published_under_lock"] = _lock_is_held(camp)
        replace(src, dst)

    def dropping(h2h):
        seen["scratch_dropped_under_lock"] = _lock_is_held(camp)
        seen["published_before_scratch_dropped"] = "previous" not in table_path.read_text()
        drop(h2h)

    monkeypatch.setattr(scorer.os, "replace", publishing)
    monkeypatch.setattr(scorer, "_drop_scratch", dropping)
    assert scorer.main() == 0
    assert seen == {
        "markers_read_under_lock": True,
        "published_from": "scores.tsv.partial",
        "new_rows_before_publish": 3,
        "old_table_until_publish": "the previous pass\n",
        "published_under_lock": True,
        "scratch_dropped_under_lock": True,
        "published_before_scratch_dropped": True,
    }
    assert not _lock_is_held(camp)
    table = (camp / "h2h" / "scores.tsv").read_text().splitlines()
    assert len(table) == 3 and not (camp / "h2h" / "scores.tsv.partial").exists()
    assert [row.split("\t")[table[0].split("\t").index("idx")] for row in table[1:]] == ["1", "1002"]
    cell = camp / "h2h" / "cells" / "lib1" / "full" / "pyfin"
    assert sorted(p.name for p in cell.iterdir()) == ["cmp.stats"]     # scratch gone, cache kept


def test_a_pass_that_fails_keeps_the_previous_table_and_releases_the_lock(tmp_path, monkeypatch):
    def broken(job):
        raise RuntimeError("gffcompare died")

    scorer, camp, _ = _scorer(tmp_path, monkeypatch, broken)
    (camp / "h2h").mkdir()
    (camp / "h2h" / "scores.tsv").write_text("the previous pass\n")
    with pytest.raises(RuntimeError, match="gffcompare died"):
        scorer.main()
    assert (camp / "h2h" / "scores.tsv").read_text() == "the previous pass\n"
    assert not _lock_is_held(camp)


def test_a_table_that_cannot_be_written_whole_is_not_published(tmp_path, monkeypatch):
    scorer, camp, _ = _scorer(tmp_path, monkeypatch, _fake_score)
    (camp / "h2h").mkdir()
    (camp / "h2h" / "scores.tsv").write_text("the previous pass\n")

    class _Dies(scorer.csv.DictWriter):
        def writerows(self, rows):
            self.writerow(rows[0])
            raise OSError("No space left on device")

    monkeypatch.setattr(scorer.csv, "DictWriter", _Dies)
    with pytest.raises(OSError, match="No space left"):
        scorer.main()
    assert (camp / "h2h" / "scores.tsv").read_text() == "the previous pass\n"
    assert not _lock_is_held(camp)


def test_the_compactor_refuses_while_a_pass_runs_and_works_under_the_lock(tmp_path, monkeypatch):
    compactor = _load_caller("gw_compact", monkeypatch)
    monkeypatch.setattr(compactor, "CAMP", tmp_path)
    monkeypatch.setattr("sys.argv", ["gw_compact.py", "--dry-run"])
    seen = []
    monkeypatch.setattr(compactor, "compact", lambda dry_run: seen.append(_lock_is_held(tmp_path)))
    with cells.scoring_lock(tmp_path, wait=True):
        with pytest.raises(SystemExit, match="locked"):
            compactor.main()
    assert seen == []
    compactor.main()
    assert seen == [True]


def test_the_scoring_job_fails_when_the_scorer_fails():
    # the wrapper must hand the scorer's status to SLURM and delete nothing itself
    lines = (_DIR / "h2h_score.sbatch").read_text().splitlines()
    code = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    run = next(i for i, ln in enumerate(code) if "gw_headtohead.py" in ln)
    assert code[run + 1].strip() == "RC=$?"
    assert any('exit "$RC"' in ln and '[ "$RC" -eq 0 ] ||' in ln for ln in code[run + 2:])
    assert not any("-delete" in ln or "rm " in ln for ln in code)


# ---------------------------------------------------------------- report statistics
def test_tied_values_share_their_average_rank():
    assert report.ranks([10, 20, 20, 30]) == [0.0, 1.5, 1.5, 3.0]
    assert report.ranks([5, 5, 5]) == [1.0, 1.0, 1.0]


def test_spearman_is_undefined_for_a_constant_side_and_exact_with_ties():
    assert math.isnan(report.spearman([1, 1, 1, 1], [1, 2, 3, 4]))
    # ranks x: .5 .5 2, y: 0 1 2 -> cov 1.5, var 1.5 and 2 -> 1.5 / sqrt(3)
    assert report.spearman([1, 1, 2], [1, 2, 3]) == pytest.approx(1.5 / math.sqrt(3))
    assert report.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)


def test_spearman_leaves_out_pairs_with_an_unknown_value():
    nan = float("nan")
    assert report.spearman([1, 2, nan, 4], [1, 2, -50, 4]) == pytest.approx(1.0)
    assert math.isnan(report.spearman([1, nan, nan, 4], [1, 2, 3, 4]))   # two left: too few


def test_the_library_mean_weights_libraries_not_cells():
    pairs = [("A", 1.0)] * 26 + [("B", -10.0)]
    assert sum(d for _, d in pairs) / len(pairs) == pytest.approx(16 / 27)   # the per-cell mean
    assert report.lib_mean(pairs) == pytest.approx(-4.5)
    assert math.isnan(report.lib_mean([]))


def _score_table(tmp_path, rows):
    h2h = tmp_path / "h2h"
    h2h.mkdir()
    cols = ["sample", "cond", "tool", "status", "tx_sn", "tx_pr", "tx_f1"]
    for name, body in (("scores.tsv", rows), ("passthrough.tsv", [])):
        (h2h / name).write_text(
            "\t".join(cols) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in body))


def test_the_overall_table_puts_each_tool_in_its_own_columns(tmp_path, monkeypatch, capsys):
    # asymmetric on purpose: a swapped column or a flipped sign cannot pass
    rows = []
    for lib, n in (("A", 3), ("B", 1)):
        for k in range(n):
            delta = 2.0 if lib == "A" else -6.0
            rows.append((lib, f"c{k}", "pyfin", "ok", 30.0, 60.0, 40.0 + delta))
            rows.append((lib, f"c{k}", "stringtie3", "ok", 20.0, 70.0, 40.0))
    _score_table(tmp_path, rows)
    monkeypatch.setattr(report, "CAMP", tmp_path)
    monkeypatch.setattr(report, "TOOLS", ("pyfin", "stringtie3", "bambu"))
    report.main()       # library sizes are unknown here: sections 5 and 6 say so and carry on
    out = capsys.readouterr().out
    assert "2 个文库的 BAM 大小读不到" in out and "| A | nan | 3 |" in out
    line = next(ln for ln in out.splitlines() if ln.startswith("| stringtie3 | 4 | "))
    # competitor Sn/Pr/F1, pyfin Sn/Pr/F1, per-cell delta (3*2-6)/4, wins, library mean (2-6)/2
    assert line == "| stringtie3 | 4 | 20.0 | 70.0 | 40.0 | 30.0 | 60.0 | 40.0 | +0.00 | 3/4 | -2.00 |"
    assert "| bambu | 0 | — |" in out                       # no shared cell: a row, not a crash
    assert "合计 4 格:pyfin 领先 3/4,平均 ΔF1 +0.00" in out and "文库等权 ΔF1 -2.00" in out


def test_the_bootstrap_interval_keeps_a_library_together():
    # two libraries, each internally constant: every resample is A+A, A+B or B+B
    lo, hi = report.lib_bootstrap([("A", 0.0)] * 5 + [("B", 10.0)] * 5, n=400, seed=3)
    assert (lo, hi) == (0.0, 10.0)
