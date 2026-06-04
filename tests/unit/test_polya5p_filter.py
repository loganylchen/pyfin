"""Unit tests for the polyA + 5'-proximity candidate-retention filter.

Covers ``fin.analysis.quantification``:
  - ``polya_read_passes`` — per-read polyA gate (qc PASS & length > min).
  - ``compute_polya5p_support`` — count of assigned reads that BOTH have a
    confident polyA tail AND map 5'-flush to the candidate (strand-aware).
  - ``polya5p_drops`` — the production drop set. UNLIKE the other filters this
    gates GTF candidates too; fusion is exempt; empty polya_map no-ops.
And the runner finalize wiring (``_apply_polya5p_filter``) with krill patched.
"""
from __future__ import annotations

import math
import sys
import types

from fin.analysis.quantification import (
    QuantResult,
    compute_polya5p_support,
    polya5p_drops,
    polya_read_passes,
)


def _qr(
    cid,
    *,
    source="novel",
    strand="+",
    start=100,
    end=900,
    assigned=(),
):
    return QuantResult(
        candidate_id=cid,
        abundance=1.0,
        confidence=1.0,
        num_assigned_reads=len(assigned),
        source=source,
        chrom="chr1",
        strand=strand,
        start=start,
        end=end,
        exons=((start, end),),
        assigned_read_ids=tuple(assigned),
    )


# --------------------------------------------------------------------------
# polya_read_passes
# --------------------------------------------------------------------------
def test_pass_when_qc_pass_and_len_above_min():
    assert polya_read_passes("r", {"r": (20.0, "PASS")}, 10.0) is True


def test_fail_when_qc_not_pass():
    assert polya_read_passes("r", {"r": (50.0, "ADAPTER")}, 10.0) is False


def test_fail_when_len_below_or_equal_min():
    assert polya_read_passes("r", {"r": (10.0, "PASS")}, 10.0) is False
    assert polya_read_passes("r", {"r": (9.9, "PASS")}, 10.0) is False


def test_fail_when_len_none_or_nan():
    assert polya_read_passes("r", {"r": (None, "PASS")}, 10.0) is False
    assert polya_read_passes("r", {"r": (math.nan, "PASS")}, 10.0) is False


def test_fail_when_read_absent():
    assert polya_read_passes("r", {}, 10.0) is False


# --------------------------------------------------------------------------
# compute_polya5p_support
# --------------------------------------------------------------------------
def test_support_counts_polya_and_5p_flush_plus_strand():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2", "r3"))
    polya = {"r1": (20.0, "PASS"), "r2": (20.0, "PASS"), "r3": (20.0, "PASS")}
    ends = {
        "r1": (100, 900),   # 5' flush -> supports
        "r2": (110, 500),   # 5' within 25 (3' irrelevant) -> supports
        "r3": (500, 900),   # 5' off by 400 -> no
    }
    assert compute_polya5p_support(qr, polya, ends, 25, 10.0) == 2


def test_support_requires_polya_even_if_5p_flush():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2"))
    polya = {"r1": (20.0, "PASS"), "r2": (5.0, "PASS")}  # r2 too short
    ends = {"r1": (100, 900), "r2": (100, 900)}
    assert compute_polya5p_support(qr, polya, ends, 25, 10.0) == 1


def test_support_minus_strand_uses_end_as_5p():
    # On '-' strand the genomic 5' end is `end`.
    qr = _qr("c", strand="-", start=100, end=900, assigned=("r1", "r2"))
    polya = {"r1": (20.0, "PASS"), "r2": (20.0, "PASS")}
    ends = {
        "r1": (100, 898),   # read 5'=898 vs cand 5'=900 -> within 25
        "r2": (500, 700),   # read 5'=700 vs 900 -> off by 200
    }
    assert compute_polya5p_support(qr, polya, ends, 25, 10.0) == 1


def test_support_window_boundary_inclusive():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2"))
    polya = {"r1": (20.0, "PASS"), "r2": (20.0, "PASS")}
    ends = {"r1": (125, 900), "r2": (126, 900)}  # 25 inclusive, 26 excluded
    assert compute_polya5p_support(qr, polya, ends, 25, 10.0) == 1


def test_support_read_without_span_not_counted():
    qr = _qr("c", assigned=("r1", "r2"))
    polya = {"r1": (20.0, "PASS"), "r2": (20.0, "PASS")}
    ends = {"r1": (100, 900)}  # r2 has no span
    assert compute_polya5p_support(qr, polya, ends, 25, 10.0) == 1


# --------------------------------------------------------------------------
# polya5p_drops
# --------------------------------------------------------------------------
def _drops(results, polya, ends, min_reads):
    return polya5p_drops(
        {qr.candidate_id: qr for qr in results},
        polya,
        ends,
        window=25,
        min_polya_len=10.0,
        min_reads=min_reads,
    )


def test_drop_novel_below_threshold():
    keep = _qr("keep", assigned=("r1",))
    drop = _qr("drop", assigned=("r2",))
    polya = {"r1": (20.0, "PASS"), "r2": (5.0, "PASS")}
    ends = {"r1": (100, 900), "r2": (100, 900)}
    assert _drops([keep, drop], polya, ends, 1) == {"drop"}


def test_gtf_is_gated_too():
    g = _qr("g", source="gtf", assigned=("r1",))
    polya = {"r1": (5.0, "PASS")}  # fails
    ends = {"r1": (100, 900)}
    assert _drops([g], polya, ends, 1) == {"g"}


def test_fusion_exempt():
    f = _qr("f", source="fusion", assigned=("r1",))
    polya = {"r1": (5.0, "PASS")}  # would fail, but fusion exempt
    ends = {"r1": (100, 900)}
    assert _drops([f], polya, ends, 1) == set()


def test_empty_polya_map_noop():
    qr = _qr("c", assigned=("r1",))
    ends = {"r1": (100, 900)}
    assert _drops([qr], {}, ends, 1) == set()


def test_zero_min_reads_disables():
    qr = _qr("c", assigned=("r1",))
    polya = {"r1": (5.0, "PASS")}  # fails support
    ends = {"r1": (100, 900)}
    assert _drops([qr], polya, ends, 0) == set()


def test_higher_threshold_drops_more():
    qr = _qr("c", assigned=("r1", "r2"))
    polya = {"r1": (20.0, "PASS"), "r2": (20.0, "PASS")}
    ends = {"r1": (100, 900), "r2": (110, 900)}  # both support -> support=2
    assert _drops([qr], polya, ends, 2) == set()
    assert _drops([qr], polya, ends, 3) == {"c"}


# --------------------------------------------------------------------------
# runner finalize wiring (krill patched)
# --------------------------------------------------------------------------
def _runner_with(monkeypatch, polya_map, *, min_polya5p_reads=1,
                 enable_score_filter=True):
    from fin.pipeline.config import PipelineConfig
    from fin.pipeline.runner import PipelineRunner

    cfg = PipelineConfig(
        bam_path="x.bam",
        signal_path="x.blow5",
        min_polya5p_reads=min_polya5p_reads,
        enable_score_filter=enable_score_filter,
        # disable the other novel filters so we isolate polya5p
        min_isoform_fraction=0.0,
        min_fulllen_fraction=0.0,
    )
    runner = PipelineRunner(cfg)
    ends = {"r1": (100, 900), "r2": (100, 900)}
    seqs = {"r1": "ACGT", "r2": "ACGT"}
    monkeypatch.setattr(
        runner, "_fetch_read_seqs_and_ends", lambda: (seqs, ends)
    )
    monkeypatch.setattr(
        "fin.scoring.polya.compute_polya",
        lambda *a, **k: polya_map,
    )
    return runner


def test_apply_filter_drops_unsupported(monkeypatch):
    runner = _runner_with(monkeypatch, {"r1": (20.0, "PASS"), "r2": (5.0, "PASS")})
    agg = {
        "keep": _qr("keep", assigned=("r1",)),
        "drop": _qr("drop", assigned=("r2",)),
    }
    out = runner._apply_polya5p_filter(agg)
    assert set(out) == {"keep"}


def test_apply_filter_noop_when_krill_empty(monkeypatch):
    runner = _runner_with(monkeypatch, {})
    agg = {"c": _qr("c", assigned=("r2",))}  # r2 would fail
    out = runner._apply_polya5p_filter(agg)
    assert set(out) == {"c"}


# --------------------------------------------------------------------------
# compute_polya graceful degradation (fake krill injected into sys.modules)
# --------------------------------------------------------------------------
def _fake_krill(*, init_raises_gpu=False, init_raises_all=False, align_raises=False,
                fail_on_call=None, results=None):
    """Build a fake ``krill`` module with controllable failure modes.

    ``fail_on_call`` (0-indexed) raises only on that align call, letting earlier
    batches succeed first (to exercise the partial-failure path).
    """
    mod = types.ModuleType("krill")
    state = {"calls": 0}

    class _Aligner:
        def __init__(self, *, pore, use_gpu, hmm_confidence, polya):
            if init_raises_all or (init_raises_gpu and use_gpu):
                raise RuntimeError("no GPU")
            self.use_gpu = use_gpu

    def _align(signal_path, reads_variants, aligner=None):
        call = state["calls"]
        state["calls"] += 1
        if align_raises or call == fail_on_call:
            raise RuntimeError("bad signal")
        # Echo a PASS result for each requested read by default.
        if results is not None:
            return results
        return {rid: [{"polya_length": 30.0, "polya_qc": "PASS"}]
                for rid in reads_variants}

    mod.Aligner = _Aligner
    mod.align_reads_variants = _align
    return mod


def test_compute_polya_import_failure_returns_empty(monkeypatch):
    from fin.scoring.polya import compute_polya

    # sys.modules[name] = None makes `import name` raise ImportError.
    monkeypatch.setitem(sys.modules, "krill", None)
    assert compute_polya({"r": "ACGT"}, "x.blow5") == {}


def test_compute_polya_gpu_init_falls_back_to_cpu(monkeypatch):
    from fin.scoring.polya import compute_polya

    fake = _fake_krill(
        init_raises_gpu=True,
        results={"r": [{"polya_length": 30.0, "polya_qc": "PASS"}]},
    )
    monkeypatch.setitem(sys.modules, "krill", fake)
    # use_gpu=True must retry on CPU and still score the read.
    out = compute_polya({"r": "ACGT"}, "x.blow5", use_gpu=True)
    assert out == {"r": (30.0, "PASS")}


def test_compute_polya_first_batch_error_returns_empty(monkeypatch):
    from fin.scoring.polya import compute_polya

    fake = _fake_krill(align_raises=True)
    monkeypatch.setitem(sys.modules, "krill", fake)
    # align_reads_variants raises -> filter disabled (empty map, no exception).
    assert compute_polya({"r": "ACGT"}, "x.blow5") == {}


def test_compute_polya_partial_batch_failure_disables_filter(monkeypatch):
    from fin.scoring.polya import compute_polya

    # 2 reads, batch=1 -> two calls; first succeeds, second fails. The already-
    # scored read MUST NOT leak through (else its candidates would be filtered
    # while the failed read's candidates are silently dropped). Expect {}.
    fake = _fake_krill(fail_on_call=1)
    monkeypatch.setitem(sys.modules, "krill", fake)
    out = compute_polya({"r1": "ACGT", "r2": "ACGT"}, "x.blow5", batch=1)
    assert out == {}


def test_compute_polya_cpu_init_failure_returns_empty(monkeypatch):
    from fin.scoring.polya import compute_polya

    fake = _fake_krill(init_raises_all=True)
    monkeypatch.setitem(sys.modules, "krill", fake)
    assert compute_polya({"r": "ACGT"}, "x.blow5", use_gpu=False) == {}


# --------------------------------------------------------------------------
# config default
# --------------------------------------------------------------------------
def test_config_default_on_at_one():
    from fin.pipeline.config import PipelineConfig

    cfg = PipelineConfig(bam_path="x.bam")
    assert cfg.min_polya5p_reads == 1
    assert cfg.polya5p_window_bp == 25
    assert cfg.min_polya_length == 10.0
