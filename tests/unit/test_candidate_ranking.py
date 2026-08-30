"""Frozen candidate-ranking model: scoring, filter semantics, exemptions."""
import json
from pathlib import Path

import pytest

from fin.analysis.candidate_evidence import compute_candidate_evidence
from fin.analysis.candidate_ranking import (
    RANKER_V1,
    featurize,
    ranking_filter,
    score_evidence_rows,
)
from fin.analysis.quantification import QuantResult


def _row(**over):
    base = {
        "candidate_id": "x", "abundance": 10.0, "num_reads": 10,
        "soft_hard_ratio": 1.0, "confidence": 1.0, "family_share": 0.9,
        "family_rank": 1, "family_dominant_share": 0.9,
        "is_subchain_of_sibling": 0, "is_superchain_of_sibling": 0,
        "weakest_junction_support": 20.0, "n_junctions_below3": 0,
        "canonical_fraction": 1.0, "end5_support_frac": 0.5,
        "end3_support_frac": 0.8, "fulllen_frac": 0.4, "is_mono": 0,
        "n_exons": 5, "tx_length": 1500,
    }
    base.update(over)
    return base


def test_model_constants_match_committed_provenance_json():
    """Source constants must equal the committed provenance file, bit for bit.

    Unconditional: a missing or drifted provenance file is a failure, never a
    skip - the operating threshold is defined on this exact score scale.
    """
    path = (Path(__file__).parents[2]
            / "experiments/prod_validation/models/candidate_ranker_v1.json")
    m = json.loads(path.read_text())
    assert tuple(m["features"]) == RANKER_V1["features"]
    assert list(RANKER_V1["mean"]) == m["mean"]
    assert list(RANKER_V1["std"]) == m["std"]
    assert list(RANKER_V1["weights"]) == m["weights"]
    assert RANKER_V1["bias"] == m["bias"]
    assert RANKER_V1["score_threshold"] == m["score_threshold"]
    assert "NOT probability-calibrated" in m["kind"]


def test_featurize_handles_strings_and_sentinels():
    numeric = featurize(_row())
    strings = featurize({k: str(v) for k, v in _row().items()})
    assert numeric == strings
    missing = featurize(_row(weakest_junction_support=-1.0,
                              canonical_fraction=-1.0,
                              end5_support_frac=-1.0,
                              end3_support_frac=-1.0,
                              fulllen_frac=-1.0,
                              n_junctions_below3=-1))
    idx = RANKER_V1["features"].index("ends_missing")
    assert missing[idx] == 1.0
    assert missing[RANKER_V1["features"].index("log_weakest_junction")] == 0.0


def test_scores_move_with_evidence_strength():
    strong = score_evidence_rows([_row(candidate_id="s")])["s"]
    weak = score_evidence_rows([
        _row(candidate_id="w", weakest_junction_support=0.0,
             canonical_fraction=0.0, num_reads=2, abundance=2.0,
             family_share=0.05, family_rank=4, is_mono=0)
    ])["w"]
    assert strong > weak


def _qr(cid, source="novel"):
    return QuantResult(candidate_id=cid, abundance=5.0, confidence=1.0,
                       num_assigned_reads=5, source=source, chrom="chr1",
                       strand="+", start=0, end=100, exons=((0, 100),))


def test_ranking_filter_drops_only_low_scoring_novel():
    results = {"lo": _qr("lo"), "hi": _qr("hi"), "ref": _qr("ref", source="gtf")}
    rows = [
        _row(candidate_id="lo", weakest_junction_support=0.0, num_reads=2,
             abundance=2.0, canonical_fraction=0.0, family_share=0.05,
             family_rank=5),
        _row(candidate_id="hi"),
        _row(candidate_id="ref", weakest_junction_support=0.0, num_reads=2,
             abundance=2.0, canonical_fraction=0.0, family_share=0.05,
             family_rank=5),
    ]
    kept, scores, dropped = ranking_filter(results, rows)
    assert dropped == ["lo"]
    assert set(kept) == {"hi", "ref"}       # gtf exempt even at low score
    assert scores["lo"] < RANKER_V1["score_threshold"] <= scores["hi"]


def test_candidate_without_evidence_row_is_kept():
    results = {"a": _qr("a")}
    kept, _scores, dropped = ranking_filter(results, [])
    assert dropped == []
    assert set(kept) == {"a"}


def test_threshold_override():
    results = {"hi": _qr("hi")}
    rows = [_row(candidate_id="hi")]
    _kept, scores, dropped = ranking_filter(results, rows, threshold=1e9)
    assert dropped == ["hi"]


def test_config_validates_ranking_mode(tmp_path):
    from fin.pipeline.config import PipelineConfig
    bam = tmp_path / "x.bam"
    bam.touch()
    ok = PipelineConfig(bam_path=str(bam), ranking_mode="filter")
    ok.validate()
    bad = PipelineConfig(bam_path=str(bam), ranking_mode="bogus")
    with pytest.raises(ValueError):
        bad.validate()


def test_mass_conservation_after_filter_and_refit():
    """Reads of a filtered candidate become orphaned or reassigned, never lost."""
    from fin.analysis.abundance_refit import (
        ResponsibilityLedger,
        refit_survivor_abundance,
    )
    results = {"keep": _qr("keep"), "cut": _qr("cut")}
    rows = [
        _row(candidate_id="keep"),
        _row(candidate_id="cut", weakest_junction_support=0.0, num_reads=2,
             abundance=2.0, canonical_fraction=0.0, family_share=0.05,
             family_rank=5),
    ]
    kept, _s, dropped = ranking_filter(results, rows)
    assert dropped == ["cut"]
    ledger = ResponsibilityLedger(
        weights={
            "r_shared": {"keep": 0.5, "cut": 0.5},
            "r_cut_only": {"cut": 1.0},
        },
        input_read_ids=("r_shared", "r_cut_only"),
    )
    refitted, diag = refit_survivor_abundance(kept, [ledger])
    # shared read renormalizes fully onto the survivor; exclusive read orphans
    assert refitted["keep"].abundance == pytest.approx(1.0)
    assert diag["selection_orphaned_reads"] == 1
    assert diag["mass_balance_error"] < 1e-9


def test_runner_raises_when_explicit_filter_lacks_evidence(monkeypatch, tmp_path):
    """An explicitly requested filter must never silently no-op."""
    from fin.pipeline.runner import PipelineRunner

    runner = PipelineRunner.__new__(PipelineRunner)
    from fin.pipeline.config import PipelineConfig
    bam = tmp_path / "x.bam"
    bam.touch()
    runner.config = PipelineConfig(
        bam_path=str(bam), ranking_mode="filter",
        work_dir=str(tmp_path), output_gtf=None, output_tsv=None,
    )
    monkeypatch.setattr(
        PipelineRunner, "_compute_candidate_evidence_rows",
        lambda self, aggregated, strict=False: None,
    )
    monkeypatch.setattr(
        "fin.pipeline.selection.select_global",
        lambda config, aggregated, fn: (aggregated, []),
    )
    runner._apply_polya5p_filter = lambda *a, **k: None
    runner._gtf_reader = None
    runner._genome_fasta = None
    runner._signal_reader = None
    with pytest.raises(RuntimeError, match="complete candidate evidence"):
        runner._finalize_and_write({"a": _qr("a")}, None, None)


def test_collector_marks_failure_incomplete(tmp_path, monkeypatch):
    from fin.analysis.candidate_evidence import collect_ranking_bam_evidence

    missing = tmp_path / "missing.bam"
    out = collect_ranking_bam_evidence(str(missing))
    assert out.complete is False
    assert out.error


def test_filtered_candidate_with_snap_redirect_mass_flows_to_representative():
    """Ledger weights on a filtered candidate that also had a snap redirect
    must follow the redirect chain of the SURVIVORS, never resurrect."""
    from fin.analysis.abundance_refit import (
        ResponsibilityLedger,
        refit_survivor_abundance,
    )
    results = {"keep": _qr("keep")}
    ledger = ResponsibilityLedger(
        weights={"r1": {"cut": 1.0}, "r2": {"cut": 0.4, "absorbed": 0.6}},
        input_read_ids=("r1", "r2"),
    )
    # 'absorbed' was snap-merged into 'keep'; 'cut' was ranking-filtered.
    refitted, diag = refit_survivor_abundance(
        results, [ledger], snap_redirects={"absorbed": "keep"},
    )
    # r1: only 'cut' -> orphaned. r2: 'absorbed' weight flows to keep,
    # 'cut' weight is lost, then renormalizes to 1.0 on keep.
    assert diag["selection_orphaned_reads"] == 1
    assert refitted["keep"].abundance == pytest.approx(1.0)
    assert diag["mass_balance_error"] < 1e-9


def test_cli_rejects_invalid_ranking_options():
    from click.testing import CliRunner
    from fin.cli import main

    bad = CliRunner().invoke(main, ["--ranking-mode", "bogus"])
    assert bad.exit_code != 0
    assert "Invalid value" in bad.output and "ranking-mode" in bad.output
    bad_thr = CliRunner().invoke(main, ["--ranking-threshold", "notafloat"])
    assert bad_thr.exit_code != 0
    help_out = CliRunner().invoke(main, ["--help"])
    assert "--ranking-mode" in help_out.output
    assert "--endpoint-refine" in help_out.output
    assert "--m2-contrast-stats-jsonl" in help_out.output


@pytest.mark.parametrize(
    "mode,threshold", [("off", None), ("filter", None), ("filter", "-0.75")]
)
def test_cli_valid_ranking_values_reach_config(monkeypatch, tmp_path, mode, threshold):
    """Valid values must actually land in PipelineConfig, not just parse."""
    from click.testing import CliRunner
    import fin.cli as cli_mod
    import fin.pipeline.runner as runner_mod

    seen = {}

    class _StubRunner:
        def __init__(self, config):
            seen["config"] = config

        def run(self):
            return {}

        def cleanup(self):
            pass

    # cli.main imports PipelineRunner inside the function body, so the patch
    # must target the defining module, not the cli namespace.
    monkeypatch.setattr(runner_mod, "PipelineRunner", _StubRunner)
    bam = tmp_path / "x.bam"; bam.touch()
    fq = tmp_path / "r.fq"; fq.touch()
    sig = tmp_path / "s.blow5"; sig.touch()
    gen = tmp_path / "g.fa"; gen.touch()
    args = ["--bam", str(bam), "--genome", str(gen), "--fastq", str(fq),
            "--signal", str(sig), "--output-dir", str(tmp_path / "out"),
            "--ranking-mode", mode]
    if threshold is not None:
        args += ["--ranking-threshold", threshold]
    res = CliRunner().invoke(cli_mod.main, args)
    cfg = seen.get("config")
    assert cfg is not None, f"runner not constructed: {res.output[-800:]}"
    assert cfg.ranking_mode == mode
    assert cfg.ranking_threshold == (None if threshold is None else float(threshold))


def test_training_receipts_frozen_in_source_match_committed_json():
    """Unconditional provenance: source constants == committed model JSON.

    The three training receipts are frozen in the source module so the
    provenance chain holds even when the local evidence artifacts have been
    pruned; no branch may skip this comparison.
    """
    from fin.analysis.candidate_ranking import RANKER_V1_RECEIPTS
    root = Path(__file__).parents[2]
    prov = json.loads(
        (root / "experiments/prod_validation/models/candidate_ranker_v1.json")
        .read_text())["training"]
    for key in ("evidence_sha256", "evidence_run_source_sha256",
                "fit_script_sha256"):
        assert key in RANKER_V1_RECEIPTS
        assert len(RANKER_V1_RECEIPTS[key]) == 64
        assert RANKER_V1_RECEIPTS[key] == prov[key], key
    for key in ("sample", "rows", "grouped_cv_auc"):
        assert key in prov and prov[key] not in (None, "")


def test_fit_script_hash_recomputed_when_artifact_present():
    """When the artifacts exist locally, the frozen receipts must recompute."""
    import hashlib
    from fin.analysis.candidate_ranking import RANKER_V1_RECEIPTS
    root = Path(__file__).parents[2]
    fit = root / "experiments/prod_validation/fit_candidate_ranker.py"
    assert RANKER_V1_RECEIPTS["fit_script_sha256"] == hashlib.sha256(
        fit.read_bytes()).hexdigest()
    ev = (root / "experiments/prod_validation/gencode/_goal_opt/evidence_v1"
          / "prec_tune/candidate_evidence.tsv")
    if ev.exists():
        assert RANKER_V1_RECEIPTS["evidence_sha256"] == hashlib.sha256(
            ev.read_bytes()).hexdigest()
        manifest = ev.parent / "run_manifest.json"
        if manifest.exists():
            assert RANKER_V1_RECEIPTS["evidence_run_source_sha256"] == json.loads(
                manifest.read_text())["source_sha256"]


def test_fusion_candidates_are_exempt_like_gtf():
    results = {"fus": _qr("fus", source="fusion")}
    rows = [_row(candidate_id="fus", weakest_junction_support=0.0, num_reads=2,
                 abundance=2.0, canonical_fraction=0.0, family_share=0.05,
                 family_rank=5)]
    kept, scores, dropped = ranking_filter(results, rows)
    assert dropped == []
    assert scores["fus"] < RANKER_V1["score_threshold"]  # low score, still kept


def test_one_pass_collector_mid_iteration_failure(monkeypatch, tmp_path):
    import pysam
    from fin.analysis.candidate_evidence import collect_ranking_bam_evidence

    class FakeBam:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def fetch(self, **kw):
            read = type("R", (), {
                "is_unmapped": False, "is_secondary": False,
                "is_supplementary": False, "query_name": "r1",
                "reference_start": 10, "reference_end": 60,
                "reference_name": "chr1", "is_reverse": False,
                "cigartuples": ((0, 50),),
            })()
            yield read
            raise OSError("truncated mid-scan")

    monkeypatch.setattr(pysam, "AlignmentFile", FakeBam)
    out = collect_ranking_bam_evidence(str(tmp_path / "x.bam"))
    assert out.complete is False
    assert "OSError" in out.error
    assert out.read_ends == {"r1": (10, 60)}  # partial data retained for audit


def test_fit_script_hash_matches_provenance():
    import hashlib
    root = Path(__file__).parents[2]
    prov = json.loads(
        (root / "experiments/prod_validation/models/candidate_ranker_v1.json")
        .read_text())
    actual = hashlib.sha256(
        (root / "experiments/prod_validation/fit_candidate_ranker.py")
        .read_bytes()).hexdigest()
    assert prov["training"]["fit_script_sha256"] == actual


def test_end_to_end_with_evidence_module():
    """compute_candidate_evidence rows feed the scorer without adaptation."""
    qr = QuantResult(candidate_id="e2e", abundance=8.0, confidence=0.9,
                     num_assigned_reads=8, source="novel", chrom="chr1",
                     strand="+", start=0, end=500,
                     exons=((0, 100), (200, 500)), family_id="F",
                     assigned_read_ids=("r1",))
    rows = compute_candidate_evidence({"e2e": qr})
    scores = score_evidence_rows(rows)
    assert "e2e" in scores and isinstance(scores["e2e"], float)
