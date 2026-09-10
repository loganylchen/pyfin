"""The selection audit (--selection-audit) records candidate fates and changes nothing else."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

from fin.analysis.quantification import QuantResult
from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner
from fin.pipeline.selection import (
    AUDIT_COLUMNS,
    AUDIT_GLOBAL_FILE,
    AUDIT_SHARD_DIR,
    SelectionOutcome,
    append_selection_audit_global,
    audit_rows_for_candidates,
    audit_rows_for_outcomes,
    write_selection_audit_shard,
)

GENOME = "AAAAA" + "GT" + "CCCCCC" + "AG" + "AAAAA"
CANON_INTRON = (5, 15)
NONCANON_INTRON = (6, 15)


def _cand(cid, introns, source, reads=1):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=15, sequence="ACGT", source=source,
        supporting_read_ids={f"r_{cid}_{i}" for i in range(reads)},
        chrom="chr1", strand="+", start=0, end=20,
    )


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=20, strand="+")


def _read_tsv(path):
    with open(path) as fh:
        return [r for r in csv.reader(fh, delimiter="\t")]


class TestRowBuilders:
    def test_candidate_rows_carry_source_and_support(self):
        rows = audit_rows_for_candidates("chr1:0-20:+", "discovered",
                                         [_cand("a", [CANON_INTRON], "gtf", reads=3)])
        assert rows == [("chr1:0-20:+", "discovered", "a", "gtf", 3, "present", "", "")]
        assert len(rows[0]) == len(AUDIT_COLUMNS)

    def test_outcome_rows_carry_action_reason_and_fold_target(self):
        rows = audit_rows_for_outcomes("iv", "interval_cascade", [
            SelectionOutcome("x", "drop", "m2_support_gate"),
            SelectionOutcome("y", "fold", "containment", fold_into="p"),
        ])
        assert rows == [("iv", "interval_cascade", "x", "", "", "drop", "m2_support_gate", ""),
                        ("iv", "interval_cascade", "y", "", "", "fold", "containment", "p")]


class TestWriters:
    def test_shard_is_written_atomically_and_named_by_interval(self, tmp_path):
        write_selection_audit_shard(tmp_path, _interval(), [("iv", "s", "c", "", 0, "present", "", "")])
        d = tmp_path / AUDIT_SHARD_DIR
        assert sorted(p.name for p in d.iterdir()) == ["chr1_0_20_+.tsv"]
        assert _read_tsv(d / "chr1_0_20_+.tsv") == [["iv", "s", "c", "", "0", "present", "", ""]]

    def test_empty_rows_write_nothing(self, tmp_path):
        write_selection_audit_shard(tmp_path, _interval(), [])
        append_selection_audit_global(tmp_path, [])
        assert not (tmp_path / AUDIT_SHARD_DIR).exists()
        assert not (tmp_path / AUDIT_GLOBAL_FILE).exists()

    def test_global_file_gets_one_header_across_appends(self, tmp_path):
        append_selection_audit_global(tmp_path, [("*", "a", "c1", "", "", "drop", "r", "")])
        append_selection_audit_global(tmp_path, [("*", "b", "c2", "", "", "keep", "", "")])
        rows = _read_tsv(tmp_path / AUDIT_GLOBAL_FILE)
        assert rows[0] == list(AUDIT_COLUMNS)
        assert [r[1] for r in rows[1:]] == ["a", "b"]


class TestGateSeamsReturnOutcomes:
    def test_canonical_gate_seam_returns_the_drop(self):
        cs = CandidateSet(interval=_interval(),
                          candidates=[_cand("good", [CANON_INTRON], "novel"),
                                      _cand("bad", [NONCANON_INTRON], "novel")],
                          read_ids={"r"})
        runner = PipelineRunner(PipelineConfig(bam_path="/fake/reads.bam"))
        outcomes = runner._apply_canonical_gate(cs, GENOME)
        assert [(o.candidate_id, o.action) for o in outcomes] == [("bad", "drop")]
        assert [c.candidate_id for c in cs.candidates] == ["good"]


def _qr(cid, source="gtf"):
    return QuantResult(candidate_id=cid, gene_id=cid, chrom="chr1", strand="+", start=0, end=100,
                       source=source, abundance=2.0, confidence=1.0, num_assigned_reads=2,
                       exons=((0, 40), (60, 100)))


def _run(tmp_path, audit: bool):
    out = tmp_path / ("on" if audit else "off")
    out.mkdir()
    cfg = PipelineConfig(bam_path="/fake/reads.bam", work_dir=str(out),
                         output_gtf=str(out / "assembly.gtf"),
                         output_tsv=str(out / "scores.tsv"),
                         selection_audit=audit)
    runner = PipelineRunner(cfg)
    aggregated = {"tx1": _qr("tx1"), "tx2": _qr("tx2", "novel")}
    with patch("fin.pipeline.runner.generate_isolated_intervals",
               return_value={"intervals": []}), \
         patch("fin.pipeline.runner.aggregate_across_intervals",
               return_value=aggregated):
        runner.run()
    return out


class TestRunnerGlobalAudit:
    def test_off_writes_no_audit_files(self, tmp_path):
        out = _run(tmp_path, audit=False)
        assert not (out / AUDIT_GLOBAL_FILE).exists()
        assert not (out / AUDIT_SHARD_DIR).exists()

    def test_on_records_global_stages_and_final_set(self, tmp_path):
        out = _run(tmp_path, audit=True)
        rows = _read_tsv(out / AUDIT_GLOBAL_FILE)
        assert rows[0] == list(AUDIT_COLUMNS)
        stages = [r[1] for r in rows[1:]]
        assert "final" in stages
        final = {r[2]: r for r in rows[1:] if r[1] == "final"}
        assert set(final) == {"tx1", "tx2"}
        assert final["tx1"][3] == "gtf" and final["tx1"][4] == "2" and final["tx1"][5] == "keep"

    def test_audit_leaves_gtf_and_scores_byte_identical(self, tmp_path):
        off = _run(tmp_path, audit=False)
        on = _run(tmp_path, audit=True)
        assert (off / "assembly.gtf").read_bytes() == (on / "assembly.gtf").read_bytes()
        assert (off / "scores.tsv").read_bytes() == (on / "scores.tsv").read_bytes()


class TestProcessIntervalAudit:
    """process_interval writes one shard covering discovery -> gates -> survivors."""

    def test_shard_records_discovered_gate_drop_and_survivor(self, tmp_path):
        cfg = PipelineConfig(bam_path="/fake/reads.bam", work_dir=str(tmp_path),
                             canonical_gate=True, quant_mode="argmax",
                             selection_audit=True)
        runner = PipelineRunner(cfg)
        runner._genome_fasta = {"chr1": GENOME}
        cs = CandidateSet(interval=_interval(),
                          candidates=[_cand("good", [CANON_INTRON], "novel"),
                                      _cand("bad", [NONCANON_INTRON], "novel")],
                          read_ids={"r1"})
        survivor = _qr("good", "novel")
        with patch("fin.pipeline.runner.discover_candidates", return_value=cs), \
             patch.object(PipelineRunner, "_quant_argmax_keep", return_value=[survivor]):
            results, _ledger = runner.process_interval(_interval())
        assert [r.candidate_id for r in results] == ["good"]
        rows = _read_tsv(tmp_path / AUDIT_SHARD_DIR / "chr1_0_20_+.tsv")
        by_stage = {}
        for r in rows:
            by_stage.setdefault(r[1], []).append(r)
        assert sorted(r[2] for r in by_stage["discovered"]) == ["bad", "good"]
        assert [(r[2], r[5]) for r in by_stage["canonical_gate"]] == [("bad", "drop")]
        assert [(r[2], r[5]) for r in by_stage["interval_survivor"]] == [("good", "keep")]

    def test_no_candidates_still_leaves_a_discovery_record(self, tmp_path):
        cfg = PipelineConfig(bam_path="/fake/reads.bam", work_dir=str(tmp_path),
                             canonical_gate=True, selection_audit=True)
        runner = PipelineRunner(cfg)
        runner._genome_fasta = {"chr1": GENOME}
        cs = CandidateSet(interval=_interval(),
                          candidates=[_cand("bad", [NONCANON_INTRON], "novel")],
                          read_ids={"r1"})
        with patch("fin.pipeline.runner.discover_candidates", return_value=cs):
            assert runner.process_interval(_interval()) is None
        rows = _read_tsv(tmp_path / AUDIT_SHARD_DIR / "chr1_0_20_+.tsv")
        assert [(r[1], r[2], r[5]) for r in rows] == [
            ("discovered", "bad", "present"), ("canonical_gate", "bad", "drop")]
