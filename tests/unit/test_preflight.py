"""Input contracts are checked before any work, and each failure names its fix."""
from __future__ import annotations

import sys
from pathlib import Path

import pysam
import pytest

from fin import preflight
from fin.preflight import PreflightError, PreflightReport

READ_SEQ = "ACGT" * 10


def _bam(tmp_path, *, contigs=(("chr1", 1000),), n_reads=60, spliced=True, with_seq=True,
         index=True, name="reads.bam"):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": n} for c, n in contigs]}
    path = tmp_path / name
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for k in range(n_reads):
            rec = pysam.AlignedSegment(out.header)
            rec.query_name = f"read{k:03d}"
            rec.reference_id = 0
            rec.reference_start = 10 + k
            rec.mapping_quality = 60
            rec.cigartuples = [(0, 20), (3, 100), (0, 20)] if spliced else [(0, 40)]
            rec.query_sequence = READ_SEQ if with_seq else None
            out.write(rec)
    if index:
        pysam.index(str(path))
    return str(path)


def _genome(tmp_path, contigs=(("chr1", 1000),), name="genome.fa"):
    path = tmp_path / name
    with path.open("w") as fh:
        for contig, length in contigs:
            fh.write(f">{contig}\n" + "A" * length + "\n")
    return str(path)


def _gtf(tmp_path, contig="chr1", feature="exon"):
    path = tmp_path / "anno.gtf"
    path.write_text(
        f'{contig}\tsrc\ttranscript\t11\t170\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        f'{contig}\tsrc\t{feature}\t11\t30\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    return str(path)


RNA002 = {"sequencing_kit": "sqk-rna002", "sample_rate": "3012", "experiment_type": "rna"}
RNA004 = {"sequencing_kit": "sqk-rna004", "sample_rate": 4000, "experiment_type": "rna"}


DNA = {"sequencing_kit": "sqk-lsk114", "sample_rate": 5000, "experiment_type": "genomic_dna"}


_POLYA_REFUSAL = preflight._polya_refusal   # the real probe, for the tests about it


@pytest.fixture(autouse=True)
def _backend_has_every_polya_model(monkeypatch):
    """What the installed krill supports is not what these tests are about."""
    monkeypatch.setattr(preflight, "_polya_refusal", lambda pore: None)


class _FakeSignal:
    """Stands in for an opened signal file: its read ids and one header per read group."""

    def __init__(self, ids, groups):
        self.ids = set(ids)
        self._groups = list(groups)

    def infos(self):
        return [dict(info, read_group=n) for n, info in enumerate(self._groups)]

    def close(self):
        pass


def _signal(tmp_path, monkeypatch, ids, *groups):
    path = tmp_path / "reads.blow5"
    path.touch()
    monkeypatch.setattr(preflight, "_open_signal",
                        lambda _path, _fmt: _FakeSignal(ids, groups or (RNA002,)))
    return str(path)


def _all_ids(n=60):
    return [f"read{k:03d}" for k in range(n)]


def _run(tmp_path, monkeypatch, **overrides):
    defaults = dict(
        bam_path=lambda: _bam(tmp_path),
        genome_path=lambda: _genome(tmp_path),
        signal_path=lambda: _signal(tmp_path, monkeypatch, _all_ids()),
        signal_format=lambda: "auto",
        m2_metric=lambda: "off",
    )
    # a default is built only when the test did not supply the input itself
    kwargs = {key: make() for key, make in defaults.items() if key not in overrides}
    kwargs.update(overrides)
    return preflight.run_preflight(**kwargs)


# ---------------------------------------------------------------- signal format
@pytest.mark.parametrize("name, expected", [
    ("a.blow5", "slow5"), ("a.slow5", "slow5"), ("a.BLOW5", "slow5"), ("a.pod5", "pod5"),
])
def test_signal_format_is_detected_from_the_name(tmp_path, name, expected):
    assert preflight.resolve_signal_format(str(tmp_path / name), "auto") == expected


def test_a_directory_of_pod5_files_is_pod5(tmp_path):
    (tmp_path / "run").mkdir()
    (tmp_path / "run" / "batch0.pod5").touch()
    assert preflight.resolve_signal_format(str(tmp_path / "run"), "auto") == "pod5"


def test_an_unrecognised_name_needs_an_explicit_format(tmp_path):
    with pytest.raises(PreflightError, match="--signal-format"):
        preflight.resolve_signal_format(str(tmp_path / "signal.dat"), "auto")
    assert preflight.resolve_signal_format(str(tmp_path / "signal.dat"), "pod5") == "pod5"


def test_a_format_contradicting_the_file_name_is_refused(tmp_path):
    with pytest.raises(PreflightError, match="POD5 file.*--signal-format pod5"):
        preflight.resolve_signal_format(str(tmp_path / "reads.pod5"), "slow5")


def test_pod5_is_refused_with_the_conversion_command(tmp_path, monkeypatch):
    pod5 = tmp_path / "reads.pod5"
    pod5.touch()
    for requested in ("auto", "pod5"):
        with pytest.raises(PreflightError, match="SLOW5/BLOW5 only") as err:
            _run(tmp_path, monkeypatch, signal_path=str(pod5), signal_format=requested)
        assert "blue-crab p2s" in str(err.value)


# ---------------------------------------------------------------- backend
def test_missing_krill_is_reported_with_the_install_command(monkeypatch):
    monkeypatch.setitem(sys.modules, "krill", None)  # makes `import krill` raise ImportError
    with pytest.raises(PreflightError) as err:
        preflight.check_backend("summed_llr")
    message = str(err.value)
    assert preflight.KRILL_INSTALL in message and "--m2-metric off" in message
    assert "PyPI" in message


def test_sequence_only_runs_do_not_need_krill(monkeypatch):
    monkeypatch.setitem(sys.modules, "krill", None)
    preflight.check_backend("off")


class _Krill015:
    """A backend with a poly(A) model for RNA002 only, as krill 0.1.5 is."""

    calls: list = []

    @classmethod
    def Aligner(cls, **kwargs):
        cls.calls.append(kwargs)
        if kwargs.get("polya") and kwargs["pore"] != "rna002":
            raise ValueError("polya=True requires pore='rna002'")
        return object()


def test_the_backend_is_asked_for_the_aligner_the_run_builds(monkeypatch):
    monkeypatch.setitem(sys.modules, "krill", _Krill015)
    monkeypatch.setattr(_Krill015, "calls", [])
    assert _POLYA_REFUSAL("rna002") is None
    assert _POLYA_REFUSAL("rna004") == "polya=True requires pore='rna002'"
    assert _Krill015.calls == [
        dict(pore=pore, use_gpu=False, hmm_confidence=True, polya=True)
        for pore in ("rna002", "rna004")
    ]


def test_no_backend_is_not_a_polya_question(monkeypatch):
    monkeypatch.setitem(sys.modules, "krill", None)
    assert _POLYA_REFUSAL("rna004") is None


def test_a_backend_that_fails_to_import_is_a_warning_not_a_crash(tmp_path, monkeypatch):
    # installed but broken: with --m2-metric off nothing has imported it yet, and
    # the run itself would carry on without poly(A) calls
    class _Broken:
        def find_spec(self, name, path=None, target=None):
            if name == "krill":
                raise ValueError("libkrill.so: undefined symbol")
            return None

    monkeypatch.delitem(sys.modules, "krill", raising=False)
    monkeypatch.setattr(sys, "meta_path", [_Broken(), *sys.meta_path])
    assert _POLYA_REFUSAL("rna002") == "libkrill.so: undefined symbol"
    monkeypatch.setattr(preflight, "_polya_refusal", _POLYA_REFUSAL)
    report = _run(tmp_path, monkeypatch)            # m2_metric="off"
    assert len(report.warnings) == 1 and "undefined symbol" in report.warnings[0]


def test_a_backend_without_a_polya_model_for_the_pore_is_a_warning(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "krill", _Krill015)
    monkeypatch.setattr(preflight, "_polya_refusal", _POLYA_REFUSAL)
    report = _run(tmp_path, monkeypatch,
                  signal_path=_signal(tmp_path, monkeypatch, _all_ids(), RNA004))
    assert report.pore == "rna004"
    assert len(report.warnings) == 1
    warning = report.warnings[0]
    assert "poly(A)" in warning and "'rna004'" in warning and "3' endpoint" in warning
    assert "end sharpness" in warning
    assert "polyA+5' filter (off in the default profile) is skipped" in warning
    assert preflight.KRILL_UPGRADE in warning and "--upgrade" in preflight.KRILL_UPGRADE
    # the same backend on the pore it does model: nothing to say
    report = _run(tmp_path, monkeypatch,
                  signal_path=_signal(tmp_path, monkeypatch, _all_ids(), RNA002))
    assert (report.pore, report.warnings) == ("rna002", [])


# ---------------------------------------------------------------- the whole run
def test_matching_inputs_pass_without_warnings(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, gtf_path=_gtf(tmp_path))
    assert report.warnings == []
    assert (report.signal_format, report.sampled_reads, report.spliced_reads) == ("slow5", 60, 60)
    assert report.signal_hits == 60


def test_a_bam_without_an_index_names_the_fix(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="samtools index"):
        _run(tmp_path, monkeypatch, bam_path=_bam(tmp_path, index=False, name="noidx.bam"))


def test_a_bam_without_read_sequences_is_refused(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="no read\\s+sequence"):
        _run(tmp_path, monkeypatch, bam_path=_bam(tmp_path, with_seq=False, name="noseq.bam"))


def test_an_unspliced_bam_is_a_warning_not_an_error(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, bam_path=_bam(tmp_path, spliced=False, name="flat.bam"))
    assert len(report.warnings) == 1 and "splice-aware" in report.warnings[0]
    assert "GENOME alignment" not in report.warnings[0]


def test_a_transcriptome_alignment_is_called_out(tmp_path, monkeypatch):
    contigs = (("ENST00000000001.1", 1000),)
    report = _run(
        tmp_path, monkeypatch,
        bam_path=_bam(tmp_path, contigs=contigs, spliced=False, name="tx.bam"),
        genome_path=_genome(tmp_path, contigs, name="tx.fa"),
    )
    assert "GENOME alignment" in report.warnings[0]


def test_contig_naming_mismatch_with_the_genome_is_refused(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="'chr1' vs '1'") as err:
        _run(tmp_path, monkeypatch, genome_path=_genome(tmp_path, (("1", 1000),), name="g1.fa"))
    assert "chr1" in str(err.value) and "genome contigs: 1" in str(err.value)


def test_a_different_build_with_the_same_names_is_refused(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="different assembly builds"):
        _run(tmp_path, monkeypatch, genome_path=_genome(tmp_path, (("chr1", 999),), name="g2.fa"))


def test_a_genome_without_an_index_gets_one(tmp_path, monkeypatch):
    genome = _genome(tmp_path)
    assert not Path(genome + ".fai").exists()
    _run(tmp_path, monkeypatch, genome_path=genome)
    assert Path(genome + ".fai").exists()


def test_an_annotation_whose_matching_contig_comes_late_is_accepted(tmp_path, monkeypatch):
    gtf = tmp_path / "late.gtf"
    other = "".join(
        f'chr0\tsrc\texon\t{k * 10 + 1}\t{k * 10 + 5}\t.\t+\t.\tgene_id "g{k}"; transcript_id "t{k}";\n'
        for k in range(250_000)
    )
    gtf.write_text("# header\n" + other
                   + 'chr1\tsrc\texon\t11\t30\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    assert _run(tmp_path, monkeypatch, gtf_path=str(gtf)).warnings == []


def test_annotation_naming_mismatch_is_refused(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="share no contig name"):
        _run(tmp_path, monkeypatch, gtf_path=_gtf(tmp_path, contig="1"))


def test_an_annotation_without_exons_warns(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, gtf_path=_gtf(tmp_path, feature="CDS"))
    assert len(report.warnings) == 1 and "'exon'" in report.warnings[0]


def test_signal_from_another_run_is_refused(tmp_path, monkeypatch):
    signal = _signal(tmp_path, monkeypatch, ["other-1", "other-2"])
    with pytest.raises(PreflightError, match="different runs") as err:
        _run(tmp_path, monkeypatch, signal_path=signal)
    assert "read000" in str(err.value) and "other-" in str(err.value)


def test_signal_for_reads_later_in_the_bam_is_partial_coverage_not_a_mismatch(tmp_path, monkeypatch):
    bam = _bam(tmp_path, n_reads=500, name="big.bam")
    # only the very last read of the BAM has signal: there is no prefix that finds it
    signal = _signal(tmp_path, monkeypatch, _all_ids(500)[-1:], RNA004)
    report = _run(tmp_path, monkeypatch, bam_path=bam, signal_path=signal)
    assert report.signal_hits == 0 and report.pore == "rna004"
    assert len(report.warnings) == 1 and "covers only part of the alignment" in report.warnings[0]


def test_partly_missing_signal_is_a_warning(tmp_path, monkeypatch):
    signal = _signal(tmp_path, monkeypatch, _all_ids()[:45])
    report = _run(tmp_path, monkeypatch, signal_path=signal)
    assert report.signal_hits == 45
    assert "only 45 of the first 60" in report.warnings[0]


def test_an_unreadable_signal_file_is_reported(tmp_path, monkeypatch):
    signal = tmp_path / "reads.blow5"
    signal.touch()

    def broken(_path, _fmt):
        raise RuntimeError("not a slow5 file")

    monkeypatch.setattr(preflight, "_open_signal", broken)
    with pytest.raises(PreflightError, match="cannot read .* as SLOW5: not a slow5 file"):
        _run(tmp_path, monkeypatch, signal_path=str(signal))


def test_fusion_without_the_reads_is_a_warning(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, fusion=True)
    assert len(report.warnings) == 1 and "--fastq" in report.warnings[0]
    fastq = tmp_path / "reads.fq"
    fastq.touch()
    assert _run(tmp_path, monkeypatch, fusion=True, fastq_path=str(fastq)).warnings == []
    assert _run(tmp_path, monkeypatch).warnings == []   # assembly alone never needs it


def test_only_the_sample_is_read(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, bam_path=_bam(tmp_path, n_reads=500, name="big.bam"),
                  signal_path=_signal(tmp_path, monkeypatch, _all_ids(500)))
    assert report.sampled_reads == preflight.SAMPLE_READS == report.signal_hits


# ---------------------------------------------------------------- pore model
@pytest.mark.parametrize("info, expected, evidence", [
    (RNA002, "rna002", "kit"),
    (RNA004, "rna004", "kit"),
    ({"sequencing_kit": "SQK-RNA004-XL", "sample_rate": 3012}, "rna004", "kit"),   # the kit wins
    ({"sequencing_kit": "sqk-rna001", "sample_rate": 3012}, "rna002", "RNA002"),
    ({"sequencing_kit": None, "sample_rate": "4000"}, "rna004", "4000 Hz"),
    ({"sequencing_kit": "", "sample_rate": 3012.0}, "rna002", "3012 Hz"),
    ({"sequencing_kit": "mystery", "sample_rate": 5000}, None, "unrecognised"),
    ({}, None, "neither"),
])
def test_pore_is_read_from_the_signal_header(info, expected, evidence):
    pore, why = preflight.detect_pore(info)
    assert pore == expected and evidence in why


def test_auto_follows_the_header(tmp_path, monkeypatch):
    for info, expected in ((RNA002, "rna002"), (RNA004, "rna004")):
        report = _run(tmp_path, monkeypatch,
                      signal_path=_signal(tmp_path, monkeypatch, _all_ids(), info))
        assert (report.pore, report.warnings) == (expected, [])


def test_auto_falls_back_to_rna002_loudly_when_the_header_is_silent(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, signal_path=_signal(tmp_path, monkeypatch, _all_ids(), {}))
    assert report.pore == "rna002"
    assert len(report.warnings) == 1 and "--pore rna004" in report.warnings[0]


def test_an_explicit_pore_is_kept_but_a_contradicting_header_warns(tmp_path, monkeypatch):
    signal = _signal(tmp_path, monkeypatch, _all_ids(), RNA004)
    report = _run(tmp_path, monkeypatch, signal_path=signal, pore="rna002")
    assert report.pore == "rna002"
    assert len(report.warnings) == 1 and "looks like rna004" in report.warnings[0]
    agreed = _run(tmp_path, monkeypatch, signal_path=signal, pore="rna004")
    assert (agreed.pore, agreed.warnings) == ("rna004", [])


def test_an_explicit_pore_needs_no_header(tmp_path, monkeypatch):
    report = _run(tmp_path, monkeypatch, pore="rna004",
                  signal_path=_signal(tmp_path, monkeypatch, _all_ids(), {}))
    assert (report.pore, report.warnings) == ("rna004", [])


def test_a_dna_run_is_refused(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="genomic_dna.*direct-RNA"):
        _run(tmp_path, monkeypatch, signal_path=_signal(tmp_path, monkeypatch, _all_ids(), DNA))


SILENT = {"sequencing_kit": None, "sample_rate": None, "experiment_type": "rna"}


def test_group_zero_does_not_decide_alone(tmp_path, monkeypatch):
    # a first group that is not RNA at all is ignored when another group is
    signal = _signal(tmp_path, monkeypatch, _all_ids(), DNA, RNA004)
    report = _run(tmp_path, monkeypatch, signal_path=signal)
    assert (report.pore, report.warnings) == ("rna004", [])


def test_read_groups_of_two_chemistries_are_never_resolved_by_picking_one(tmp_path, monkeypatch):
    # whatever the sampled reads happen to be: the file does not say which reads are which
    for groups in ((RNA002, RNA004), (RNA004, RNA002), (RNA002, DNA, RNA004)):
        signal = _signal(tmp_path, monkeypatch, _all_ids(), *groups)
        with pytest.raises(PreflightError, match="different chemistries.*--pore rna002 or --pore rna004"):
            _run(tmp_path, monkeypatch, signal_path=signal)
        report = _run(tmp_path, monkeypatch, signal_path=signal, pore="rna004")
        assert report.pore == "rna004"
        assert len(report.warnings) == 1 and "different" in report.warnings[0]


def test_a_group_that_names_nothing_is_reported_even_when_another_names_a_model(tmp_path, monkeypatch):
    signal = _signal(tmp_path, monkeypatch, _all_ids(), RNA002, SILENT)
    report = _run(tmp_path, monkeypatch, signal_path=signal)
    assert report.pore == "rna002"
    assert len(report.warnings) == 1 and "1 read group(s)" in report.warnings[0]
    assert "scored as rna002" in report.warnings[0]
    explicit = _run(tmp_path, monkeypatch, signal_path=signal, pore="rna004")
    assert explicit.pore == "rna004" and len(explicit.warnings) == 2   # silent group + contradiction


def test_a_file_that_is_dna_in_every_group_is_refused(tmp_path, monkeypatch):
    with pytest.raises(PreflightError, match="genomic_dna.*direct-RNA"):
        _run(tmp_path, monkeypatch, signal_path=_signal(tmp_path, monkeypatch, _all_ids(), DNA, DNA))


def test_the_real_fixture_is_one_rna002_read_group():
    fixture = Path(__file__).resolve().parents[1] / "data" / "sirv1_mini" / "reads.blow5"
    signal = preflight._open_signal(str(fixture), "slow5")
    try:
        assert len(signal.ids) == 91
        infos = signal.infos()
        assert len(infos) == 1 and infos[0]["read_group"] == 0
        assert preflight.detect_pore(infos[0]) == ("rna002", "sequencing kit 'sqk-rna002'")
    finally:
        signal.close()


def test_an_index_without_per_contig_counts_is_not_an_empty_bam(tmp_path, monkeypatch):
    # CRAM: the .crai index stores no mapped-read counts
    genome = _genome(tmp_path)
    cram = tmp_path / "reads.cram"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(cram), "wc", header=header, reference_filename=genome) as out:
        for k in range(60):
            rec = pysam.AlignedSegment(out.header)
            rec.query_name = f"read{k:03d}"
            rec.reference_id = 0
            rec.reference_start = 10 + k
            rec.mapping_quality = 60
            rec.cigartuples = [(0, 20), (3, 100), (0, 20)]
            rec.query_sequence = "A" * 40
            out.write(rec)
    pysam.index(str(cram))
    monkeypatch.setenv("REF_PATH", str(tmp_path))
    with pysam.AlignmentFile(str(cram), "r", reference_filename=genome) as check:
        assert not any(s.mapped for s in check.get_index_statistics())
    report = _run(tmp_path, monkeypatch, bam_path=str(cram), genome_path=genome)
    assert (report.sampled_reads, report.spliced_reads, report.warnings) == (60, 60, [])
    with pytest.raises(PreflightError, match="'chr1' vs '1'"):
        _run(tmp_path, monkeypatch, bam_path=str(cram),
             genome_path=_genome(tmp_path, (("1", 1000),), name="g1.fa"))


def test_report_collects_warnings():
    report = PreflightReport(signal_format="slow5")
    report.warn("one")
    assert report.warnings == ["one"]
