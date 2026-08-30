import pytest

from fin.pipeline.config import (
    PIPELINE_PROFILES,
    PipelineConfig,
    resolve_profile_values,
)


def test_config_minimal_instantiation():
    """Minimal required field is bam_path; all others have defaults."""
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.bam_path == "/tmp/x.bam"
    assert c.gtf_path is None
    assert c.use_gpu is True


def test_config_new_scoring_fields_defaults():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.score_alpha == 0.5
    assert c.prior_weight_cap == 10.0
    assert c.use_prior is True


def test_config_fusion_fields_defaults():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.fusion_enabled is False
    assert c.fusion_min_support == 2
    assert c.fusion_max_dist == 500
    assert c.fusion_flank_bp == 500


def test_config_dtw_subsampling_default():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.max_reads_per_interval_for_dtw == 2000


def test_config_containment_collapse_defaults_off():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.containment_collapse is False
    assert c.containment_3p_tol_bp == 20
    assert c.containment_min_abundance_ratio == 1.0


def test_config_mono_exon_gate_defaults_off():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.drop_mono_exon_novel is False
    assert c.min_mono_exon_reads == 0
    assert c.min_mono_exon_length == 0


def test_config_junction_snap_defaults_off():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.junction_snap is False
    assert c.junction_snap_tolerance == 6
    assert c.junction_snap_min_support == 2
    assert c.junction_snap_min_ratio == 2.0


def test_config_junction_support_gate_default_on():
    # Lever 2 promoted to production default (>=2 reads per novel junction).
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.novel_junction_min_reads == 2
    assert c.novel_junction_reads_tol == 2


def test_config_junction_dominance_gate_defaults_off():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.junction_dominance_filter is False
    assert c.junction_dominance_min_reads == 2
    assert c.junction_dominance_window_bp == 20
    assert c.junction_dominance_tol_bp == 2


def test_config_output_paths_default_none():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.output_gtf is None
    assert c.output_tsv is None
    assert c.output_bedpe is None


def test_config_override_fusion_fields():
    c = PipelineConfig(
        bam_path="/tmp/x.bam",
        fusion_enabled=True,
        fusion_min_support=3,
        fusion_max_dist=1000,
        fusion_flank_bp=250,
    )
    assert c.fusion_enabled is True
    assert c.fusion_min_support == 3
    assert c.fusion_max_dist == 1000
    assert c.fusion_flank_bp == 250


def test_config_override_scoring():
    c = PipelineConfig(bam_path="/tmp/x.bam", score_alpha=0.3, use_prior=False)
    assert c.score_alpha == 0.3
    assert c.use_prior is False


def test_config_backward_compat_no_new_args():
    """Existing instantiation style must still work."""
    c = PipelineConfig(
        bam_path="/tmp/x.bam",
        gtf_path="/tmp/x.gtf",
        genome_fasta_path="/tmp/genome.fa",
        fastq_path="/tmp/reads.fq",
        signal_path="/tmp/signals.slow5",
        signal_format="slow5",
        work_dir="./work",
        use_gpu=True,
        output_gtf="/tmp/out.gtf",
    )
    assert c.output_gtf == "/tmp/out.gtf"


def test_config_isoform_fraction_uses_family_locus_by_default():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.isoform_fraction_locus == "family"


def test_validate_warns_on_inert_m2_em_iteration_knobs(tmp_path, caplog):
    import logging
    bam = tmp_path / "x.bam"
    bam.touch()
    cfg = PipelineConfig(bam_path=str(bam), em_sigma=2.5, em_max_iter=50,
                         em_max_iter_override=3)
    with caplog.at_level(logging.WARNING, logger="fin.pipeline.config"):
        cfg.validate()
    text = caplog.text
    assert "Inert under quant_mode=m2_em" in text
    assert "em_sigma=2.5" in text
    assert "em_max_iter=50" in text
    assert "em_max_iter_override=3" in text


def test_validate_silent_when_m2_em_knobs_default(tmp_path, caplog):
    import logging
    bam = tmp_path / "x.bam"
    bam.touch()
    cfg = PipelineConfig(bam_path=str(bam))
    with caplog.at_level(logging.WARNING, logger="fin.pipeline.config"):
        cfg.validate()
    assert "Inert under quant_mode=m2_em" not in caplog.text


def test_config_post_selection_refit_raw_default_off():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.post_selection_refit is False
    assert c.post_selection_refit_effective is False
    assert c.post_selection_refit_disable_reason is None


def test_config_m2_metric_defaults():
    c = PipelineConfig(bam_path="/tmp/x.bam")
    assert c.m2_metric == "mean"
    assert c.m2_summed_llr_margin == 2.0
    assert c.m2_summed_llr_flank == 6


def test_named_profiles_separate_sirv_and_real_finalize_gates():
    assert PIPELINE_PROFILES["sirv"]["min_abundance"] == 3.0
    assert PIPELINE_PROFILES["sirv"]["strict_novel_abundance_floor"] is False
    assert PIPELINE_PROFILES["sirv"]["floor_gtf_abundance"] is True
    assert PIPELINE_PROFILES["sirv"]["min_fulllen_fraction"] == 0.1
    assert PIPELINE_PROFILES["sirv"]["min_polya5p_reads"] == 0
    assert PIPELINE_PROFILES["sirv"]["post_selection_refit"] is True
    assert PIPELINE_PROFILES["sirv"]["m2_metric"] == "auto"
    assert PIPELINE_PROFILES["sirv"]["m2_summed_llr_margin"] == 1.0
    assert PIPELINE_PROFILES["sirv"]["m2_summed_llr_flank"] == 4
    assert PIPELINE_PROFILES["real-drna"]["min_abundance"] == 1.0
    assert PIPELINE_PROFILES["real-drna"]["strict_novel_abundance_floor"] is True
    assert PIPELINE_PROFILES["real-drna"]["max_soft_mass_ratio"] == 0.0
    assert PIPELINE_PROFILES["real-drna"]["min_fulllen_fraction"] == 0.0
    assert PIPELINE_PROFILES["real-drna"]["min_polya5p_reads"] == 0
    assert PIPELINE_PROFILES["real-drna"]["post_selection_refit"] is True
    assert PIPELINE_PROFILES["real-drna"]["drop_mono_exon_novel"] is True
    assert PIPELINE_PROFILES["real-drna"]["min_mono_exon_reads"] == 5
    assert PIPELINE_PROFILES["real-drna"]["junction_snap"] is True
    assert PIPELINE_PROFILES["real-drna"]["junction_snap_tolerance"] == 6
    assert PIPELINE_PROFILES["real-drna"]["junction_snap_min_support"] == 2
    assert PIPELINE_PROFILES["real-drna"]["junction_snap_min_ratio"] == 2.0
    assert PIPELINE_PROFILES["real-drna"]["m2_metric"] == "summed_llr"
    assert PIPELINE_PROFILES["real-drna"]["m2_summed_llr_margin"] == 1.0
    assert PIPELINE_PROFILES["real-drna"]["m2_summed_llr_flank"] == 8
    assert PIPELINE_PROFILES["real-drna-precision"]["min_novel_reads"] == 2
    assert PIPELINE_PROFILES["real-drna-precision"]["min_mono_exon_reads"] == 5


def test_profile_resolution_preserves_explicit_overrides():
    values = {
        "min_abundance": 99.0,
        "min_fulllen_fraction": 99.0,
        "m2_metric": "off",
    }
    resolved = resolve_profile_values(
        "sirv", values, explicit_fields={"min_abundance", "m2_metric"}
    )
    assert resolved["min_abundance"] == 99.0
    assert resolved["m2_metric"] == "off"
    assert resolved["min_fulllen_fraction"] == 0.1


def test_profile_resolution_rejects_unknown_name():
    with pytest.raises(ValueError, match="unknown profile"):
        resolve_profile_values("unknown", {})


def test_programmatic_from_profile_matches_cli_precedence():
    c = PipelineConfig.from_profile(
        "sirv", bam_path="/tmp/x.bam", min_abundance=2.0, m2_metric="off"
    )
    assert c.profile == "sirv"
    assert c.min_abundance == 2.0
    assert c.floor_gtf_abundance is True
    assert c.min_polya5p_reads == 0
    assert c.m2_metric == "off"
    assert c.profile_overrides == ("m2_metric", "min_abundance")


@pytest.mark.parametrize(
    "field,value",
    [
        ("junction_snap_tolerance", 0),
        ("junction_snap_min_support", 0),
        ("junction_snap_min_ratio", 0.5),
    ],
)
def test_config_validate_rejects_invalid_junction_snap_values(
    tmp_path, field, value
):
    bam = tmp_path / "input.bam"
    bam.touch()
    config = PipelineConfig(bam_path=str(bam))
    setattr(config, field, value)
    with pytest.raises(ValueError, match=field):
        config.validate()


def test_programmatic_sirv_auto_metric_routes_by_guide(tmp_path):
    unguided = PipelineConfig.from_profile("sirv", bam_path="/tmp/x.bam")
    assert unguided.m2_metric == "summed_llr"
    assert unguided.m2_metric_route == "auto-unguided"

    gtf = tmp_path / "guide.gtf"
    gtf.write_text(
        "chr1\tt\ttranscript\t1\t10\t.\t+\t.\ttranscript_id \"a\";\n"
        "chr1\tt\ttranscript\t20\t30\t.\t+\t.\ttranscript_id \"b\";\n"
    )
    guided = PipelineConfig.from_profile(
        "sirv", bam_path="/tmp/x.bam", gtf_path=str(gtf)
    )
    assert guided.m2_metric == "mean"
    assert guided.m2_metric_route == "auto-guided"


def test_config_validate_enables_named_profile_refit(tmp_path):
    bam = tmp_path / "x.bam"
    bam.touch()
    c = PipelineConfig.from_profile("real-drna", bam_path=str(bam))
    c.validate()
    assert c.post_selection_refit_effective is True
    assert c.post_selection_refit_disable_reason is None


def test_config_validate_warn_disables_profile_refit_for_other_mode(
    tmp_path, caplog
):
    bam = tmp_path / "x.bam"
    bam.touch()
    c = PipelineConfig.from_profile(
        "real-drna", bam_path=str(bam), quant_mode="argmax"
    )
    c.validate()
    assert c.post_selection_refit is True
    assert c.post_selection_refit_effective is False
    assert "unsupported quant_mode" in c.post_selection_refit_disable_reason
    assert "refit disabled" in caplog.text.lower()


def test_config_validate_warn_disables_profile_refit_for_feedback(
    tmp_path, caplog
):
    bam = tmp_path / "x.bam"
    bam.touch()
    c = PipelineConfig.from_profile(
        "real-drna", bam_path=str(bam), abundance_feedback=True
    )
    c.validate()
    assert c.post_selection_refit_effective is False
    assert "full EM rerun" in c.post_selection_refit_disable_reason
    assert "refit disabled" in caplog.text.lower()


def test_config_validate_rejects_explicit_refit_with_feedback(tmp_path):
    bam = tmp_path / "x.bam"
    bam.touch()
    c = PipelineConfig(
        bam_path=str(bam),
        post_selection_refit=True,
        abundance_feedback=True,
    )
    with pytest.raises(ValueError, match="abundance-feedback EM"):
        c.validate()


def test_config_validate_rejects_invalid_isoform_fraction_locus(tmp_path):
    bam = tmp_path / "x.bam"
    bam.touch()
    c = PipelineConfig(
        bam_path=str(bam), isoform_fraction_locus="gene"
    )
    with pytest.raises(ValueError, match="isoform_fraction_locus"):
        c.validate()


def test_config_validate_rejects_invalid_m2_metric(tmp_path):
    bam = tmp_path / "input.bam"
    bam.touch()
    c = PipelineConfig(bam_path=str(bam), m2_metric="invalid")
    with pytest.raises(ValueError, match="m2_metric"):
        c.validate()
