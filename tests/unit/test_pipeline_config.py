from fin.pipeline.config import PipelineConfig


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
