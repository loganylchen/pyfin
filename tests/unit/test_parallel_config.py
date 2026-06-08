"""Tests for the interval-level parallelism config knobs (--threads/--gpu-workers).

Covers PipelineConfig defaults and the validate() invariants:
  threads >= 1; 0 <= gpu_workers <= threads; gpu_workers == 0 when use_gpu is False.
Paths are left empty so validate() skips the file-existence checks and reaches the
parallelism invariants.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.pipeline.config import PipelineConfig  # noqa: E402


def _cfg(**kw):
    # Empty paths -> validate() skips FileNotFoundError checks.
    base = dict(bam_path="", genome_fasta_path="", fastq_path="", signal_path="")
    base.update(kw)
    return PipelineConfig(**base)


class TestParallelDefaults:
    def test_threads_default_one(self):
        assert _cfg().threads == 1

    def test_gpu_workers_default_zero(self):
        assert _cfg().gpu_workers == 0


class TestParallelValidate:
    def test_serial_default_validates(self):
        _cfg().validate()  # no raise

    def test_cpu_parallel_validates(self):
        _cfg(threads=4, gpu_workers=0).validate()

    def test_partition_validates(self):
        _cfg(threads=4, gpu_workers=2, use_gpu=True).validate()

    def test_all_gpu_validates(self):
        _cfg(threads=4, gpu_workers=4, use_gpu=True).validate()

    def test_threads_zero_rejected(self):
        with pytest.raises(ValueError, match="threads must be >= 1"):
            _cfg(threads=0).validate()

    def test_gpu_workers_negative_rejected(self):
        with pytest.raises(ValueError, match="gpu_workers must be in"):
            _cfg(threads=2, gpu_workers=-1).validate()

    def test_gpu_workers_exceeds_threads_rejected(self):
        with pytest.raises(ValueError, match="gpu_workers must be in"):
            _cfg(threads=2, gpu_workers=3, use_gpu=True).validate()

    def test_gpu_workers_without_gpu_rejected(self):
        with pytest.raises(ValueError, match="must be 0 when use_gpu is False"):
            _cfg(threads=4, gpu_workers=2, use_gpu=False).validate()
