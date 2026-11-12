"""Importing PyFIN, or asking it for its version, must not write anywhere."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def _run(tmp_path, *code_or_args, module=False):
    env_path = str(REPO)
    command = [sys.executable, "-m", "fin.cli", *code_or_args] if module else [
        sys.executable, "-c", *code_or_args]
    return subprocess.run(command, cwd=tmp_path, capture_output=True, text=True,
                          env={"PYTHONPATH": env_path, "PATH": "/usr/bin:/bin"})


def test_importing_the_package_creates_nothing_and_prints_nothing(tmp_path):
    result = _run(tmp_path, "import fin, fin.io.interval_manager, fin.analysis.clustering, "
                            "fin.io.io_read_manager, fin.pipeline.runner; print('ok')")
    assert result.returncode == 0, result.stderr
    assert result.stdout == "ok\n"
    assert list(tmp_path.iterdir()) == []


def test_importing_works_in_a_directory_that_cannot_be_written(tmp_path):
    tmp_path.chmod(0o555)
    try:
        result = _run(tmp_path, "import fin; print(fin.__version__)")
    finally:
        tmp_path.chmod(0o755)
    assert result.returncode == 0, result.stderr


def test_version_and_help_leave_the_working_directory_alone(tmp_path):
    for flag in ("--version", "--help", "--help-all"):
        result = _run(tmp_path, flag, module=True)
        assert result.returncode == 0, result.stderr
    assert list(tmp_path.iterdir()) == []
