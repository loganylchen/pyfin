"""What gets installed: the commands the documentation names, and safe pins."""
from __future__ import annotations

import re
from pathlib import Path

from click.testing import CliRunner

import fin
from fin import versions
from fin.cli import main

ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = (ROOT / "pyproject.toml").read_text()


def _table(name: str) -> str:
    match = re.search(rf"^\[{re.escape(name)}\]\n(.*?)(?=^\[|\Z)", PYPROJECT, re.S | re.M)
    assert match, f"[{name}] missing from pyproject.toml"
    return match.group(1)


def test_both_documented_commands_are_declared():
    scripts = dict(re.findall(r'^(\w+)\s*=\s*"([^"]+)"', _table("project.scripts"), re.M))
    assert scripts == {"pyfin": "fin.cli:main", "fin": "fin.cli:main"}


def test_version_has_one_source():
    project = _table("project")
    assert re.search(r'^dynamic\s*=\s*\["version"\]', project, re.M)
    assert not re.search(r"^version\s*=", project, re.M)
    assert 'version = {attr = "fin.__version__"}' in _table("tool.setuptools.dynamic")
    assert "version=" not in (ROOT / "setup.py").read_text()


def test_pyslow5_floor_excludes_versions_that_drop_every_read():
    for text in (PYPROJECT, (ROOT / "requirements.txt").read_text()):
        assert re.search(r"pyslow5>=1\.", text) and "pyslow5>=0" not in text
    assert "mappy" in (ROOT / "requirements.txt").read_text()


def test_krill_is_never_a_declared_dependency():
    # a different project owns the name on PyPI; krill comes from its own index
    assert not re.search(r'^\s*"krill', PYPROJECT, re.M)
    assert not re.search(r"^krill", (ROOT / "requirements.txt").read_text(), re.M)


def test_version_flag_reports_pyfin_and_its_libraries(monkeypatch):
    monkeypatch.setattr(versions, "library_versions", lambda: {"krill": "9.9", "mappy": None})
    result = CliRunner().invoke(main, ["--version"])
    assert result.exit_code == 0
    assert result.output.strip() == f"PyFIN {fin.__version__} (krill 9.9, mappy not installed)"


def test_library_versions_covers_every_result_bearing_library():
    found = versions.library_versions()
    assert tuple(found) == versions.LIBRARIES
    assert found["pysam"] and found["numpy"]
