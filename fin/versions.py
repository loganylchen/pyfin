"""Versions of PyFIN and of the libraries that decide what a run computes."""
from __future__ import annotations

from importlib import metadata
from typing import Dict, Optional

# Alignment, signal access and signal scoring: a change in any of these can
# change an assembly, so a run records them and --version prints them.
LIBRARIES = ("krill", "mappy", "pysam", "pyslow5", "pod5", "numpy", "scipy")


def library_versions() -> Dict[str, Optional[str]]:
    """Installed version of each library, ``None`` when it is not installed."""
    versions: Dict[str, Optional[str]] = {}
    for name in LIBRARIES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def version_text() -> str:
    from fin import __version__

    libraries = ", ".join(
        f"{name} {version or 'not installed'}" for name, version in library_versions().items()
    )
    return f"PyFIN {__version__} ({libraries})"
