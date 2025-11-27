"""
fin: A Python tool for detecting RNA modifications using nanopore Direct RNA-seq data.

This package provides tools to compare signal differences between native RNA
and whole-transcriptomic in vitro transcribed products to identify RNA modifications.
"""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

__author__ = "fin authors"
__email__ = ""
__description__ = "A Python tool for detecting RNA modifications using nanopore Direct RNA-seq data"

__all__ = [
    "core",
    "io",
    "utils",
]
