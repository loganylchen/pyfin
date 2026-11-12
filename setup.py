"""Build shim. All package metadata lives in pyproject.toml.

PyFIN is a pure-Python package: signal scoring (eventalign and DTW) is provided
by the separately installed ``krill`` backend, so nothing is compiled here.
"""
from setuptools import setup

setup()
