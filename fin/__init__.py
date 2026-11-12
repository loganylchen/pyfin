"""
fin - A Python package for nanopore Direct RNA-seq data analysis
"""

__version__ = "0.1.0"
__author__ = "loganylchen"

# Import submodules
from . import io
from .utils.log_config import setup_logger, get_package_logger

# Try to import _dtw extension (optional CUDA DTW)
try:
    from ._dtw import dtw_distance, dtw_pairwise, cleanup, is_available as dtw_is_available

    _DTW_IMPORT_ERROR = None
    _DTW_AVAILABLE = True
except ImportError as e:
    _DTW_AVAILABLE = False
    dtw_distance = None
    dtw_pairwise = None
    cleanup = None
    def dtw_is_available() -> bool:
        return False
    _DTW_IMPORT_ERROR = str(e)

# Try CuPy for GPU-accelerated EM
try:
    import cupy  # noqa: F401  # availability probe
    _CUPY_AVAILABLE = True
except ImportError:
    _CUPY_AVAILABLE = False

# Importing the package configures nothing and writes nothing: the command line
# sets up logging (fin.cli), and a program embedding PyFIN sets up its own.
import logging as _logging

_logging.getLogger(__name__).addHandler(_logging.NullHandler())

__all__ = [
    "io",
    "setup_logger",
    "get_package_logger",
]

if _DTW_AVAILABLE:
    __all__.extend(["dtw_distance", "dtw_pairwise", "cleanup", "dtw_is_available"])
