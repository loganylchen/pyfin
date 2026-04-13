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
    dtw_is_available = lambda: False
    _DTW_IMPORT_ERROR = str(e)

# Try CuPy for GPU-accelerated EM
try:
    import cupy
    _CUPY_AVAILABLE = True
except ImportError:
    _CUPY_AVAILABLE = False

# Initialize package logger
package_logger = get_package_logger(__name__, level="INFO")

# Log package initialization
package_logger.info(f"Package initialized - version {__version__}")

if not _DTW_AVAILABLE:
    package_logger.info(
        f"_dtw CUDA extension not available (optional). Error: {_DTW_IMPORT_ERROR}"
    )

if _CUPY_AVAILABLE:
    package_logger.info("CuPy GPU acceleration available for EM algorithm")

__all__ = [
    "io",
    "setup_logger",
    "get_package_logger",
]

if _DTW_AVAILABLE:
    __all__.extend(["dtw_distance", "dtw_pairwise", "cleanup", "dtw_is_available"])
