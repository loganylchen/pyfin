"""
fin - A Python package for nanopore Direct RNA-seq data analysis
"""

__version__ = "0.1.0"
__author__ = "loganylchen"

# Import submodules
from . import io
from .utils.log_config import setup_logger, get_package_logger

# Try to import optional extensions
try:
    from ._eventalign import getevents, set_model, MODEL_RNA002, MODEL_RNA004
    from ._eventalign import _EVENTALIGN_AVAILABLE

    _EVENTALIGN_IMPORT_ERROR = None
except ImportError as e:
    _EVENTALIGN_AVAILABLE = False
    getevents = None
    set_model = None
    MODEL_RNA002 = 1
    MODEL_RNA004 = 2
    _EVENTALIGN_IMPORT_ERROR = str(e)

# Try to import _dtw extension (optional CUDA DTW)
try:
    from ._dtw import cuda_dtw, cuda_dtw_batch
    _DTW_AVAILABLE = True
except ImportError:
    _DTW_AVAILABLE = False
    cuda_dtw = None
    cuda_dtw_batch = None


# Initialize package logger
package_logger = get_package_logger(__name__, level="INFO")

# Log package initialization
package_logger.info(f"Package initialized - version {__version__}")

if not _EVENTALIGN_AVAILABLE:
    package_logger.warning(
        f"_eventalign extension not available - event detection disabled. "
        f"Error: {_EVENTALIGN_IMPORT_ERROR}"
    )
if not _DTW_AVAILABLE:
    package_logger.info("_dtw CUDA extension not available (optional)")

__all__ = [
    "io",
    "setup_logger",
    "get_package_logger",
]

if _EVENTALIGN_AVAILABLE:
    __all__.extend(["getevents", "set_model", "MODEL_RNA002", "MODEL_RNA004"])

if _DTW_AVAILABLE:
    __all__.extend(["cuda_dtw", "cuda_dtw_batch"])
