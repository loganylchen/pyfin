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
    from ._f5c import detect_events, is_available as f5c_is_available

    _F5C_AVAILABLE = True
except ImportError:
    _F5C_AVAILABLE = False
    detect_events = None
    f5c_is_available = lambda: False

try:
    from ._f5c import eventalign
    from ._f5c._eventalign import eventalign as _eventalign_raw
    from ._f5c._eventalign import profile_hmm_eventalign

    _EVENTALIGN_AVAILABLE = True
except ImportError:
    _EVENTALIGN_AVAILABLE = False
    eventalign = None
    _eventalign_raw = None
    profile_hmm_eventalign = None


# Initialize package logger
package_logger = get_package_logger(__name__, level="INFO")

# Log package initialization
package_logger.info(f"Package initialized - version {__version__}")

if not _F5C_AVAILABLE:
    package_logger.warning("f5c extension not available - event detection disabled")
if not _EVENTALIGN_AVAILABLE:
    package_logger.warning("eventalign extension not available")

__all__ = ["io", "setup_logger", "get_package_logger"]
if _F5C_AVAILABLE:
    __all__.extend(["detect_events", "f5c_is_available"])
if _EVENTALIGN_AVAILABLE:
    __all__.extend(["eventalign", "profile_hmm_eventalign"])
