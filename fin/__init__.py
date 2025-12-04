"""
fin - A Python package for nanopore Direct RNA-seq data analysis
"""

__version__ = "0.1.0"
__author__ = "loganylchen"

# Import submodules
from . import io
from .utils.log_config import setup_logger, get_package_logger

# Initialize package logger
package_logger = get_package_logger(__name__, level='INFO')

# Log package initialization
package_logger.info(f"Package initialized - version {__version__}")

__all__ = ["io", "setup_logger", "get_package_logger"]
