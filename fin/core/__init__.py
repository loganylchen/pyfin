"""
Core module for fin.

This module contains the main logic for RNA modification detection,
signal processing, and statistical analysis.
"""

from . import signal_processor
from . import modification_detector
from . import eventalign
from . import event_detection
from . import f5c_wrapper
from . import completeness
from . import integration_matrix
from . import isoform
from . import isoform_detector
from . import dtw_gpu

__all__ = [
    "signal_processor",
    "modification_detector",
    "eventalign",
    "event_detection",
    "f5c_wrapper",
    "completeness",
    "integration_matrix",
    "isoform",
    "isoform_detector",
    "dtw_gpu",
]
