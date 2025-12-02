"""
Logging configuration for the py-fin package

This module provides centralized logging configuration with:
- Time format: YYYY-MM-DD HH:MM:SS
- File name and line number
- Log levels: DEBUG, INFO, WARNING, ERROR, CRITICAL
- Optional file and console handlers
"""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, List


def setup_logger(
    name: str,
    level: str = 'INFO',
    log_file: Optional[str] = None,
    console: bool = True,
    file_mode: str = 'a',
    format_string: Optional[str] = None
) -> logging.Logger:
    """
    Set up a logger with datetime and file name formatting

    Args:
        name: Logger name (use __name__ for module-level logging)
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional path to log file
        console: Whether to log to console
        file_mode: File mode ('w' for overwrite, 'a' for append)
        format_string: Custom format string (optional)

    Returns:
        Configured logger instance

    Example:
        >>> logger = setup_logger(__name__, level='DEBUG', log_file='fin.log')
        >>> logger.info("Process started")
    """
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))

    # Clear existing handlers to avoid duplicates
    logger.handlers.clear()

    # Define format
    if format_string is None:
        # Default format with datetime and filename
        format_string = (
            "%(asctime)s | "
            "%(levelname)-8s | "
            "%(filename)s:%(lineno)d | "
            "%(message)s"
        )

    formatter = logging.Formatter(format_string, datefmt='%Y-%m-%d %H:%M:%S')

    # Console handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # File handler
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode=file_mode)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_package_logger(name: str, **kwargs) -> logging.Logger:
    """
    Get a logger instance configured for the py-fin package

    This is a convenient wrapper around setup_logger with sensible defaults
    for the py-fin package.

    Args:
        name: Logger name (use __name__)
        **kwargs: Additional arguments passed to setup_logger

    Returns:
        Configured logger instance

    Example:
        >>> from fin.utils.log_config import get_package_logger
        >>> logger = get_package_logger(__name__)
        >>> logger.info("Starting analysis")
    """
    # Set default log file if not specified
    if 'log_file' not in kwargs:
        # Create logs directory in current working directory
        log_dir = Path.cwd() / 'logs'
        log_dir.mkdir(exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        kwargs['log_file'] = str(log_dir / f'fin_{timestamp}.log')

    logger = setup_logger(name, **kwargs)

    # Disable propagation for non-root loggers to prevent duplicate logging
    # This ensures messages aren't sent to parent loggers' handlers
    if name != 'fin':  # Allow the package root logger to propagate
        logger.propagate = False

    return logger


def list_log_files(log_dir: str = 'logs') -> List[str]:
    """
    List all log files in the logs directory

    Args:
        log_dir: Directory containing log files

    Returns:
        List of log file paths
    """
    log_path = Path(log_dir)
    if not log_path.exists():
        return []

    return sorted([str(f) for f in log_path.glob('*.log')])
