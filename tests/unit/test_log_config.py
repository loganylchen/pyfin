#!/usr/bin/env python3
"""
Unit tests for fin.utils.log_config module.

Tests the logging configuration functions:
- setup_logger: Configure logger with handlers
- get_package_logger: Get package-specific logger
- list_log_files: List log files in directory

Run with:
    pytest tests/unit/test_log_config.py -v
"""

import pytest
import logging
import tempfile
import os
import sys
from pathlib import Path
from datetime import datetime

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.utils.log_config import setup_logger, get_package_logger, list_log_files


class TestSetupLogger:
    """Tests for setup_logger function."""

    def test_setup_logger_returns_logger_instance(self):
        """Test that setup_logger returns a Logger instance."""
        logger = setup_logger("test_logger", console=False)
        assert isinstance(logger, logging.Logger)
        assert logger.name == "test_logger"

    def test_setup_logger_sets_correct_level(self):
        """Test that setup_logger sets the correct logging level."""
        logger = setup_logger("test_debug", level="DEBUG", console=False)
        assert logger.level == logging.DEBUG

        logger = setup_logger("test_info", level="INFO", console=False)
        assert logger.level == logging.INFO

        logger = setup_logger("test_warning", level="WARNING", console=False)
        assert logger.level == logging.WARNING

        logger = setup_logger("test_error", level="ERROR", console=False)
        assert logger.level == logging.ERROR

        logger = setup_logger("test_critical", level="CRITICAL", console=False)
        assert logger.level == logging.CRITICAL

    def test_setup_logger_level_case_insensitive(self):
        """Test that level parameter is case insensitive."""
        logger1 = setup_logger("test_case1", level="debug", console=False)
        logger2 = setup_logger("test_case2", level="DEBUG", console=False)
        logger3 = setup_logger("test_case3", level="Debug", console=False)
        
        assert logger1.level == logger2.level == logger3.level == logging.DEBUG

    def test_setup_logger_console_handler(self):
        """Test that console handler is added when console=True."""
        logger = setup_logger("test_console", console=True)
        
        stream_handlers = [h for h in logger.handlers 
                          if isinstance(h, logging.StreamHandler)]
        assert len(stream_handlers) >= 1

    def test_setup_logger_no_console_handler(self):
        """Test that no console handler is added when console=False."""
        # Clear any existing handlers first
        logger = setup_logger("test_no_console", console=False)
        
        # Check that there are no stream handlers going to stdout
        stdout_handlers = [h for h in logger.handlers 
                         if isinstance(h, logging.StreamHandler) 
                         and h.stream == sys.stdout]
        assert len(stdout_handlers) == 0

    def test_setup_logger_file_handler(self):
        """Test that file handler is added when log_file is specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            logger = setup_logger("test_file", log_file=log_file, console=False)
            
            file_handlers = [h for h in logger.handlers 
                           if isinstance(h, logging.FileHandler)]
            assert len(file_handlers) == 1

    def test_setup_logger_writes_to_file(self):
        """Test that logger actually writes to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            logger = setup_logger("test_write", log_file=log_file, console=False)
            
            test_message = "Test log message"
            logger.info(test_message)
            
            # Flush handlers
            for handler in logger.handlers:
                handler.flush()
            
            # Read the log file
            with open(log_file, 'r') as f:
                content = f.read()
            
            assert test_message in content

    def test_setup_logger_creates_log_directory(self):
        """Test that setup_logger creates directory if it doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "nested", "dir", "test.log")
            logger = setup_logger("test_mkdir", log_file=log_file, console=False)
            
            # Check that directory was created
            assert os.path.exists(os.path.dirname(log_file))
            
            # Check that file can be written
            logger.info("Test message")
            for handler in logger.handlers:
                handler.flush()
            assert os.path.exists(log_file)

    def test_setup_logger_append_mode(self):
        """Test that log file append mode works correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            
            # Write first message
            logger1 = setup_logger("test_append1", log_file=log_file, 
                                  file_mode='a', console=False)
            logger1.info("Message 1")
            for h in logger1.handlers:
                h.flush()
            
            # Write second message
            logger2 = setup_logger("test_append2", log_file=log_file, 
                                  file_mode='a', console=False)
            logger2.info("Message 2")
            for h in logger2.handlers:
                h.flush()
            
            # Both messages should be in file
            with open(log_file, 'r') as f:
                content = f.read()
            
            assert "Message 1" in content
            assert "Message 2" in content

    def test_setup_logger_overwrite_mode(self):
        """Test that log file overwrite mode works correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            
            # Write first message
            logger1 = setup_logger("test_overwrite1", log_file=log_file, 
                                  file_mode='w', console=False)
            logger1.info("Message 1")
            for h in logger1.handlers:
                h.flush()
            
            # Write second message with overwrite
            logger2 = setup_logger("test_overwrite2", log_file=log_file, 
                                  file_mode='w', console=False)
            logger2.info("Message 2")
            for h in logger2.handlers:
                h.flush()
            
            # Only second message should be in file
            with open(log_file, 'r') as f:
                content = f.read()
            
            assert "Message 1" not in content
            assert "Message 2" in content

    def test_setup_logger_custom_format(self):
        """Test that custom format string works."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            custom_format = "CUSTOM | %(levelname)s | %(message)s"
            
            logger = setup_logger("test_format", log_file=log_file,
                                 format_string=custom_format, console=False)
            logger.info("Test")
            for h in logger.handlers:
                h.flush()
            
            with open(log_file, 'r') as f:
                content = f.read()
            
            assert "CUSTOM | INFO | Test" in content

    def test_setup_logger_clears_existing_handlers(self):
        """Test that setup_logger clears existing handlers."""
        logger_name = "test_clear_handlers"
        
        # Setup logger twice
        logger1 = setup_logger(logger_name, console=False)
        initial_handlers = len(logger1.handlers)
        
        logger2 = setup_logger(logger_name, console=False)
        
        # Should have same number of handlers, not doubled
        assert len(logger2.handlers) == initial_handlers


class TestGetPackageLogger:
    """Tests for get_package_logger function."""

    def test_get_package_logger_returns_logger(self):
        """Test that get_package_logger returns a Logger instance."""
        logger = get_package_logger("test_pkg_logger")
        assert isinstance(logger, logging.Logger)

    def test_get_package_logger_creates_log_file(self):
        """Test that get_package_logger creates a log file."""
        # This test may create files in the current directory
        logger = get_package_logger("test_pkg_creates_file")
        
        # Check that a file handler was added
        file_handlers = [h for h in logger.handlers 
                        if isinstance(h, logging.FileHandler)]
        assert len(file_handlers) >= 1

    def test_get_package_logger_accepts_custom_level(self):
        """Test that get_package_logger accepts custom log level."""
        logger = get_package_logger("test_pkg_level", level="DEBUG")
        assert logger.level == logging.DEBUG

    def test_get_package_logger_accepts_custom_log_file(self):
        """Test that get_package_logger accepts custom log file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "custom.log")
            logger = get_package_logger("test_pkg_custom", log_file=log_file)
            
            logger.info("Custom log message")
            for h in logger.handlers:
                h.flush()
            
            assert os.path.exists(log_file)
            with open(log_file, 'r') as f:
                content = f.read()
            assert "Custom log message" in content


class TestListLogFiles:
    """Tests for list_log_files function."""

    def test_list_log_files_empty_directory(self):
        """Test list_log_files with empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = list_log_files(tmpdir)
            assert result == []

    def test_list_log_files_nonexistent_directory(self):
        """Test list_log_files with non-existent directory."""
        result = list_log_files("/nonexistent/path/to/logs")
        assert result == []

    def test_list_log_files_finds_log_files(self):
        """Test list_log_files finds .log files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create some log files
            log_files = ["test1.log", "test2.log", "test3.log"]
            for lf in log_files:
                Path(os.path.join(tmpdir, lf)).touch()
            
            # Also create a non-log file
            Path(os.path.join(tmpdir, "not_a_log.txt")).touch()
            
            result = list_log_files(tmpdir)
            
            assert len(result) == 3
            for lf in log_files:
                assert any(lf in r for r in result)

    def test_list_log_files_returns_sorted_list(self):
        """Test list_log_files returns sorted list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create log files in random order
            for name in ["c.log", "a.log", "b.log"]:
                Path(os.path.join(tmpdir, name)).touch()
            
            result = list_log_files(tmpdir)
            
            # Should be sorted
            names = [os.path.basename(r) for r in result]
            assert names == sorted(names)


class TestLoggerFormat:
    """Tests for log message formatting."""

    def test_log_message_contains_timestamp(self):
        """Test that log messages contain timestamp."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            logger = setup_logger("test_timestamp", log_file=log_file, console=False)
            
            logger.info("Test message")
            for h in logger.handlers:
                h.flush()
            
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Check for date pattern YYYY-MM-DD
            import re
            assert re.search(r'\d{4}-\d{2}-\d{2}', content)

    def test_log_message_contains_level(self):
        """Test that log messages contain level."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            logger = setup_logger("test_level", log_file=log_file, console=False)
            
            logger.info("Info message")
            logger.warning("Warning message")
            for h in logger.handlers:
                h.flush()
            
            with open(log_file, 'r') as f:
                content = f.read()
            
            assert "INFO" in content
            assert "WARNING" in content

    def test_log_message_contains_filename(self):
        """Test that log messages contain filename."""
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = os.path.join(tmpdir, "test.log")
            logger = setup_logger("test_filename", log_file=log_file, console=False)
            
            logger.info("Test message")
            for h in logger.handlers:
                h.flush()
            
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Should contain the test file name
            assert "test_log_config.py" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
