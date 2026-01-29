"""
pytest configuration for PyFIN eventalign tests.
"""

import sys
from pathlib import Path

# Add project root to Python path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def pytest_configure(config):
    """Configure custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )


def pytest_collection_modifyitems(config, items):
    """Add markers based on test file names."""
    for item in items:
        # Mark integration tests
        if "vs_f5c" in item.nodeid or "integration" in item.nodeid:
            item.add_marker("integration")
        # Mark unit tests
        if "test_event_detection" in item.nodeid or "test_scaling" in item.nodeid:
            item.add_marker("unit")
