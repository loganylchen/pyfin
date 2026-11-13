"""pytest configuration shared by tests/unit and tests/integration."""

import sys
from pathlib import Path

# Run against this checkout even when another PyFIN is installed.
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def pytest_collection_modifyitems(config, items):
    """Everything with "integration" in its node id carries the integration
    marker, so `-m "not integration"` selects the fast tests."""
    for item in items:
        if "integration" in item.nodeid:
            item.add_marker("integration")
