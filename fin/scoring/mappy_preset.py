"""Resolve the mappy preset used for M1 (read↔candidate sequence alignment).

M1 alignment historically hard-coded ``preset="map-ont"``. To allow ablation
of the alignment preset (e.g. comparing ``map-ont`` against ``sr``) without
touching call sites, the preset is read once from the ``MAPPY_PRESET``
environment variable, defaulting to ``map-ont`` so existing behaviour is
unchanged.

Only M1 sequence-alignment call sites use this. Signal-side tiebreak alignment
(krill) deliberately keeps ``map-ont`` so that the preset under test is the
single isolated variable.
"""
from __future__ import annotations

import os


def get_m1_preset() -> str:
    """Return the mappy preset for M1 alignment (``MAPPY_PRESET`` or map-ont)."""
    return os.environ.get("MAPPY_PRESET", "map-ont") or "map-ont"
