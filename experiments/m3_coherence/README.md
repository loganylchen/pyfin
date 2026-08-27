# M3 coherence prototype

This directory preserves the retired read-by-read junction-window DTW
prototype and its focused tests. M3 was removed from the production pipeline,
configuration, and CLI because current named profiles do not use it and its
missing-comparison semantics were not production-ready.

The generic `fin.analysis.assignments.em_with_coherence` engine remains in the
library. Production M1/M2 assignment calls it with `beta=0`; this prototype can
be revisited without keeping an inactive M3 product surface.
