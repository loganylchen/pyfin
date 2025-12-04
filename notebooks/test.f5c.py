import numpy as np
from fin import detect_events

# Generate sample raw signal (simulate nanopore data)
raw_signal = np.random.normal(loc=100.0, scale=10.0, size=1000).astype(np.float32)

# Detect events (DNA mode)
events = detect_events(raw_signal, is_rna=False)

# Print results
print(f"Detected {len(events)} events")
print("First event:", events[0])
print(events)
# Output example:
# Detected X events
# First event: {'mean': 98.765, 'stdv': 9.876, 'start': 0, 'length': 42.0}