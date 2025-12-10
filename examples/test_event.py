"""
Example: Nanopore Event Detection with Visualization

This example demonstrates:
1. Generating synthetic nanopore signal with adapter region
2. Detecting events using f5c's event detection algorithm
3. Visualizing raw signal, detected events, and trimmed regions
"""

import numpy as np
import matplotlib.pyplot as plt
from fin import detect_events
import pickle

def generate_synthetic_signal(n_samples=2000, seed=42):
    """
    Generate synthetic nanopore signal with:
    - Adapter region (first 200-300 samples): high current, high noise
    - Main signal region: lower current with step changes (events)
    - Tail region (last 50 samples): noisy
    """
    np.random.seed(seed)
    
    # Adapter region (200-300 samples): high unstable current
    adapter_length = 750
    adapter = np.random.normal(loc=150.0, scale=20.0, size=adapter_length)
    # Add some drift to adapter
    adapter += np.linspace(0, -30, adapter_length)
    
    # Main signal region: simulate events with different current levels
    main_length = n_samples - adapter_length - 50
    main_signal = []
    
    # Create events with varying durations and levels
    current_pos = 0
    current_levels = [100.0, 95.0, 105.0, 90.0, 98.0, 102.0, 88.0, 96.0]
    
    while current_pos < main_length:
        # Random event duration (30-150 samples)
        event_length = int(np.random.uniform(30, 150))
        if current_pos + event_length > main_length:
            event_length = main_length - current_pos
        
        # Pick a current level
        level = current_levels[len(main_signal) % len(current_levels)]
        
        # Generate event with noise
        event = np.random.normal(loc=level, scale=2.5, size=event_length)
        main_signal.extend(event)
        current_pos += event_length
    
    main_signal = np.array(main_signal[:main_length])
    
    # Tail region: noisy
    tail = np.random.normal(loc=110.0, scale=15.0, size=50)
    
    # Concatenate all regions
    signal = np.concatenate([adapter, main_signal, tail])
    
    return signal.astype(np.float32), adapter_length


def visualize_event_detection(raw_signal, events, adapter_length, 
                               trim_start=200, trim_end=10):
    """
    Create comprehensive visualization of event detection results.
    """
    fig, axes = plt.subplots(3, 1, figsize=(15, 10))
    
    # Color scheme
    adapter_color = 'lightcoral'
    tail_color = 'lightblue'
    signal_color = 'darkblue'
    event_color = 'green'
    
    # --- Panel 1: Full raw signal with regions highlighted ---
    ax1 = axes[0]
    ax1.plot(raw_signal, color=signal_color, linewidth=0.5, alpha=0.7, label='Raw Signal')
    
    # Highlight adapter region (first ~250 samples in our synthetic data)
    ax1.axvspan(0, adapter_length, alpha=0.3, color=adapter_color, label='Adapter Region')
    
    # Highlight trimmed regions (based on f5c defaults)
    ax1.axvspan(0, trim_start, alpha=0.5, color='red', label=f'Trimmed Start ({trim_start})')
    ax1.axvspan(len(raw_signal) - trim_end, len(raw_signal), alpha=0.5, 
                color='orange', label=f'Trimmed End ({trim_end})')
    
    ax1.set_xlabel('Sample Index', fontsize=11)
    ax1.set_ylabel('Current (pA)', fontsize=11)
    ax1.set_title('Raw Nanopore Signal with Adapter and Trimmed Regions', 
                  fontsize=13, fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # --- Panel 2: Detected events overlaid on signal ---
    ax2 = axes[1]
    # Plot raw signal (after trimming region)
    signal_start = trim_start
    signal_end = len(raw_signal) - trim_end
    ax2.plot(range(signal_start, signal_end), 
             raw_signal[signal_start:signal_end], 
             color=signal_color, linewidth=0.5, alpha=0.5, label='Signal')
    
    # Overlay detected events
    for i, event in enumerate(events):
        start_idx = int(event['start'])
        end_idx = int(event['start'] + event['length'])
        event_mean = event['mean']
        
        # Draw horizontal line for event mean
        ax2.plot([start_idx, end_idx], [event_mean, event_mean], 
                color=event_color, linewidth=2, alpha=0.8)
        
        # Draw vertical boundaries
        if i == 0 or i == len(events) - 1:
            ax2.axvline(start_idx, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    ax2.set_xlabel('Sample Index (after trimming)', fontsize=11)
    ax2.set_ylabel('Current (pA)', fontsize=11)
    ax2.set_title(f'Detected Events (n={len(events)}) Overlaid on Signal', 
                  fontsize=13, fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # --- Panel 3: Event statistics ---
    ax3 = axes[2]
    
    # Extract event properties
    event_means = [e['mean'] for e in events]
    event_lengths = [e['length'] for e in events]
    event_stdvs = [e['stdv'] for e in events]
    
    # Create subplots within panel 3
    ax3_1 = plt.subplot(3, 3, 7)
    ax3_1.hist(event_means, bins=30, color='green', alpha=0.7, edgecolor='black')
    ax3_1.set_xlabel('Mean Current (pA)', fontsize=9)
    ax3_1.set_ylabel('Count', fontsize=9)
    ax3_1.set_title('Event Mean Distribution', fontsize=10)
    ax3_1.grid(True, alpha=0.3)
    
    ax3_2 = plt.subplot(3, 3, 8)
    ax3_2.hist(event_lengths, bins=30, color='blue', alpha=0.7, edgecolor='black')
    ax3_2.set_xlabel('Duration (samples)', fontsize=9)
    ax3_2.set_ylabel('Count', fontsize=9)
    ax3_2.set_title('Event Duration Distribution', fontsize=10)
    ax3_2.grid(True, alpha=0.3)
    
    ax3_3 = plt.subplot(3, 3, 9)
    ax3_3.hist(event_stdvs, bins=30, color='orange', alpha=0.7, edgecolor='black')
    ax3_3.set_xlabel('Std Dev (pA)', fontsize=9)
    ax3_3.set_ylabel('Count', fontsize=9)
    ax3_3.set_title('Event Std Dev Distribution', fontsize=10)
    ax3_3.grid(True, alpha=0.3)
    
    # Hide the main ax3 (used for layout)
    ax3.axis('off')
    
    plt.tight_layout()
    return fig


def main():
    """Main execution function."""
    
    print("=" * 70)
    print("Nanopore Event Detection Example".center(70))
    print("=" * 70)
    
    # 1. Generate synthetic signal
    print("\n[1] Generating synthetic nanopore signal...")
    # raw_signal, adapter_length = generate_synthetic_signal(n_samples=2000, seed=42)
    with open('signal.pkl','rb') as f:
        raw_signal= pickle.load(f)
    adapter_length=200
    print(f"    ✓ Generated {len(raw_signal)} samples")
    print(f"    ✓ Adapter region: 0-{adapter_length} samples")
    print(f"    ✓ Signal range: {raw_signal.min():.2f} - {raw_signal.max():.2f} pA")
    
    # 2. Detect events
    print("\n[2] Detecting events using f5c algorithm...")
    print("    ✓ Trimming: start=200, end=10 (f5c defaults)")
    print("    ✓ Using RNA detection parameters")
    
    try:
        events = detect_events(raw_signal, is_rna=True)
        print(f"    ✓ Detected {len(events)} events")
    except Exception as e:
        print(f"    ✗ Error: {e}")
        return
    
    # 3. Display event statistics
    print("\n[3] Event Statistics:")
    print("-" * 70)
    
    if events:
        event_means = [e['mean'] for e in events]
        event_lengths = [e['length'] for e in events]
        event_stdvs = [e['stdv'] for e in events]
        
        print(f"    Total events:     {len(events)}")
        print(f"    Mean current:     {np.mean(event_means):.2f} ± {np.std(event_means):.2f} pA")
        print(f"    Event duration:   {np.mean(event_lengths):.1f} ± {np.std(event_lengths):.1f} samples")
        print(f"    Event std dev:    {np.mean(event_stdvs):.2f} ± {np.std(event_stdvs):.2f} pA")
        print(f"    First event:      start={events[0]['start']}, length={events[0]['length']:.1f}")
        print(f"    Last event:       start={events[-1]['start']}, length={events[-1]['length']:.1f}")
        
        # Check if adapter was successfully trimmed
        first_event_start = events[0]['start']
        if first_event_start >= 200:
            print(f"    ✓ Adapter successfully trimmed (first event at sample {first_event_start})")
        else:
            print(f"    ⚠ Warning: First event at sample {first_event_start} (expected ≥200)")
    else:
        print("    ✗ No events detected!")
        return
    
    # 4. Visualize results
    print("\n[4] Creating visualization...")
    fig = visualize_event_detection(raw_signal, events, adapter_length, 
                                     trim_start=200, trim_end=10)
    
    # Save figure
    output_path = 'event_detection_example.png'
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"    ✓ Figure saved to: {output_path}")
    
    plt.show()
    
    print("\n" + "=" * 70)
    print("Example completed successfully!".center(70))
    print("=" * 70)


if __name__ == "__main__":
    main()