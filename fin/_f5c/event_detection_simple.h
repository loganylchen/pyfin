/*
 * Simple event detection header
 * Standalone implementation adapted from f5c/scrappie
 */

#ifndef EVENT_DETECTION_SIMPLE_H
#define EVENT_DETECTION_SIMPLE_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

// Event structure
typedef struct
{
    float mean;
    float stdv;
    uint64_t start;
    float length;
} event_t;

// Event table structure
typedef struct
{
    size_t n;
    size_t start;
    size_t end;
    event_t *event;
} event_table;

// Raw signal table structure
typedef struct
{
    size_t n;
    size_t start;
    size_t end;
    float *raw;
} raw_table;

// Detector parameters structure
typedef struct
{
    size_t window_length1;
    size_t window_length2;
    float threshold1;
    float threshold2;
    float peak_height;
} detector_param_t;

// Main event detection function using f5c algorithm
// Takes raw signal table and detector parameters
event_table detect_events_simple(raw_table const rt, detector_param_t const edparam);

// High-level wrapper function for numpy arrays with adapter trimming
// Uses RNA defaults for event detection
// Events are returned in 3'->5' order (reversed) to match RNA pore transit direction
event_table getevents_simple(size_t nsample, float *rawptr);

// Free event table memory
void free_event_table(event_table *et);

#endif // EVENT_DETECTION_SIMPLE_H
