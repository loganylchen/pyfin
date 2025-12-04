/*
 * Simple event detection header
 * Standalone implementation without external dependencies
 */

#ifndef EVENT_DETECTION_SIMPLE_H
#define EVENT_DETECTION_SIMPLE_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

// Forward declare the event structure
typedef struct
{
    float mean;
    float stdv;
    uint64_t start;
    float length;
} event_t;

typedef struct
{
    size_t n;
    size_t start;
    size_t end;
    event_t *event;
} event_table;

typedef struct
{
    size_t n;
    size_t start;
    size_t end;
    float *raw;
} raw_table;

// Event detection function
event_table detect_events_simple(raw_table const rt, int is_rna);

void getevents_simple(size_t nsample, float *rawptr, int is_rna);

// Free event table
void free_event_table(event_table *et);

#endif // EVENT_DETECTION_SIMPLE_H
