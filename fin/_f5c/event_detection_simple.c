/*
 * Simple event detection implementation
 * Standalone version adapted from scrappie/f5c without external dependencies
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include "event_detection_simple.h"

// Detector parameters
typedef struct
{
    size_t window_length1;
    size_t window_length2;
    float threshold1;
    float threshold2;
    float peak_height;
} detector_param_t;

static detector_param_t RNA_DEFAULTS = {
    .window_length1 = 7,
    .window_length2 = 14,
    .threshold1 = 2.5f,
    .threshold2 = 9.0f,
    .peak_height = 1.0f};

// Helper function for sorting floats
static int compare_float(const void *a, const void *b)
{
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

// Calculate median of float array
static float medianf(float *x, size_t n)
{
    if (n == 0)
        return 0.0f;

    float *tmp = (float *)malloc(n * sizeof(float));
    if (!tmp)
        return 0.0f;

    memcpy(tmp, x, n * sizeof(float));
    qsort(tmp, n, sizeof(float), compare_float);

    float result;
    if (n % 2 == 0)
        result = (tmp[n / 2 - 1] + tmp[n / 2]) / 2.0f;
    else
        result = tmp[n / 2];

    free(tmp);
    return result;
}

// Calculate quantile of float array
static void quantilef(float *x, size_t n, const float *p, size_t np)
{
    if (n == 0 || np == 0)
        return;

    float *tmp = (float *)malloc(n * sizeof(float));
    if (!tmp)
        return;

    memcpy(tmp, x, n * sizeof(float));
    qsort(tmp, n, sizeof(float), compare_float);

    for (size_t i = 0; i < np; i++)
    {
        if (p[i] < 0.0f || p[i] > 1.0f)
            continue;

        float pos = p[i] * (n - 1);
        size_t idx = (size_t)pos;

        if (idx >= n - 1)
        {
            x[i] = tmp[n - 1];
        }
        else
        {
            float frac = pos - idx;
            x[i] = tmp[idx] * (1.0f - frac) + tmp[idx + 1] * frac;
        }
    }

    free(tmp);
}

// Calculate Median Absolute Deviation (MAD)
static float madf(float const *x, size_t n, float const *med)
{
    if (n == 0)
        return 0.0f;

    // MAD scaling factor for normal distribution
    const float mad_scaling_factor = 1.4826f;

    float _med = (med == NULL) ? medianf((float *)x, n) : *med;

    float *absdiff = (float *)malloc(n * sizeof(float));
    if (!absdiff)
        return 0.0f;

    for (size_t i = 0; i < n; i++)
    {
        absdiff[i] = fabsf(x[i] - _med);
    }

    float mad = medianf(absdiff, n);
    free(absdiff);

    return mad * mad_scaling_factor;
}

// Trim raw signal by MAD-based segmentation
// This removes adapter regions with low signal variation
static raw_table trim_raw_by_mad(raw_table rt, int chunk_size, float perc)
{
    assert(chunk_size > 1);
    assert(perc >= 0.0f && perc <= 1.0f);

    const size_t nsample = rt.end - rt.start;
    const size_t nchunk = nsample / chunk_size;

    // Truncate end to be consistent with f5c/scrappie
    rt.end = rt.start + nchunk * chunk_size;

    if (nchunk == 0)
        return rt;

    float *madarr = (float *)malloc(nchunk * sizeof(float));
    if (!madarr)
        return rt;

    // Calculate MAD for each chunk
    for (size_t i = 0; i < nchunk; i++)
    {
        madarr[i] = madf(rt.raw + rt.start + i * chunk_size, chunk_size, NULL);
    }

    // Calculate threshold as percentile of MADs
    quantilef(madarr, nchunk, &perc, 1);
    const float thresh = perc;

    // Trim from start: remove chunks with MAD <= threshold
    for (size_t i = 0; i < nchunk; i++)
    {
        if (madarr[i] > thresh)
            break;
        rt.start += chunk_size;
    }

    // Trim from end: remove chunks with MAD <= threshold
    for (size_t i = nchunk; i > 0; i--)
    {
        if (madarr[i - 1] > thresh)
            break;
        rt.end -= chunk_size;
    }

    free(madarr);

    // Ensure we still have valid data
    if (rt.start >= rt.end)
    {
        rt.start = 0;
        rt.end = 0;
    }

    return rt;
}

// Trim and segment raw signal to remove adapters and open pore regions
static raw_table trim_and_segment_raw(raw_table rt, int trim_start, int trim_end,
                                      int varseg_chunk, float varseg_thresh)
{
    if (!rt.raw)
        return rt;

    // First, apply MAD-based segmentation to remove adapter regions
    rt = trim_raw_by_mad(rt, varseg_chunk, varseg_thresh);
    if (!rt.raw)
        return rt;

    // Then apply additional fixed trimming
    rt.start += trim_start;
    rt.end -= trim_end;

    // Ensure we have valid range
    if (rt.start >= rt.end)
    {
        rt.start = 0;
        rt.end = 0;
        rt.n = 0;
    }
    else
    {
        rt.n = rt.end - rt.start;
    }

    return rt;
}

// Compute sum and sum of squares for the raw signal
static void compute_sum_sumsq(float const *raw, size_t start, double *sums, double *sumsqs, size_t n)
{
    sums[0] = 0.0;
    sumsqs[0] = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        float val = raw[start + i];
        sums[i + 1] = sums[i] + val;
        sumsqs[i + 1] = sumsqs[i] + val * val;
    }
}

// Compute t-statistic for a given window length
static float *compute_tstat(double const *sums, double const *sumsqs, size_t n, size_t window_length)
{
    float *tstat = (float *)calloc(n, sizeof(float));
    if (!tstat)
        return NULL;

    double inv_winlen = 1.0 / window_length;
    double var_winlen = (double)window_length;

    for (size_t i = 0; i < window_length; i++)
    {
        tstat[i] = 0.0f;
    }

    for (size_t i = window_length; i < n; i++)
    {
        double sum_X = sums[i] - sums[i - window_length];
        double sum_X2 = sumsqs[i] - sumsqs[i - window_length];
        double mean = sum_X * inv_winlen;
        double var = sum_X2 * inv_winlen - mean * mean;
        double t = (mean - 0.0) / sqrt(var / var_winlen);
        tstat[i] = (float)fabs(t);
    }

    return tstat;
}

// Simple peak detector for event boundaries
static size_t *peak_detector(float const *signal, size_t n, float threshold,
                             float peak_height, size_t window_length, size_t *n_peaks)
{
    // Conservative estimate: at most one peak per window_length
    size_t max_peaks = n / window_length + 10;
    size_t *peaks = (size_t *)calloc(max_peaks, sizeof(size_t));
    if (!peaks)
    {
        *n_peaks = 0;
        return NULL;
    }

    size_t peak_count = 0;

    // Find local maxima above threshold
    for (size_t i = window_length; i < n - window_length; i++)
    {
        if (signal[i] > threshold)
        {
            // Check if this is a peak
            bool is_peak = true;
            float current = signal[i];

            for (size_t j = 1; j <= window_length / 2; j++)
            {
                if (signal[i - j] >= current || signal[i + j] >= current)
                {
                    is_peak = false;
                    break;
                }
            }

            if (is_peak && current > peak_height)
            {
                peaks[peak_count++] = i;

                // Skip ahead to avoid detecting the same event multiple times
                i += window_length;
            }
        }
    }

    *n_peaks = peak_count;
    return peaks;
}

// Create an event from raw signal boundaries
static event_t create_event(size_t start, size_t end, double const *sums,
                            double const *sumsqs)
{
    assert(start < end);

    event_t event = {0};
    event.length = (float)(end - start);
    event.start = start;

    float sum = (float)(sums[end] - sums[start]);
    event.mean = sum / event.length;

    float sum_sq = (float)(sumsqs[end] - sumsqs[start]);
    float var = sum_sq / event.length - event.mean * event.mean;
    event.stdv = sqrtf(fmaxf(var, 0.0f));

    return event;
}

// Simple dual-window peak detector (short and long)
static size_t *detect_peaks(float const *tstat1, float const *tstat2, size_t n,
                            float peak_height, size_t short_win, size_t long_win,
                            size_t *n_peaks)
{
    // Use a simple approach: detect peaks in both statistics
    // and combine them

    size_t n_peaks1, n_peaks2;
    size_t *peaks1 = peak_detector(tstat1, n, 2.0f, peak_height, short_win, &n_peaks1);
    size_t *peaks2 = peak_detector(tstat2, n, 3.0f, peak_height, long_win, &n_peaks2);

    if (!peaks1 || !peaks2)
    {
        free(peaks1);
        free(peaks2);
        *n_peaks = 0;
        return NULL;
    }

    // Conservative estimate for merged peaks
    size_t max_peaks = n_peaks1 + n_peaks2 + 10;
    size_t *merged_peaks = (size_t *)calloc(max_peaks, sizeof(size_t));
    if (!merged_peaks)
    {
        free(peaks1);
        free(peaks2);
        *n_peaks = 0;
        return NULL;
    }

    // Copy peaks from short window detector
    size_t merged_count = 0;
    for (size_t i = 0; i < n_peaks1; i++)
    {
        merged_peaks[merged_count++] = peaks1[i];
    }

    // Add peaks from long window detector that aren't too close
    for (size_t i = 0; i < n_peaks2; i++)
    {
        size_t peak_pos = peaks2[i];
        bool too_close = false;

        // Check if this peak is within short window of an existing peak
        for (size_t j = 0; j < merged_count; j++)
        {
            if (abs((int)peak_pos - (int)merged_peaks[j]) < (int)short_win)
            {
                too_close = true;
                break;
            }
        }

        if (!too_close)
        {
            merged_peaks[merged_count++] = peak_pos;
        }
    }

    // Sort peaks by position
    for (size_t i = 0; i < merged_count - 1; i++)
    {
        for (size_t j = 0; j < merged_count - i - 1; j++)
        {
            if (merged_peaks[j] > merged_peaks[j + 1])
            {
                size_t temp = merged_peaks[j];
                merged_peaks[j] = merged_peaks[j + 1];
                merged_peaks[j + 1] = temp;
            }
        }
    }

    free(peaks1);
    free(peaks2);

    *n_peaks = merged_count;
    return merged_peaks;
}

// Main event detection function
event_table detect_events_simple(raw_table const rt, int is_rna)
{
    event_table et = {0};

    if (!rt.raw || rt.n == 0)
    {
        return et;
    }

    // Calculate the working range
    size_t work_start = rt.start;
    size_t work_end = rt.end;
    size_t work_n = work_end - work_start;

    if (work_n == 0)
    {
        return et;
    }

    // Allocate memory for sums and sum squares
    double *sums = (double *)calloc(work_n + 1, sizeof(double));
    double *sumsqs = (double *)calloc(work_n + 1, sizeof(double));

    if (!sums || !sumsqs)
    {
        free(sums);
        free(sumsqs);
        return et;
    }

    // Compute sums from the working range
    compute_sum_sumsq(rt.raw, work_start, sums, sumsqs, work_n);

    // Get detector parameters
    detector_param_t *params = is_rna ? &RNA_DEFAULTS : &DNA_DEFAULTS;

    // Compute t-statistics for both window sizes
    float *tstat1 = compute_tstat(sums, sumsqs, work_n, params->window_length1);
    float *tstat2 = compute_tstat(sums, sumsqs, work_n, params->window_length2);

    if (!tstat1 || !tstat2)
    {
        free(tstat1);
        free(tstat2);
        free(sums);
        free(sumsqs);
        return et;
    }

    // Detect peaks
    size_t n_peaks;
    size_t *peaks = detect_peaks(tstat1, tstat2, work_n, params->peak_height,
                                 params->window_length1, params->window_length2,
                                 &n_peaks);

    if (!peaks || n_peaks == 0)
    {
        free(peaks);
        free(tstat2);
        free(tstat1);
        free(sumsqs);
        free(sums);
        return et;
    }

    // Create events from peaks
    // First event starts at 0 and goes to first peak
    size_t max_events = n_peaks + 2;
    et.event = (event_t *)calloc(max_events, sizeof(event_t));
    if (!et.event)
    {
        free(peaks);
        free(tstat2);
        free(tstat1);
        free(sumsqs);
        free(sums);
        return et;
    }

    size_t event_idx = 0;

    // First event (0 to first peak)
    et.event[event_idx++] = create_event(0, peaks[0], sums, sumsqs);

    // Middle events (between peaks)
    for (size_t i = 0; i < n_peaks - 1 && event_idx < max_events; i++)
    {
        et.event[event_idx++] = create_event(peaks[i], peaks[i + 1], sums, sumsqs);
    }

    // Last event (last peak to end)
    et.event[event_idx++] = create_event(peaks[n_peaks - 1], work_n, sums, sumsqs);

    et.n = event_idx;
    et.start = 0;
    et.end = et.n;

    // Adjust event start positions to account for trimming
    for (size_t i = 0; i < et.n; i++)
    {
        et.event[i].start += work_start;
    }

    // Cleanup
    free(peaks);
    free(tstat2);
    free(tstat1);
    free(sumsqs);
    free(sums);

    return et;
}

// Wrapper function for numpy array with adapter trimming
event_table getevents_simple(size_t nsample, float *rawptr, int is_rna)
{
    // Parameters from f5c/scrappie defaults
    int trim_start = 200;      // Trim 200 samples from start
    int trim_end = 10;         // Trim 10 samples from end
    int varseg_chunk = 100;    // Chunk size for MAD calculation
    float varseg_thresh = 0.0; // Percentile threshold (0.0 = median)

    // Create raw table
    raw_table rt = {nsample, 0, nsample, rawptr};

    // Trim adapters and segment the signal
    rt = trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);

    // Detect events on the trimmed signal
    return detect_events_simple(rt, is_rna);
}

void free_event_table(event_table *et)
{
    if (et)
    {
        free(et->event);
        et->event = NULL;
        et->n = 0;
    }
}
