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
#include <stdbool.h>
#include "event_detection_simple.h"

// RNA detector parameters (from f5c)
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
// x: input array to calculate quantiles from
// n: length of x
// p: array of percentiles (e.g., 0.0, 0.5, 1.0), modified in place with results
// np: number of percentiles to calculate
static void quantilef(const float *x, size_t n, float *p, size_t np)
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
            p[i] = tmp[n - 1]; // Write result to p[i], not x[i]
        }
        else
        {
            float frac = pos - idx;
            p[i] = tmp[idx] * (1.0f - frac) + tmp[idx + 1] * frac; // Write result to p[i]
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

    // Update n to reflect trimmed range
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
// Cumulative sum: element i is sum up to but excluding element i
static void compute_sum_sumsq(const float *data, double *sum, double *sumsq, size_t d_length)
{
    assert(d_length > 0);

    sum[0] = 0.0;
    sumsq[0] = 0.0;
    for (size_t i = 0; i < d_length; i++)
    {
        sum[i + 1] = sum[i] + data[i];
        sumsq[i + 1] = sumsq[i] + data[i] * data[i];
    }
}

// Compute windowed t-statistic from summary information
// Uses two-sample t-test between adjacent windows
static float *compute_tstat(const double *sum, const double *sumsq, size_t d_length, size_t w_length)
{
    assert(d_length > 0);
    assert(w_length > 0);

    float *tstat = (float *)calloc(d_length, sizeof(float));
    if (!tstat)
        return NULL;

    const float eta = FLT_MIN; // Prevent division by zero
    const float w_lengthf = (float)w_length;

    // Quick return: t-test not defined for < 2 points or insufficient data
    if (d_length < 2 * w_length || w_length < 2)
    {
        return tstat;
    }

    // Fudge boundaries
    for (size_t i = 0; i < w_length; i++)
    {
        tstat[i] = 0.0f;
        tstat[d_length - i - 1] = 0.0f;
    }

    // Compute t-statistic for each position
    for (size_t i = w_length; i <= d_length - w_length; i++)
    {
        // First window: [i-w_length, i)
        double sum1 = sum[i];
        double sumsq1 = sumsq[i];
        if (i > w_length)
        {
            sum1 -= sum[i - w_length];
            sumsq1 -= sumsq[i - w_length];
        }

        // Second window: [i, i+w_length)
        float sum2 = (float)(sum[i + w_length] - sum[i]);
        float sumsq2 = (float)(sumsq[i + w_length] - sumsq[i]);

        // Calculate means
        float mean1 = sum1 / w_lengthf;
        float mean2 = sum2 / w_lengthf;

        // Combined variance from both windows
        float combined_var = sumsq1 / w_lengthf - mean1 * mean1 +
                             sumsq2 / w_lengthf - mean2 * mean2;

        // Prevent problem due to very small variances
        combined_var = fmaxf(combined_var, eta);

        // Student's t-statistic for two samples of equal size
        const float delta_mean = mean2 - mean1;
        tstat[i] = fabsf(delta_mean) / sqrtf(combined_var / w_lengthf);
    }

    return tstat;
}

// Detector state for stateful peak detection
typedef struct
{
    int DEF_PEAK_POS;
    float DEF_PEAK_VAL;
    float *signal;
    size_t signal_length;
    float threshold;
    size_t window_length;
    size_t masked_to;
    int peak_pos;
    float peak_value;
    bool valid_peak;
} Detector;
typedef Detector *DetectorPtr;

// Stateful short-long peak detector from f5c/scrappie
// Short detector dominates long detector to prevent false positives
static size_t *short_long_peak_detector(DetectorPtr short_detector,
                                        DetectorPtr long_detector,
                                        const float peak_height)
{
    assert(short_detector->signal_length == long_detector->signal_length);

    const size_t ndetector = 2;
    DetectorPtr detectors[] = {short_detector, long_detector};

    size_t *peaks = (size_t *)calloc(short_detector->signal_length, sizeof(size_t));
    if (!peaks)
        return NULL;

    size_t peak_count = 0;
    for (size_t i = 0; i < short_detector->signal_length; i++)
    {
        for (unsigned int k = 0; k < ndetector; k++)
        {
            DetectorPtr detector = detectors[k];

            // Skip if we've been masked out
            if (detector->masked_to >= i)
            {
                continue;
            }

            float current_value = detector->signal[i];

            if (detector->peak_pos == detector->DEF_PEAK_POS)
            {
                // CASE 1: No peak recorded yet
                if (current_value < detector->peak_value)
                {
                    // Record a deeper minimum
                    detector->peak_value = current_value;
                }
                else if (current_value - detector->peak_value > peak_height)
                {
                    // We've seen a qualifying maximum
                    detector->peak_value = current_value;
                    detector->peak_pos = i;
                }
            }
            else
            {
                // CASE 2: In an existing peak, waiting to see if it's good
                if (current_value > detector->peak_value)
                {
                    // Update the peak
                    detector->peak_value = current_value;
                    detector->peak_pos = i;
                }

                // Short detector dominates long detector
                if (detector == short_detector)
                {
                    if (detector->peak_value > detector->threshold)
                    {
                        long_detector->masked_to = detector->peak_pos + detector->window_length;
                        long_detector->peak_pos = long_detector->DEF_PEAK_POS;
                        long_detector->peak_value = long_detector->DEF_PEAK_VAL;
                        long_detector->valid_peak = false;
                    }
                }

                // Have we convinced ourselves we've seen a peak?
                if (detector->peak_value - current_value > peak_height &&
                    detector->peak_value > detector->threshold)
                {
                    detector->valid_peak = true;
                }

                // Check distance if this is a good peak
                if (detector->valid_peak &&
                    (i - detector->peak_pos) > detector->window_length / 2)
                {
                    // Emit the boundary and reset
                    peaks[peak_count] = detector->peak_pos;
                    peak_count++;
                    detector->peak_pos = detector->DEF_PEAK_POS;
                    detector->peak_value = current_value;
                    detector->valid_peak = false;
                }
            }
        }
    }

    return peaks;
}

// Create an event from raw signal boundaries
static event_t create_event(size_t start, size_t end, double const *sums,
                            double const *sumsqs, size_t nsample)
{
    assert(start < nsample);
    assert(end <= nsample);

    event_t event = {0};
    event.start = (uint64_t)start;
    event.length = (float)(end - start);
    event.mean = (float)(sums[end] - sums[start]) / event.length;

    const float deltasqr = (sumsqs[end] - sumsqs[start]);
    const float var = deltasqr / event.length - event.mean * event.mean;
    event.stdv = sqrtf(fmaxf(var, 0.0f));

    return event;
}

// Create event table from detected peaks
static event_table create_events(size_t const *peaks, double const *sums,
                                 double const *sumsqs, size_t nsample)
{
    event_table et = {0};

    // Count number of events found
    size_t n = 1;
    for (size_t i = 0; i < nsample; i++)
    {
        if (peaks[i] > 0 && peaks[i] < nsample)
        {
            n++;
        }
    }

    et.event = (event_t *)calloc(n, sizeof(event_t));
    if (!et.event)
        return et;

    et.n = n;
    et.end = et.n;

    // First event -- starts at zero
    et.event[0] = create_event(0, peaks[0], sums, sumsqs, nsample);

    // Other events -- peak[i-1] -> peak[i]
    for (size_t ev = 1; ev < n - 1; ev++)
    {
        et.event[ev] = create_event(peaks[ev - 1], peaks[ev], sums, sumsqs, nsample);
    }

    // Last event -- ends at nsample
    et.event[n - 1] = create_event(peaks[n - 2], nsample, sums, sumsqs, nsample);

    return et;
}

// Main event detection function using f5c algorithm
event_table detect_events_simple(raw_table const rt, detector_param_t const edparam)
{
    event_table et = {0};

    if (!rt.raw || rt.n == 0)
    {
        return et;
    }

    double *sums = (double *)calloc(rt.n + 1, sizeof(double));
    double *sumsqs = (double *)calloc(rt.n + 1, sizeof(double));

    if (!sums || !sumsqs)
    {
        free(sums);
        free(sumsqs);
        return et;
    }

    compute_sum_sumsq(rt.raw + rt.start, sums, sumsqs, rt.n);

    float *tstat1 = compute_tstat(sums, sumsqs, rt.n, edparam.window_length1);
    float *tstat2 = compute_tstat(sums, sumsqs, rt.n, edparam.window_length2);

    if (!tstat1 || !tstat2)
    {
        free(tstat1);
        free(tstat2);
        free(sums);
        free(sumsqs);
        return et;
    }

    // Initialize short detector
    Detector short_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .signal = tstat1,
        .signal_length = rt.n,
        .threshold = edparam.threshold1,
        .window_length = edparam.window_length1,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false};

    // Initialize long detector
    Detector long_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .signal = tstat2,
        .signal_length = rt.n,
        .threshold = edparam.threshold2,
        .window_length = edparam.window_length2,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false};

    size_t *peaks = short_long_peak_detector(&short_detector, &long_detector,
                                             edparam.peak_height);

    if (!peaks)
    {
        free(tstat2);
        free(tstat1);
        free(sumsqs);
        free(sums);
        return et;
    }

    et = create_events(peaks, sums, sumsqs, rt.n);

    // Adjust event start positions to account for trimming
    for (size_t i = 0; i < et.n; i++)
    {
        et.event[i].start += rt.start;
    }

    free(peaks);
    free(tstat2);
    free(tstat1);
    free(sumsqs);
    free(sums);

    return et;
}

// Reverse event table in place (for RNA which transits pore 3'->5')
static void reverse_event_table(event_table *et)
{
    if (!et || !et->event || et->n < 2)
        return;

    size_t i = 0;
    size_t j = et->n - 1;

    while (i < j)
    {
        // Swap events
        event_t tmp = et->event[i];
        et->event[i] = et->event[j];
        et->event[j] = tmp;
        i++;
        j--;
    }
}

// Wrapper function for numpy array with adapter trimming
// For RNA-seq: events are reversed to match the 5'->3' sequence direction
// This is used internally by eventalign
event_table getevents_simple(size_t nsample, float *rawptr)
{
    // Parameters from f5c/scrappie defaults
    // For RNA Direct Sequencing:
    //   Signal start (index 0) = 3' end of RNA (enters pore first)
    //   Signal end (index N-1) = 5' end of RNA (exits pore last)
    int trim_start = 200;      // Trim from START of signal (3' end of RNA, adapter region)
    int trim_end = 10;         // Trim from END of signal (5' end of RNA, minimal trim)
    int varseg_chunk = 100;    // Chunk size for MAD calculation
    float varseg_thresh = 0.0; // Percentile threshold (0.0 = median)

    // Create raw table
    raw_table rt = {nsample, 0, nsample, rawptr};

    // Trim adapters and segment the signal
    rt = trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);

    // Use RNA detector parameters (only RNA supported)
    detector_param_t const *ed_params = &RNA_DEFAULTS;

    // Detect events on the trimmed signal
    event_table et = detect_events_simple(rt, *ed_params);

    // Reverse events for RNA (RNA transits pore 3'->5', so first event = 3' end)
    // After reversal: first event = 5' end, matching the 5'->3' sequence direction
    // This means users should NOT reverse their sequence - just use original 5'->3' sequence
    reverse_event_table(&et);

    return et;
}

// Wrapper function that returns events in RAW SIGNAL order (no reversal)
// Event[0] = first event in raw signal (3' end for RNA)
// Event[n-1] = last event in raw signal (5' end for RNA)
// This matches the order expected when indexing into the raw signal
event_table getevents_raw_order(size_t nsample, float *rawptr)
{
    // Parameters from f5c/scrappie defaults
    // Signal coordinates: start=3' end, end=5' end (for RNA)
    int trim_start = 200; // Trim 3' end (adapter)
    int trim_end = 10;    // Trim 5' end (minimal)
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;

    // Create raw table
    raw_table rt = {nsample, 0, nsample, rawptr};

    // Trim adapters and segment the signal
    rt = trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);

    // Use RNA detector parameters
    detector_param_t const *ed_params = &RNA_DEFAULTS;

    // Detect events on the trimmed signal - NO REVERSAL
    event_table et = detect_events_simple(rt, *ed_params);

    return et;
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
