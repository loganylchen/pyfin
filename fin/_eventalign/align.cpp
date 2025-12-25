/* @file align.cpp
**
** CPU implementation of the Adaptive banded Event Alignment algorithm
**
** This file was adapted from align.c in the f5c project
** @author: pyfin
** @@
******************************************************************************/

#include "common.h"
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

// =============================================================================
// Helper functions
// =============================================================================

//todo : can make more efficient using bit encoding
static inline uint32_t get_rank(char base) {
    if (base == 'A') { //todo: do we neeed simple alpha?
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else {
        return 0;  // For 'U' and 'N'
    }
}

// return the lexicographic rank of the kmer amongst all strings of
// length k for this alphabet
static inline uint32_t get_kmer_rank(const char* str, uint32_t k) {
    uint32_t r = 0;
    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

//copy a kmer from a reference
static inline void kmer_cpy(char* dest, const char* src, uint32_t k) {
    uint32_t i = 0;
    for (i = 0; i < k; i++) {
        dest[i] = src[i];
    }
    dest[i] = '\0';
}

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// =============================================================================
// Alignment functions (extern "C" for Python wrapper)
// =============================================================================

extern "C" {

scalings_t estimate_scalings_using_mom(const char* sequence, int32_t sequence_len,
                                        const model_t* pore_model, uint32_t kmer_size,
                                        const event_table& et) {
    scalings_t out;
    int32_t n_kmers = sequence_len - kmer_size + 1;

    // Calculate summary statistics over the events and the model implied by the read
    double event_level_sum = 0.0f;
    for (size_t i = 0; i < et.n; ++i) {
        event_level_sum += et.event[i].mean;
    }

    double kmer_level_sum = 0.0f;
    double kmer_level_sq_sum = 0.0f;
    for (int32_t i = 0; i < n_kmers; ++i) {
        int32_t kr = get_kmer_rank(&sequence[i], kmer_size);
        double l = pore_model[kr].level_mean;
        kmer_level_sum += l;
        kmer_level_sq_sum += l * l;
    }

    double shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // estimate scale
    double event_level_sq_sum = 0.0f;
    for (size_t i = 0; i < et.n; ++i) {
        event_level_sq_sum += (et.event[i].mean - shift) * (et.event[i].mean - shift);
    }

    double scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

    out.shift = (float)shift;
    out.scale = (float)scale;
    out.var = 1.0f;
#ifdef CACHED_LOG
    out.log_var = 0.0f;
#endif

    return out;
}

static inline float log_normal_pdf(float x, float gp_mean, float gp_stdv, float gp_log_stdv) {
    float log_inv_sqrt_2pi = -0.918938f; // Natural logarithm
    float a = (x - gp_mean) / gp_stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
}

static inline float log_probability_match_r9(scalings_t scaling,
                                             const model_t* models,
                                             const event_table& events, int event_idx,
                                             uint32_t kmer_rank, uint8_t strand,
                                             float sample_rate) {
    (void)strand;  // unused
    (void)sample_rate;  // unused

    float unscaledLevel = events.event[event_idx].mean;
    float scaledLevel = unscaledLevel;

    float gp_mean = scaling.scale * models[kmer_rank].level_mean + scaling.shift;
    float gp_stdv = models[kmer_rank].level_stdv;

#ifdef CACHED_LOG
    float gp_log_stdv = models[kmer_rank].level_log_stdv;
#else
    float gp_log_stdv = log(models[kmer_rank].level_stdv);
#endif

    float lp = log_normal_pdf(scaledLevel, gp_mean, gp_stdv, gp_log_stdv);
    return lp;
}

#define event_kmer_to_band(ei, ki) ((ei) + 1) + (ki + 1)
#define band_event_to_offset(bi, ei, band_lower_left, bandwidth) \
    (band_lower_left[(bi)].event_idx - (ei))
#define band_kmer_to_offset(bi, ki, band_lower_left) ((ki) - band_lower_left[(bi)].kmer_idx)
#define is_offset_valid(offset, bandwidth) ((offset) >= 0 && (offset) < (bandwidth))
#define event_at_offset(bi, offset, band_lower_left) (band_lower_left[(bi)].event_idx - (offset))
#define kmer_at_offset(bi, offset, band_lower_left) (band_lower_left[(bi)].kmer_idx + (offset))

#define move_down(curr_band) \
    { (curr_band).event_idx + 1, (curr_band).kmer_idx }
#define move_right(curr_band) \
    { (curr_band).event_idx, (curr_band).kmer_idx + 1 }

#define BAND_ARRAY(bands, r, c, bandwidth) (bands[((r) * ((bandwidth))) + (c)])
#define TRACE_ARRAY(trace, r, c, bandwidth) (trace[((r) * ((bandwidth))) + (c)])

int32_t align(AlignedPair* out_2, const char* sequence, int32_t sequence_len,
              const event_table& events, const model_t* models, uint32_t kmer_size,
              scalings_t scaling, float sample_rate) {

    (void)sample_rate;  // unused
    size_t strand_idx = 0;
    size_t k = kmer_size;

    size_t n_events = events.n;
    size_t n_kmers = sequence_len - k + 1;

    // backtrack markers
    const uint8_t FROM_D = 0;
    const uint8_t FROM_U = 1;
    const uint8_t FROM_L = 2;

    // qc
    double min_average_log_emission = -5.0;
    int max_gap_threshold = 50;

    // banding
    int bandwidth = ALN_BANDWIDTH;
    int half_bandwidth = ALN_BANDWIDTH / 2;

    // transition penalties
    double events_per_kmer = (double)n_events / n_kmers;
    double p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);

    // dp matrix
    size_t n_rows = n_events + 1;
    size_t n_cols = n_kmers + 1;
    size_t n_bands = n_rows + n_cols;

    // Precompute k-mer ranks
    size_t* kmer_ranks = (size_t*)malloc(sizeof(size_t) * n_kmers);
    if (!kmer_ranks) return 0;

    for (size_t i = 0; i < n_kmers; ++i) {
        kmer_ranks[i] = get_kmer_rank(&sequence[i], k);
    }

    float* bands = (float*)malloc(sizeof(float) * n_bands * bandwidth);
    uint8_t* trace = (uint8_t*)malloc(sizeof(uint8_t) * n_bands * bandwidth);

    if (!bands || !trace) {
        free(kmer_ranks);
        if (bands) free(bands);
        if (trace) free(trace);
        return 0;
    }

    for (size_t i = 0; i < n_bands * bandwidth; i++) {
        bands[i] = -INFINITY;
        trace[i] = 0;
    }

    // Keep track of the event/kmer index for the lower left corner of the band
    EventKmerPair* band_lower_left = (EventKmerPair*)malloc(sizeof(EventKmerPair) * n_bands);
    if (!band_lower_left) {
        free(kmer_ranks);
        free(bands);
        free(trace);
        return 0;
    }

    // initialize range of first two
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    int start_cell_offset = band_kmer_to_offset(0, -1, band_lower_left);
    BAND_ARRAY(bands, 0, start_cell_offset, bandwidth) = 0.0f;

    // band 1: first event is trimmed
    int first_trim_offset = band_event_to_offset(1, 0, band_lower_left, bandwidth);
    BAND_ARRAY(bands, 1, first_trim_offset, bandwidth) = lp_trim;
    TRACE_ARRAY(trace, 1, first_trim_offset, bandwidth) = FROM_U;

    // fill in remaining bands
    for (size_t band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        float ll = BAND_ARRAY(bands, band_idx - 1, 0, bandwidth);
        float ur = BAND_ARRAY(bands, band_idx - 1, bandwidth - 1, bandwidth);
        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;

        bool right = false;
        if (ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if (right) {
            band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1]);
        } else {
            band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1]);
        }

        // If the trim state is within the band, fill it in here
        int trim_offset = band_kmer_to_offset(band_idx, -1, band_lower_left);
        if (is_offset_valid(trim_offset, bandwidth)) {
            int64_t event_idx = event_at_offset(band_idx, trim_offset, band_lower_left);
            if (event_idx >= 0 && event_idx < (int64_t)n_events) {
                BAND_ARRAY(bands, band_idx, trim_offset, bandwidth) = lp_trim * (event_idx + 1);
                TRACE_ARRAY(trace, band_idx, trim_offset, bandwidth) = FROM_U;
            } else {
                BAND_ARRAY(bands, band_idx, trim_offset, bandwidth) = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        int kmer_min_offset = band_kmer_to_offset(band_idx, 0, band_lower_left);
        int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers, band_lower_left);
        int event_min_offset = band_event_to_offset(band_idx, n_events - 1, band_lower_left, bandwidth);
        int event_max_offset = band_event_to_offset(band_idx, -1, band_lower_left, bandwidth);

        int min_offset = MAX(kmer_min_offset, event_min_offset);
        min_offset = MAX(min_offset, 0);

        int max_offset = MIN(kmer_max_offset, event_max_offset);
        max_offset = MIN(max_offset, bandwidth);

        for (int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = event_at_offset(band_idx, offset, band_lower_left);
            int kmer_idx = kmer_at_offset(band_idx, offset, band_lower_left);

            size_t kmer_rank = kmer_ranks[kmer_idx];

            int offset_up = band_event_to_offset(band_idx, event_idx - 1, band_lower_left, bandwidth);
            int offset_left = band_kmer_to_offset(band_idx, kmer_idx - 1, band_lower_left);
            int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1, band_lower_left);

            float up = is_offset_valid(offset_up, bandwidth) ?
                           BAND_ARRAY(bands, band_idx - 1, offset_up, bandwidth) : -INFINITY;
            float left = is_offset_valid(offset_left, bandwidth) ?
                             BAND_ARRAY(bands, band_idx - 1, offset_left, bandwidth) : -INFINITY;
            float diag = is_offset_valid(offset_diag, bandwidth) ?
                             BAND_ARRAY(bands, band_idx - 2, offset_diag, bandwidth) : -INFINITY;

            float lp_emission = log_probability_match_r9(scaling, models, events, event_idx,
                                                          kmer_rank, strand_idx, 4000.0f);
            float score_d = diag + lp_step + lp_emission;
            float score_u = up + lp_stay + lp_emission;
            float score_l = left + lp_skip;

            float max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;

            BAND_ARRAY(bands, band_idx, offset, bandwidth) = max_score;
            TRACE_ARRAY(trace, band_idx, offset, bandwidth) = from;
        }
    }

    //
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;

    int outIndex = 0;

    float max_score = -INFINITY;
    int curr_event_idx = 0;
    int curr_kmer_idx = n_kmers - 1;

    // Find best score between an event and the last k-mer
    for (size_t event_idx = 0; event_idx < n_events; ++event_idx) {
        int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
        int offset = band_event_to_offset(band_idx, event_idx, band_lower_left, bandwidth);
        if (is_offset_valid(offset, bandwidth)) {
            float s = BAND_ARRAY(bands, band_idx, offset, bandwidth) +
                      (n_events - event_idx) * lp_trim;
            if (s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }

    int curr_gap = 0;
    int max_gap = 0;
    while (curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        // emit alignment
        out_2[outIndex].ref_pos = curr_kmer_idx;
        out_2[outIndex].read_pos = curr_event_idx;
        outIndex++;

        // qc stats
        size_t kmer_rank = get_kmer_rank(&sequence[curr_kmer_idx], k);
        float tempLogProb = log_probability_match_r9(scaling, models, events, curr_event_idx,
                                                      kmer_rank, 0, 4000.0f);
        sum_emission += tempLogProb;
        n_aligned_events += 1;

        int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int offset = band_event_to_offset(band_idx, curr_event_idx, band_lower_left, bandwidth);

        uint8_t from = TRACE_ARRAY(trace, band_idx, offset, bandwidth);
        if (from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if (from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = MAX(curr_gap, max_gap);
        }
    }

    // Reverse the output
    int c;
    int end = outIndex - 1;
    for (c = 0; c < outIndex / 2; c++) {
        int ref_pos_temp = out_2[c].ref_pos;
        int read_pos_temp = out_2[c].read_pos;
        out_2[c].ref_pos = out_2[end].ref_pos;
        out_2[c].read_pos = out_2[end].read_pos;
        out_2[end].ref_pos = ref_pos_temp;
        out_2[end].read_pos = read_pos_temp;
        end--;
    }

    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    bool spanned = out_2[0].ref_pos == 0 && out_2[outIndex - 1].ref_pos == (int)(n_kmers - 1);

    // Debug output for alignment QC failure
    if (avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
        fprintf(stderr, "[align] Alignment QC failed:\n");
        fprintf(stderr, "  avg_log_emission: %.4f (threshold: %.4f) %s\n",
                avg_log_emission, min_average_log_emission,
                avg_log_emission < min_average_log_emission ? "FAIL" : "OK");
        fprintf(stderr, "  spanned: %s (first_ref_pos: %d, last_ref_pos: %d, n_kmers-1: %d) %s\n",
                spanned ? "true" : "false", out_2[0].ref_pos, out_2[outIndex - 1].ref_pos, (int)(n_kmers - 1),
                !spanned ? "FAIL" : "OK");
        fprintf(stderr, "  max_gap: %d (threshold: %d) %s\n",
                max_gap, max_gap_threshold,
                max_gap > max_gap_threshold ? "FAIL" : "OK");
        outIndex = 0;
    }

    free(kmer_ranks);
    free(bands);
    free(trace);
    free(band_lower_left);

    return outIndex;
}

int32_t postalign(event_alignment_t* alignment, index_pair_t* base_to_event_map,
                  double* events_per_base, const char* sequence, int32_t n_kmers,
                  const AlignedPair* event_alignment, int32_t n_events, uint32_t kmer_size) {
    // create base-to-event map
    int32_t i = 0;
    for (i = 0; i < n_kmers; i++) {
        base_to_event_map[i].start = -1;
        base_to_event_map[i].stop = -1;
    }

    int32_t max_event = 0;
    int32_t min_event = INT32_MAX;

    int32_t prev_event_idx = -1;

    for (i = 0; i < n_events; ++i) {
        int32_t k_idx = event_alignment[i].ref_pos;
        int32_t event_idx = event_alignment[i].read_pos;
        index_pair_t* elem = &base_to_event_map[k_idx];

        if (event_idx != prev_event_idx) {
            if (elem->start == -1) {
                elem->start = event_idx;
            }
            elem->stop = event_idx;
        }
        max_event = max_event > event_idx ? max_event : event_idx;
        min_event = min_event < event_idx ? min_event : event_idx;
        prev_event_idx = event_idx;
    }

    *events_per_base = (double)(max_event - min_event) / n_kmers;

    // prepare data structures for the final calibration
    int32_t alignment_index = 0;
    int32_t prev_kmer_rank = -1;

    int32_t ki;
    for (ki = 0; ki < n_kmers; ++ki) {
        index_pair_t event_range_for_kmer = base_to_event_map[ki];

        // skip kmers without events
        if (event_range_for_kmer.start == -1) {
            continue;
        }

        for (int32_t event_idx = event_range_for_kmer.start;
             event_idx <= event_range_for_kmer.stop; event_idx++) {

            int32_t kmer_rank = get_kmer_rank(&sequence[ki], kmer_size);

            event_alignment_t ea;
            ea.read_idx = -1;
            kmer_cpy(ea.ref_kmer, &sequence[ki], kmer_size);
            ea.ref_position = ki;
            ea.event_idx = event_idx;
            ea.rc = false;
            kmer_cpy(ea.model_kmer, &sequence[ki], kmer_size);
            ea.hmm_state = prev_kmer_rank != kmer_rank ? 'M' : 'E';

            if (alignment_index > n_events) {
                return alignment_index;
            }
            alignment[alignment_index] = ea;
            alignment_index++;
            prev_kmer_rank = kmer_rank;
        }
    }

    return alignment_index;
}

}  // extern "C"
