/* @file align_cuda.cu
**
** GPU implementation of the Adaptive banded Event Alignment algorithm
** Based on f5c CUDA implementation, adapted for pyfin
**
** This file contains the CUDA kernels for event alignment:
** - reverse_events_kernel: Reverses events for RNA (3'->5' to 5'->3')
** - align_kernel_pre_2d: K-mer rank precomputation and band initialization
** - align_kernel_core_2d_shm: Main adaptive banded alignment (Suzuki algorithm)
** - align_kernel_post: Backtracking and output generation
** @@
******************************************************************************/

#define HAVE_CUDA 1

#include "align_cuda.h"
#include <assert.h>
#include <math.h>
#include <cuda_runtime.h>

// CUDA error checking macro (redefined from header to avoid conflicts)
#undef CUDA_CHK
#define CUDA_CHK()                                                             \
    {                                                                           \
        cudaError_t code = cudaGetLastError();                                  \
        if (code != cudaSuccess) {                                              \
            fprintf(stderr, "[CUDA ERROR] %s:%d: %s\n", __FILE__, __LINE__,     \
                    cudaGetErrorString(code));                                 \
            exit(-1);                                                           \
        }                                                                       \
    }

// Device functions for k-mer operations
__forceinline__ __device__ uint32_t get_rank(char base) {
    if (base == 'A') {
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else {
        return 0;
    }
}

__forceinline__ __device__ uint32_t get_kmer_rank(const char* str, uint32_t k) {
    uint32_t r = 0;
    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

#define log_inv_sqrt_2pi  -0.918938f // Natural logarithm

__forceinline__ __device__ float log_normal_pdf(float x, float gp_mean, float gp_stdv, float gp_log_stdv) {
    float a = (x - gp_mean) / gp_stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
}

__forceinline__ __device__ float log_probability_match_r9(scalings_t scaling, model_t* models, event_t* event,
                                                            int event_idx, uint32_t kmer_rank) {
    float unscaledLevel = event[event_idx].mean;
    float scaledLevel = unscaledLevel;

    float gp_mean = scaling.scale * models[kmer_rank].level_mean + scaling.shift;
    float gp_stdv = models[kmer_rank].level_stdv * 1;

#ifdef CACHED_LOG
    float gp_log_stdv = models[kmer_rank].level_log_stdv;
#else
    float gp_log_stdv = logf(models[kmer_rank].level_stdv);
#endif

    float lp = log_normal_pdf(scaledLevel, gp_mean, gp_stdv, gp_log_stdv);
    return lp;
}

// Band alignment macros
#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[(bi)].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[(bi)].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

#define move_down(curr_band) { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) { curr_band.event_idx, curr_band.kmer_idx + 1 }

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define BAND_ARRAY(r, c) ( bands[((r)*(ALN_BANDWIDTH)+(c))] )
#define TRACE_ARRAY(r, c) ( trace[((r)*(ALN_BANDWIDTH)+(c))] )

#define FROM_D  0
#define FROM_U  1
#define FROM_L  2

#define max_gap_threshold  50
#define bandwidth  ALN_BANDWIDTH
#define half_bandwidth  ALN_BANDWIDTH/2

#define min_average_log_emission -5.0f
#define epsilon 1e-10f

/************************* RNA Event Reversal Kernel *************************/

__global__ void reverse_events_kernel(event_t* event_table,
                                       int32_t* n_events,
                                       ptr_t* event_ptr,
                                       int32_t n_bam_rec)
{
    // RNA is sequenced 3'->5' but alignment expects 5'->3' order
    // This kernel reverses events in-place for all reads (RNA-only workflow)
    int i = blockIdx.x * blockDim.x + threadIdx.x;  // One thread per read

    if (i < n_bam_rec) {
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events[i];

        // Reverse events in-place (swap first with last, second with second-to-last, etc.)
        for (int j = 0; j < n_event / 2; ++j) {
            event_t tmp = events[j];
            events[j] = events[n_event - 1 - j];
            events[n_event - 1 - j] = tmp;
        }
    }
}

/************************* Pre-Alignment Kernel *************************/

__global__ void align_kernel_pre_2d(char* read,
    int32_t* read_len, ptr_t* read_ptr,
    int32_t* n_events, ptr_t* event_ptr, model_t* models, uint32_t kmer_size,
    int32_t n_bam_rec, model_t* model_kmer_caches, float *bands1, uint8_t *trace1,
    EventKmerPair* band_lower_left1)
{
    // 1D block indexing: blockIdx.x = read index, threadIdx.x = bandwidth offset
    int i = blockIdx.x;
    int tid = threadIdx.x;

    if (i < n_bam_rec) {
        char* sequence = &read[read_ptr[i]];
        int32_t sequence_len = read_len[i];
        model_t* model_kmer_cache = &model_kmer_caches[read_ptr[i]];
        float *bands = &bands1[(read_ptr[i] + event_ptr[i]) * ALN_BANDWIDTH];
        uint8_t *trace = &trace1[(read_ptr[i] + event_ptr[i]) * ALN_BANDWIDTH];
        EventKmerPair* band_lower_left = &band_lower_left1[read_ptr[i] + event_ptr[i]];

        int32_t n_kmers = sequence_len - kmer_size + 1;
        float lp_trim = logf(0.01f);

        // Precompute k-mer ranks
        if (tid == 0) {
            for (int32_t k = 0; k < n_kmers; ++k) {
                char* substring = &sequence[k];
                uint32_t kmer_rank = get_kmer_rank(substring, kmer_size);
                model_kmer_cache[k] = models[kmer_rank];
            }
        }

        // Initialize bands
        if (tid < bandwidth) {
            for (int32_t k = 0; k < 3; ++k) {
                BAND_ARRAY(k, tid) = -INFINITY;
            }
        }

        if (tid == 0) {
            // Initialize range of first two bands
            band_lower_left[0].event_idx = half_bandwidth - 1;
            band_lower_left[0].kmer_idx = -1 - half_bandwidth;
            band_lower_left[1] = move_down(band_lower_left[0]);

            int start_cell_offset = band_kmer_to_offset(0, -1);
            assert(is_offset_valid(start_cell_offset));
            BAND_ARRAY(0, start_cell_offset) = 0.0f;

            // Band 1: first event is trimmed
            int first_trim_offset = band_event_to_offset(1, 0);
            assert(kmer_at_offset(1, first_trim_offset) == -1);
            assert(is_offset_valid(first_trim_offset));
            BAND_ARRAY(1, first_trim_offset) = lp_trim;
            TRACE_ARRAY(1, first_trim_offset) = FROM_U;
        }
    }
}

/************************* Core Alignment Kernel *************************/

#define band_event_to_offset_shm(bi, ei) band_lower_left_shm[bi].event_idx - (ei)
#define band_kmer_to_offset_shm(bi, ki) (ki) - band_lower_left_shm[bi].kmer_idx
#define event_at_offset_shm(bi, offset) band_lower_left_shm[(bi)].event_idx - (offset)
#define kmer_at_offset_shm(bi, offset) band_lower_left_shm[(bi)].kmer_idx + (offset)
#define BAND_ARRAY_SHM(r, c) ( bands_shm[(r)][(c)] )

__global__ void align_kernel_core_2d_shm(int32_t* read_len, ptr_t* read_ptr,
    event_t* event_table, int32_t* n_events1, ptr_t* event_ptr,
    scalings_t* scalings, int32_t n_bam_rec, model_t* model_kmer_caches, uint32_t kmer_size,
    float *band, uint8_t *traces, EventKmerPair* band_lower_lefts)
{
    // 1D block indexing: blockIdx.x = read index, threadIdx.x = bandwidth offset
    int i = blockIdx.x;
    int offset = threadIdx.x;

    __shared__ float  bands_shm[3][ALN_BANDWIDTH];
    __shared__ EventKmerPair  band_lower_left_shm[3];

    if (i < n_bam_rec && offset < ALN_BANDWIDTH) {
        int32_t sequence_len = read_len[i];
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events1[i];
        scalings_t scaling = scalings[i];
        model_t* model_kmer_cache = &model_kmer_caches[read_ptr[i]];
        float *bands = &band[(read_ptr[i] + event_ptr[i]) * ALN_BANDWIDTH];
        uint8_t *trace = &traces[(read_ptr[i] + event_ptr[i]) * ALN_BANDWIDTH];
        EventKmerPair* band_lower_left = &band_lower_lefts[read_ptr[i] + event_ptr[i]];

        int32_t n_events = n_event;
        int32_t n_kmers = sequence_len - kmer_size + 1;

        // Transition penalties
        float events_per_kmer = (float)n_events / n_kmers;
        float p_stay = 1 - (1 / (events_per_kmer + 1));

        float lp_skip = logf(epsilon);
        float lp_stay = logf(p_stay);
        float lp_step = logf(1.0f - expf(lp_skip) - expf(lp_stay));
        float lp_trim = logf(0.01f);

        int32_t n_rows = n_events + 1;
        int32_t n_cols = n_kmers + 1;
        int32_t n_bands = n_rows + n_cols;

        // Rotate shared memory buffers
        BAND_ARRAY_SHM(0, offset) = BAND_ARRAY(2, offset);
        BAND_ARRAY_SHM(1, offset) = BAND_ARRAY(1, offset);
        BAND_ARRAY_SHM(2, offset) = BAND_ARRAY(0, offset);

        band_lower_left_shm[0] = band_lower_left[2];
        band_lower_left_shm[1] = band_lower_left[1];
        band_lower_left_shm[2] = band_lower_left[0];

        __syncthreads();

        // Fill in remaining bands
        for (int32_t band_idx = 2; band_idx < n_bands; ++band_idx) {
            if (offset == 0) {
                // Determine placement of this band (Suzuki's adaptive algorithm)
                float ll = BAND_ARRAY_SHM(1, 0);
                float ur = BAND_ARRAY_SHM(1, (bandwidth - 1));
                bool ll_ob = ll == -INFINITY;
                bool ur_ob = ur == -INFINITY;

                bool right = false;
                if (ll_ob && ur_ob) {
                    right = band_idx % 2 == 1;
                } else {
                    right = ll < ur; // Suzuki's rule
                }

                if (right) {
                    band_lower_left[band_idx] = band_lower_left_shm[0] =
                        move_right(band_lower_left_shm[1]);
                } else {
                    band_lower_left[band_idx] = band_lower_left_shm[0] =
                        move_down(band_lower_left_shm[1]);
                }

                // Fill trim state
                int trim_offset = band_kmer_to_offset_shm(0, -1);
                if (is_offset_valid(trim_offset)) {
                    int32_t event_idx = event_at_offset_shm(0, trim_offset);
                    if (event_idx >= 0 && event_idx < n_events) {
                        BAND_ARRAY_SHM(0, trim_offset) = lp_trim * (event_idx + 1);
                        TRACE_ARRAY(band_idx, trim_offset) = FROM_U;
                    } else {
                        BAND_ARRAY_SHM(0, trim_offset) = -INFINITY;
                    }
                }
            }
            __syncthreads();

            // Get valid offset range
            int kmer_min_offset = band_kmer_to_offset_shm(0, 0);
            int kmer_max_offset = band_kmer_to_offset_shm(0, n_kmers);
            int event_min_offset = band_event_to_offset_shm(0, n_events - 1);
            int event_max_offset = band_event_to_offset_shm(0, -1);

            int min_offset = MAX(kmer_min_offset, event_min_offset);
            min_offset = MAX(min_offset, 0);

            int max_offset = MIN(kmer_max_offset, event_max_offset);
            max_offset = MIN(max_offset, bandwidth);

            __syncthreads();

            if (offset >= min_offset && offset < max_offset) {
                int event_idx = event_at_offset_shm(0, offset);
                int kmer_idx = kmer_at_offset_shm(0, offset);

                int offset_up = band_event_to_offset_shm(1, event_idx - 1);
                int offset_left = band_kmer_to_offset_shm(1, kmer_idx - 1);
                int offset_diag = band_kmer_to_offset_shm(2, kmer_idx - 1);

                float up = is_offset_valid(offset_up) ? BAND_ARRAY_SHM(1, offset_up) : -INFINITY;
                float left = is_offset_valid(offset_left) ? BAND_ARRAY_SHM(1, offset_left) : -INFINITY;
                float diag = is_offset_valid(offset_diag) ? BAND_ARRAY_SHM(2, offset_diag) : -INFINITY;

                float lp_emission = log_probability_match_r9(scaling, model_kmer_cache, events, event_idx, kmer_idx);

                float score_d = diag + lp_step + lp_emission;
                float score_u = up + lp_stay + lp_emission;
                float score_l = left + lp_skip;

                float max_score = score_d;
                uint8_t from = FROM_D;

                max_score = score_u > max_score ? score_u : max_score;
                from = max_score == score_u ? FROM_U : from;
                max_score = score_l > max_score ? score_l : max_score;
                from = max_score == score_l ? FROM_L : from;

                BAND_ARRAY_SHM(0, offset) = max_score;
                TRACE_ARRAY(band_idx, offset) = from;
            }

            __syncthreads();
            BAND_ARRAY(band_idx, offset) = BAND_ARRAY_SHM(0, offset);

            BAND_ARRAY_SHM(2, offset) = BAND_ARRAY_SHM(1, offset);
            BAND_ARRAY_SHM(1, offset) = BAND_ARRAY_SHM(0, offset);
            BAND_ARRAY_SHM(0, offset) = -INFINITY;

            if (offset == 0) {
                band_lower_left_shm[2] = band_lower_left_shm[1];
                band_lower_left_shm[1] = band_lower_left_shm[0];
            }

            __syncthreads();
        }
    }
}

/************************* Post-Alignment Kernel *************************/

__global__ void align_kernel_post(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs,
    int32_t* read_len, ptr_t* read_ptr,
    event_t* event_table, int32_t* n_events, ptr_t* event_ptr,
    scalings_t* scalings, int32_t n_bam_rec, model_t* model_kmer_caches, uint32_t kmer_size,
    float *bands1, uint8_t *trace1, EventKmerPair* band_lower_left1)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < n_bam_rec) {
        AlignedPair* out_2 = &event_align_pairs[event_ptr[i] * 2];
        int32_t sequence_len = read_len[i];
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events[i];
        scalings_t scaling = scalings[i];
        model_t* model_kmer_cache = &model_kmer_caches[read_ptr[i]];
        float *bands = &bands1[(read_ptr[i] + event_ptr[i]) * ALN_BANDWIDTH];
        uint8_t *trace = &trace1[(read_ptr[i] + event_ptr[i]) * ALN_BANDWIDTH];
        EventKmerPair* band_lower_left = &band_lower_left1[read_ptr[i] + event_ptr[i]];

        int32_t n_kmers = sequence_len - kmer_size + 1;
        int32_t n_rows = n_event + 1;
        int32_t n_cols = n_kmers + 1;
        int32_t n_bands = n_rows + n_cols;

        float lp_trim = logf(0.01f);

        // Backtrack to compute alignment
        double sum_emission = 0;
        double n_aligned_events = 0;

        int outIndex = 0;
        float max_score = -INFINITY;
        int curr_event_idx = 0;
        int curr_kmer_idx = n_kmers - 1;

        // Find best score between an event and the last k-mer
        for (int32_t event_idx = 0; event_idx < n_event; ++event_idx) {
            int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
            assert(band_idx < n_bands);
            int offset = band_event_to_offset(band_idx, event_idx);
            if (is_offset_valid(offset)) {
                float s = BAND_ARRAY(band_idx, offset) + (n_event - event_idx) * lp_trim;
                if (s > max_score) {
                    max_score = s;
                    curr_event_idx = event_idx;
                }
            }
        }

        int curr_gap = 0;
        int max_gap = 0;

        while (curr_kmer_idx >= 0 && curr_event_idx >= 0) {
            // Emit alignment pair
            assert(outIndex < n_event * 2);
            out_2[outIndex].ref_pos = curr_kmer_idx;
            out_2[outIndex].read_pos = curr_event_idx;
            outIndex++;

            // Calculate emission for QC
            float unscaledLevel = events[curr_event_idx].mean;
            float scaledLevel = unscaledLevel;
            model_t model = model_kmer_cache[curr_kmer_idx];
            float gp_mean = scaling.scale * model.level_mean + scaling.shift;
            float gp_stdv = model.level_stdv;

            #ifdef CACHED_LOG
                float gp_log_stdv = model.level_log_stdv;
            #else
                float gp_log_stdv = logf(gp_stdv);
            #endif

            float a = (scaledLevel - gp_mean) / gp_stdv;
            float tempLogProb = log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);

            sum_emission += tempLogProb;
            n_aligned_events++;

            int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
            int offset = band_event_to_offset(band_idx, curr_event_idx);
            assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

            uint8_t from = TRACE_ARRAY(band_idx, offset);
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

        // Note: Reversal is done BEFORE alignment by reverse_events_kernel
        // (RNA workflow reverses events first, then aligns)

        bool spanned = out_2[0].ref_pos == 0 && out_2[outIndex - 1].ref_pos == int(n_kmers - 1);

        // QC results
        double avg_log_emission = sum_emission / n_aligned_events;

        // Check alignment quality
        if (avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
            outIndex = 0;  // Alignment failed
        }

        n_event_align_pairs[i] = outIndex;
    }
}
