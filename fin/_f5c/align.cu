/*
 * GPU (CUDA) implementation of event-to-sequence alignment with soft-clipping
 * Uses HMM with Viterbi algorithm and flanking states, parallelized for GPU
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "align_common.h"

// CUDA error checking macro
#define CUDA_CHECK(call)                                                     \
    do                                                                       \
    {                                                                        \
        cudaError_t err = call;                                              \
        if (err != cudaSuccess)                                              \
        {                                                                    \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(err));                                \
            return 0;                                                        \
        }                                                                    \
    } while (0)

// Device functions (run on GPU)
__device__ static inline float d_log_normal_pdf(float x, float mean, float stdv, float log_stdv)
{
    float log_inv_sqrt_2pi = -0.918938f;
    float a = (x - mean) / stdv;
    return log_inv_sqrt_2pi - log_stdv + (-0.5f * a * a);
}

__device__ static inline float d_get_scaled_level(float raw_level, simple_scalings_t scaling)
{
    return (raw_level - scaling.shift) / scaling.scale;
}

__device__ static inline uint32_t d_get_rank(char base)
{
    if (base == 'A')
        return 0;
    else if (base == 'C')
        return 1;
    else if (base == 'G')
        return 2;
    else if (base == 'T' || base == 'U')
        return 3;
    else
        return 0;
}

__device__ static inline uint32_t d_get_kmer_rank(const char *str, uint32_t k)
{
    uint32_t r = 0;
    for (uint32_t i = 0; i < k; ++i)
    {
        r += d_get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

// Kernel: Calculate pre-flanking probabilities
__global__ void kernel_make_pre_flanking(float *pre_flank, int32_t num_events)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i == 0)
    {
        pre_flank[0] = logf(1.0f - TRANS_START_TO_CLIP);
        if (num_events > 0)
        {
            pre_flank[1] = logf(TRANS_START_TO_CLIP) + (-3.0f) + logf(1.0f - TRANS_CLIP_SELF);
        }
    }

    __syncthreads();

    // Each thread handles one position
    if (i >= 2 && i <= num_events)
    {
        pre_flank[i] = logf(TRANS_CLIP_SELF) + (-3.0f) + pre_flank[i - 1];
    }
}

// Kernel: Calculate post-flanking probabilities
__global__ void kernel_make_post_flanking(float *post_flank, int32_t num_events)
{
    if (num_events <= 0)
        return;

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i == 0)
    {
        post_flank[num_events - 1] = logf(1.0f - TRANS_START_TO_CLIP);
        if (num_events > 1)
        {
            post_flank[num_events - 2] = logf(TRANS_START_TO_CLIP) + (-3.0f) + logf(1.0f - TRANS_CLIP_SELF);
        }
    }

    __syncthreads();

    // Process in reverse order (handled sequentially for dependencies)
    for (int idx = num_events - 3; idx >= 0; --idx)
    {
        if (i == 0)
        {
            post_flank[idx] = logf(TRANS_CLIP_SELF) + (-3.0f) + post_flank[idx + 1];
        }
    }
}

// Kernel: Initialize first row of DP table
__global__ void kernel_init_dp_row(
    float *dp,
    int *traceback,
    const char *sequence,
    event_t *events,
    simple_model_t *model,
    simple_scalings_t scaling,
    float *pre_flank,
    int32_t n_kmers,
    uint32_t kmer_size)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k < n_kmers)
    {
        uint32_t rank = d_get_kmer_rank(&sequence[k], kmer_size);
        float scaled_level = d_get_scaled_level(events[0].mean, scaling);
        float emission = d_log_normal_pdf(scaled_level, model[rank].level_mean,
                                          model[rank].level_stdv, model[rank].level_log_stdv);

        dp[k] = pre_flank[0] + emission;
        traceback[k] = -1; // Start marker
    }
}

// Kernel: Fill DP table (one event per block, kmers parallelized)
__global__ void kernel_fill_dp(
    float *dp,
    int *traceback,
    const char *sequence,
    event_t *events,
    simple_model_t *model,
    simple_scalings_t scaling,
    float *pre_flank,
    int32_t event_idx,
    int32_t n_kmers,
    int32_t n_events,
    uint32_t kmer_size)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k >= n_kmers || event_idx >= n_events)
        return;

    int curr_offset = event_idx * n_kmers;
    int prev_offset = (event_idx - 1) * n_kmers;

    uint32_t rank = d_get_kmer_rank(&sequence[k], kmer_size);
    float scaled_level = d_get_scaled_level(events[event_idx].mean, scaling);
    float emission = d_log_normal_pdf(scaled_level, model[rank].level_mean,
                                      model[rank].level_stdv, model[rank].level_log_stdv);

    float best_score = -INFINITY;
    int best_prev_k = -1;

    // Transition from same kmer (stay)
    float score_stay = dp[prev_offset + k] + emission;
    if (score_stay > best_score)
    {
        best_score = score_stay;
        best_prev_k = k;
    }

    // Transition from previous kmer (move)
    if (k > 0)
    {
        float score_move = dp[prev_offset + k - 1] + emission;
        if (score_move > best_score)
        {
            best_score = score_move;
            best_prev_k = k - 1;
        }
    }

    // Allow pre-flanking for first kmer
    if (k == 0)
    {
        float score_flank = pre_flank[event_idx] + emission;
        if (score_flank > best_score)
        {
            best_score = score_flank;
            best_prev_k = -2; // Flanking marker
        }
    }

    dp[curr_offset + k] = best_score;
    traceback[curr_offset + k] = best_prev_k;
}

// GPU implementation: Enhanced alignment function with soft-clipping support
// Note: This is a host function that launches GPU kernels
extern "C" int32_t align_with_flanking_gpu(
    simple_aligned_pair_t *out,
    const char *sequence,
    int32_t seq_len,
    event_table events,
    simple_model_t *model,
    uint32_t kmer_size,
    simple_scalings_t scaling)
{
    int32_t n_kmers = seq_len - kmer_size + 1;
    int32_t n_events = events.n;

    if (n_kmers <= 0 || n_events <= 0)
        return 0;

    // Allocate device memory
    float *d_pre_flank, *d_post_flank, *d_dp;
    int *d_traceback;
    char *d_sequence;
    event_t *d_events;
    simple_model_t *d_model;

    size_t dp_size = n_events * n_kmers * sizeof(float);
    size_t tb_size = n_events * n_kmers * sizeof(int);
    int n_kmers_model = 1 << (kmer_size * 2);

    CUDA_CHECK(cudaMalloc(&d_pre_flank, (n_events + 1) * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_post_flank, n_events * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_dp, dp_size));
    CUDA_CHECK(cudaMalloc(&d_traceback, tb_size));
    CUDA_CHECK(cudaMalloc(&d_sequence, seq_len));
    CUDA_CHECK(cudaMalloc(&d_events, n_events * sizeof(event_t)));
    CUDA_CHECK(cudaMalloc(&d_model, n_kmers_model * sizeof(simple_model_t)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_sequence, sequence, seq_len, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_events, events.event, n_events * sizeof(event_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_model, model, n_kmers_model * sizeof(simple_model_t), cudaMemcpyHostToDevice));

    // Initialize DP table with -INFINITY
    CUDA_CHECK(cudaMemset(d_dp, 0xFF, dp_size)); // Bit pattern for -INFINITY
    CUDA_CHECK(cudaMemset(d_traceback, 0xFF, tb_size));

    // Calculate flanking probabilities
    int threads = 256;
    int blocks = (n_events + threads - 1) / threads;
    kernel_make_pre_flanking<<<blocks, threads>>>(d_pre_flank, n_events);
    kernel_make_post_flanking<<<blocks, threads>>>(d_post_flank, n_events);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Initialize first row
    blocks = (n_kmers + threads - 1) / threads;
    kernel_init_dp_row<<<blocks, threads>>>(
        d_dp, d_traceback, d_sequence, d_events, d_model,
        scaling, d_pre_flank, n_kmers, kmer_size);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Fill DP table (process events sequentially, kmers in parallel)
    for (int e = 1; e < n_events; ++e)
    {
        kernel_fill_dp<<<blocks, threads>>>(
            d_dp, d_traceback, d_sequence, d_events, d_model,
            scaling, d_pre_flank, e, n_kmers, n_events, kmer_size);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back to host for traceback
    float *h_dp = (float *)malloc(dp_size);
    int *h_traceback = (int *)malloc(tb_size);
    float *h_post_flank = (float *)malloc(n_events * sizeof(float));

    CUDA_CHECK(cudaMemcpy(h_dp, d_dp, dp_size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_traceback, d_traceback, tb_size, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_post_flank, d_post_flank, n_events * sizeof(float), cudaMemcpyDeviceToHost));

    // Find best ending position (on CPU)
    float best_end_score = -INFINITY;
    int best_end_kmer = -1;
    int best_end_event = -1;

    for (int e = 0; e < n_events; ++e)
    {
        for (int k = 0; k < n_kmers; ++k)
        {
            float score = h_dp[e * n_kmers + k] + h_post_flank[e];
            if (score > best_end_score)
            {
                best_end_score = score;
                best_end_kmer = k;
                best_end_event = e;
            }
        }
    }

    // Traceback (on CPU)
    int out_idx = 0;
    int curr_e = best_end_event;
    int curr_k = best_end_kmer;

    int max_path = n_events + n_kmers;
    int *path_e = (int *)malloc(max_path * sizeof(int));
    int *path_k = (int *)malloc(max_path * sizeof(int));
    int path_len = 0;

    while (curr_e >= 0 && curr_k >= 0 && path_len < max_path)
    {
        path_e[path_len] = curr_e;
        path_k[path_len] = curr_k;
        path_len++;

        int prev_k = h_traceback[curr_e * n_kmers + curr_k];
        if (prev_k == -1 || prev_k == -2)
            break;

        curr_e--;
        curr_k = prev_k;
    }

    // Reverse path
    for (int i = path_len - 1; i >= 0 && out_idx < max_path - 1; --i)
    {
        out[out_idx].ref_pos = path_k[i];
        out[out_idx].read_pos = path_e[i];
        out[out_idx].log_probability = h_dp[path_e[i] * n_kmers + path_k[i]];
        out_idx++;
    }

    // Cleanup
    free(path_e);
    free(path_k);
    free(h_dp);
    free(h_traceback);
    free(h_post_flank);

    cudaFree(d_pre_flank);
    cudaFree(d_post_flank);
    cudaFree(d_dp);
    cudaFree(d_traceback);
    cudaFree(d_sequence);
    cudaFree(d_events);
    cudaFree(d_model);

    return out_idx;
}

// CPU fallback for profile_hmm_align (needed by eventalign.c)
// This function is declared in align_common.h and called from eventalign.c
extern "C" int32_t profile_hmm_align(
    event_alignment_t **out_alignment,
    const char *sequence,
    int32_t seq_len,
    event_table events,
    simple_model_t *model,
    uint32_t kmer_size,
    simple_scalings_t scaling,
    float events_per_base)
{
    // For CUDA build, just use GPU implementation
    // Call align_with_flanking_gpu which is already defined above
    alignment_result_t *gpu_results = (alignment_result_t *)malloc(
        (events.n + seq_len) * sizeof(alignment_result_t));

    if (!gpu_results)
        return -1;

    int32_t n_aligned = align_with_flanking_gpu(
        gpu_results, sequence, seq_len, events, model, kmer_size, scaling);

    if (n_aligned <= 0)
    {
        free(gpu_results);
        return n_aligned;
    }

    // Convert alignment_result_t to event_alignment_t format
    *out_alignment = (event_alignment_t *)malloc(n_aligned * sizeof(event_alignment_t));
    if (!*out_alignment)
    {
        free(gpu_results);
        return -1;
    }

    for (int32_t i = 0; i < n_aligned; ++i)
    {
        (*out_alignment)[i].ref_position = gpu_results[i].ref_pos;
        (*out_alignment)[i].event_idx = gpu_results[i].read_pos;
        (*out_alignment)[i].hmm_state = 'M'; // Match state for GPU results

        // Copy event data if within bounds
        if (gpu_results[i].read_pos >= 0 && gpu_results[i].read_pos < (int32_t)events.n)
        {
            (*out_alignment)[i].event = events.event[gpu_results[i].read_pos];
        }

        // Get k-mer and model parameters
        if (gpu_results[i].ref_pos >= 0 && gpu_results[i].ref_pos < seq_len - (int32_t)kmer_size + 1)
        {
            // Copy k-mer from sequence
            for (uint32_t j = 0; j < kmer_size; ++j)
            {
                (*out_alignment)[i].ref_kmer[j] = sequence[gpu_results[i].ref_pos + j];
            }
            (*out_alignment)[i].ref_kmer[kmer_size] = '\0';

            // Get model parameters
            uint32_t kmer_rank = get_kmer_rank((*out_alignment)[i].ref_kmer, kmer_size);
            if (kmer_rank < model->k)
            {
                (*out_alignment)[i].model_mean = model->mean[kmer_rank];
                (*out_alignment)[i].model_stdv = model->stdv[kmer_rank];

                // Apply scaling
                (*out_alignment)[i].scaled_model_mean =
                    (*out_alignment)[i].model_mean * scaling.scale + scaling.shift;
                (*out_alignment)[i].scaled_model_stdv =
                    (*out_alignment)[i].model_stdv * scaling.var;
            }
        }
    }

    free(gpu_results);
    return n_aligned;
}
