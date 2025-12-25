/* @file cuda_framework.cu
**
** GPU framework for event alignment
** Based on f5c CUDA implementation, adapted for pyfin
**
** This file contains:
** - init_cuda_eventalign(): Allocate GPU memory for alignment
** - align_cuda(): Main GPU alignment entry point
** - free_cuda_eventalign(): Clean up GPU memory
** @@
******************************************************************************/

#define HAVE_CUDA 1

#include "align_cuda.h"
#include "common.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

// CUDA error checking macro
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

// Block sizes for CUDA kernels
#define BLOCK_LEN_READS 32
#define BLOCK_LEN_BANDWIDTH 128

/************************* GPU Memory Allocation *************************/

/**
 * Initialize CUDA memory for event alignment
 *
 * @param cuda_data Pointer to cuda_data_t structure to populate
 * @param n_bam_rec Maximum number of reads to process in batch
 * @param model Pore model data
 * @return 0 on success, -1 on failure
 */
int init_cuda_eventalign(cuda_data_t** cuda_data_ptr, int32_t n_bam_rec, model_t* model) {

    cuda_data_t* cuda_data = (cuda_data_t*)malloc(sizeof(cuda_data_t));
    if (!cuda_data) {
        fprintf(stderr, "[CUDA ERROR] Failed to allocate cuda_data_t structure\n");
        return -1;
    }

    // Initialize to NULL
    cuda_data->read_ptr = NULL;
    cuda_data->read_len = NULL;
    cuda_data->read = NULL;
    cuda_data->n_events = NULL;
    cuda_data->event_ptr = NULL;
    cuda_data->event_table = NULL;
    cuda_data->scalings = NULL;
    cuda_data->model = NULL;
    cuda_data->model_kmer_cache = NULL;
    cuda_data->n_event_align_pairs = NULL;
    cuda_data->event_align_pairs = NULL;
    cuda_data->bands = NULL;
    cuda_data->trace = NULL;
    cuda_data->band_lower_left = NULL;

    cuda_data->read_ptr_host = (ptr_t*)malloc(sizeof(ptr_t) * n_bam_rec);
    cuda_data->n_events_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    cuda_data->event_ptr_host = (ptr_t*)malloc(sizeof(ptr_t) * n_bam_rec);
    cuda_data->read_len_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    cuda_data->scalings_host = (scalings_t*)malloc(sizeof(scalings_t) * n_bam_rec);
    cuda_data->n_event_align_pairs_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);

    if (!cuda_data->read_ptr_host || !cuda_data->n_events_host ||
        !cuda_data->event_ptr_host || !cuda_data->read_len_host ||
        !cuda_data->scalings_host || !cuda_data->n_event_align_pairs_host) {
        fprintf(stderr, "[CUDA ERROR] Failed to allocate host arrays\n");
        free_cuda_eventalign(cuda_data);
        return -1;
    }

    // Allocate GPU arrays
    cudaMalloc((void**)&(cuda_data->read_ptr), n_bam_rec * sizeof(ptr_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->read_len), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->n_events), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->event_ptr), n_bam_rec * sizeof(ptr_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->scalings), n_bam_rec * sizeof(scalings_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->model), MAX_NUM_KMER * sizeof(model_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->n_event_align_pairs), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();

    // Copy model to GPU
    cudaMemcpy(cuda_data->model, model, MAX_NUM_KMER * sizeof(model_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();

    // Allocate dynamic arrays (will be sized per-batch)
    // These are allocated dynamically in align_cuda() based on actual data size

    *cuda_data_ptr = cuda_data;
    return 0;
}

/************************* GPU Memory Free *************************/

/**
 * Free CUDA memory for event alignment
 *
 * @param cuda_data Pointer to cuda_data_t structure to free
 */
void free_cuda_eventalign(cuda_data_t* cuda_data) {
    if (!cuda_data) return;

    free(cuda_data->event_ptr_host);
    free(cuda_data->n_events_host);
    free(cuda_data->read_ptr_host);
    free(cuda_data->read_len_host);
    free(cuda_data->scalings_host);
    free(cuda_data->n_event_align_pairs_host);

    if (cuda_data->read_ptr) cudaFree(cuda_data->read_ptr);
    if (cuda_data->read_len) cudaFree(cuda_data->read_len);
    if (cuda_data->n_events) cudaFree(cuda_data->n_events);
    if (cuda_data->event_ptr) cudaFree(cuda_data->event_ptr);
    if (cuda_data->model) cudaFree(cuda_data->model);
    if (cuda_data->scalings) cudaFree(cuda_data->scalings);
    if (cuda_data->n_event_align_pairs) cudaFree(cuda_data->n_event_align_pairs);

    if (cuda_data->read) cudaFree(cuda_data->read);
    if (cuda_data->event_table) cudaFree(cuda_data->event_table);
    if (cuda_data->model_kmer_cache) cudaFree(cuda_data->model_kmer_cache);
    if (cuda_data->event_align_pairs) cudaFree(cuda_data->event_align_pairs);
    if (cuda_data->bands) cudaFree(cuda_data->bands);
    if (cuda_data->trace) cudaFree(cuda_data->trace);
    if (cuda_data->band_lower_left) cudaFree(cuda_data->band_lower_left);

    free(cuda_data);
}

/************************* Main CUDA Alignment Function *************************/

/**
 * Run event alignment on GPU
 *
 * This is the main entry point for GPU-accelerated event alignment.
 * It handles data preparation, kernel launching, and result retrieval.
 *
 * @param cuda_data CUDA data structure with GPU memory
 * @param reads Array of read sequences
 * @param read_lens Array of read lengths
 * @param read_ptrs Array of read pointers (for flattened reads)
 * @param n_reads Number of reads
 * @param event_tables Array of event tables
 * @param n_events Array of event counts per read
 * @param event_ptrs Array of event pointers (for flattened events)
 * @param scalings Array of scaling parameters
 * @param kmer_size K-mer size for alignment
 * @param n_event_align_pairs Output: number of aligned pairs per read
 * @param event_align_pairs Output: aligned event-reference pairs
 * @return 0 on success, -1 on failure
 */
int align_cuda(cuda_data_t* cuda_data,
               char** reads, int32_t* read_lens, ptr_t* read_ptrs,
               int32_t n_reads,
               event_t** event_tables, int32_t* n_events, ptr_t* event_ptrs,
               scalings_t* scalings,
               uint32_t kmer_size,
               int32_t* n_event_align_pairs,
               AlignedPair** event_align_pairs)
{
    if (!cuda_data) {
        fprintf(stderr, "[CUDA ERROR] cuda_data is NULL\n");
        return -1;
    }

    if (n_reads == 0) {
        return 0;  // Nothing to do
    }

    // Calculate total sizes
    int64_t sum_read_len = 0;
    int64_t sum_n_events = 0;
    for (int32_t i = 0; i < n_reads; i++) {
        sum_read_len += (read_lens[i] + 1);  // +1 for null terminator
        sum_n_events += n_events[i];
    }

    // Flatten reads and events on host
    char* read_host = (char*)malloc(sizeof(char) * sum_read_len);
    event_t* event_table_host = (event_t*)malloc(sizeof(event_t) * sum_n_events);
    AlignedPair* event_align_pairs_host = (AlignedPair*)malloc(2 * sum_n_events * sizeof(AlignedPair));

    if (!read_host || !event_table_host || !event_align_pairs_host) {
        fprintf(stderr, "[CUDA ERROR] Failed to allocate host arrays\n");
        free(read_host);
        free(event_table_host);
        free(event_align_pairs_host);
        return -1;
    }

    // Build flattened arrays
    int64_t read_offset = 0;
    int64_t event_offset = 0;
    for (int32_t i = 0; i < n_reads; i++) {
        cuda_data->read_ptr_host[i] = read_offset;
        strcpy(&read_host[read_offset], reads[i]);
        cuda_data->read_len_host[i] = read_lens[i];
        cuda_data->scalings_host[i] = scalings[i];
        read_offset += (read_lens[i] + 1);

        cuda_data->n_events_host[i] = n_events[i];
        cuda_data->event_ptr_host[i] = event_offset;
        memcpy(&event_table_host[event_offset], event_tables[i], sizeof(event_t) * n_events[i]);
        event_offset += n_events[i];
    }

    // Allocate dynamic GPU arrays
    cudaMalloc((void**)&(cuda_data->read), sum_read_len * sizeof(char));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->event_table), sum_n_events * sizeof(event_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->model_kmer_cache), sum_read_len * sizeof(model_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->event_align_pairs), 2 * sum_n_events * sizeof(AlignedPair));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->bands), (sum_n_events + sum_read_len) * ALN_BANDWIDTH * sizeof(float));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->trace), (sum_n_events + sum_read_len) * ALN_BANDWIDTH * sizeof(uint8_t));
    CUDA_CHK();
    cudaMalloc((void**)&(cuda_data->band_lower_left), (sum_n_events + sum_read_len) * sizeof(EventKmerPair));
    CUDA_CHK();

    // Initialize trace to 0
    cudaMemset(cuda_data->trace, 0, sizeof(uint8_t) * (sum_n_events + sum_read_len) * ALN_BANDWIDTH);
    CUDA_CHK();

    // Copy data to GPU
    cudaMemcpy(cuda_data->read_ptr, cuda_data->read_ptr_host, n_reads * sizeof(ptr_t), cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_data->read, read_host, sum_read_len * sizeof(char), cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_data->read_len, cuda_data->read_len_host, n_reads * sizeof(int32_t), cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_data->n_events, cuda_data->n_events_host, n_reads * sizeof(int32_t), cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_data->event_ptr, cuda_data->event_ptr_host, n_reads * sizeof(ptr_t), cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_data->event_table, event_table_host, sizeof(event_t) * sum_n_events, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_data->scalings, cuda_data->scalings_host, sizeof(scalings_t) * n_reads, cudaMemcpyHostToDevice);
    CUDA_CHK();

    // Reverse events for RNA (RNA is sequenced 3'->5' but alignment expects 5'->3')
    dim3 grid_rev((n_reads + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 block_rev(1, BLOCK_LEN_READS);
    reverse_events_kernel<<<grid_rev, block_rev>>>(cuda_data->event_table, cuda_data->n_events,
                                                   cuda_data->event_ptr, n_reads);
    cudaDeviceSynchronize();
    CUDA_CHK();

    // Pre-alignment kernel: k-mer rank precomputation and band initialization
    assert(BLOCK_LEN_BANDWIDTH >= ALN_BANDWIDTH);
    dim3 grid_pre(1, (n_reads + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 block_pre(BLOCK_LEN_BANDWIDTH, BLOCK_LEN_READS);
    align_kernel_pre_2d<<<grid_pre, block_pre>>>(cuda_data->read,
        cuda_data->read_len, cuda_data->read_ptr,
        cuda_data->n_events, cuda_data->event_ptr, cuda_data->model, kmer_size, n_reads,
        cuda_data->model_kmer_cache, cuda_data->bands, cuda_data->trace, cuda_data->band_lower_left);
    cudaDeviceSynchronize();
    CUDA_CHK();

    // Core alignment kernel: adaptive banded event alignment
    dim3 grid_core(1, (n_reads + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 block_core(BLOCK_LEN_BANDWIDTH, BLOCK_LEN_READS);
    align_kernel_core_2d_shm<<<grid_core, block_core>>>(cuda_data->read_len, cuda_data->read_ptr,
        cuda_data->event_table, cuda_data->n_events, cuda_data->event_ptr,
        cuda_data->scalings, n_reads, cuda_data->model_kmer_cache, kmer_size,
        cuda_data->bands, cuda_data->trace, cuda_data->band_lower_left);
    cudaDeviceSynchronize();
    CUDA_CHK();

    // Post-alignment kernel: backtracking and output generation
    int32_t BLOCK_LEN = 128;
    dim3 grid_post((n_reads + BLOCK_LEN - 1) / BLOCK_LEN);
    dim3 block_post(BLOCK_LEN);
    align_kernel_post<<<grid_post, block_post>>>(cuda_data->event_align_pairs,
        cuda_data->n_event_align_pairs, cuda_data->read_len, cuda_data->read_ptr,
        cuda_data->event_table, cuda_data->n_events, cuda_data->event_ptr,
        cuda_data->scalings, n_reads, cuda_data->model_kmer_cache, kmer_size,
        cuda_data->bands, cuda_data->trace, cuda_data->band_lower_left);
    cudaDeviceSynchronize();
    CUDA_CHK();

    // Copy results back to host
    cudaMemcpy(cuda_data->n_event_align_pairs_host, cuda_data->n_event_align_pairs,
               n_reads * sizeof(int32_t), cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(event_align_pairs_host, cuda_data->event_align_pairs,
               2 * sum_n_events * sizeof(AlignedPair), cudaMemcpyDeviceToHost);
    CUDA_CHK();

    // Free dynamic GPU arrays
    cudaFree(cuda_data->read);
    cudaFree(cuda_data->event_table);
    cudaFree(cuda_data->model_kmer_cache);
    cudaFree(cuda_data->event_align_pairs);
    cudaFree(cuda_data->bands);
    cudaFree(cuda_data->trace);
    cudaFree(cuda_data->band_lower_left);

    cuda_data->read = NULL;
    cuda_data->event_table = NULL;
    cuda_data->model_kmer_cache = NULL;
    cuda_data->event_align_pairs = NULL;
    cuda_data->bands = NULL;
    cuda_data->trace = NULL;
    cuda_data->band_lower_left = NULL;

    // Copy results to output arrays
    for (int32_t i = 0; i < n_reads; i++) {
        n_event_align_pairs[i] = cuda_data->n_event_align_pairs_host[i];
        ptr_t event_idx = cuda_data->event_ptr_host[i];
        memcpy(event_align_pairs[i], &event_align_pairs_host[event_idx * 2],
               sizeof(AlignedPair) * n_event_align_pairs[i]);
    }

    free(read_host);
    free(event_table_host);
    free(event_align_pairs_host);

    return 0;
}
