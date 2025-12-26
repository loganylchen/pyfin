/* @file align_cuda.h
**
** CUDA kernel function prototypes for event alignment
** Based on f5c CUDA implementation
** @@
******************************************************************************/

#ifndef ALIGN_CUDA_H
#define ALIGN_CUDA_H

#include "common.h"

#ifdef HAVE_CUDA
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

// CUDA framework function prototypes
int init_cuda_eventalign(cuda_data_t** cuda_data, int32_t n_bam_rec, model_t* model);
void free_cuda_eventalign(cuda_data_t* cuda_data);
int align_cuda(cuda_data_t* cuda_data,
               char** reads, int32_t* read_lens, ptr_t* read_ptrs,
               int32_t n_reads,
               event_t** event_tables, int32_t* n_events, ptr_t* event_ptrs,
               scalings_t* scalings,
               uint32_t kmer_size,
               int32_t* n_event_align_pairs,
               AlignedPair** event_align_pairs);

// RNA event reversal kernel (RNA is sequenced 3'->5' but alignment expects 5'->3')
__global__ void reverse_events_kernel(event_t* event_table,
                                       int32_t* n_events,
                                       ptr_t* event_ptr,
                                       int32_t n_bam_rec);

// Pre-alignment kernel: k-mer rank precomputation and band initialization
__global__ void align_kernel_pre_2d(char* read,
    int32_t* read_len, ptr_t* read_ptr,
    int32_t* n_events, ptr_t* event_ptr, model_t* models, uint32_t kmer_size,
    int32_t n_bam_rec, model_t* model_kmer_caches, float *bands1, uint8_t *trace1,
    EventKmerPair* band_lower_left1);

// Core alignment kernel: adaptive banded event alignment (Suzuki algorithm)
__global__ void align_kernel_core_2d_shm(int32_t* read_len, ptr_t* read_ptr,
    event_t* event_table, int32_t* n_events1, ptr_t* event_ptr,
    scalings_t* scalings, int32_t n_bam_rec, model_t* model_kmer_caches, uint32_t kmer_size,
    float *band, uint8_t *traces, EventKmerPair* band_lower_lefts);

// Post-alignment kernel: backtracking and output generation
__global__ void align_kernel_post(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs,
    int32_t* read_len, ptr_t* read_ptr,
    event_t* event_table, int32_t* n_events, ptr_t* event_ptr,
    scalings_t* scalings, int32_t n_bam_rec, model_t* model_kmer_caches, uint32_t kmer_size,
    float *bands1, uint8_t *trace1, EventKmerPair* band_lower_left1);

#ifdef __cplusplus
}
#endif

// CUDA error checking macro
#define CUDA_CHK()                                                             \
    {                                                                           \
        cudaError_t code = cudaGetLastError();                                  \
        if (code != cudaSuccess) {                                              \
            fprintf(stderr, "[CUDA ERROR] %s:%d: %s\n", __FILE__, __LINE__,     \
                    cudaGetErrorString(code));                                 \
            exit(-1);                                                           \
        }                                                                       \
    }

#endif // HAVE_CUDA

#endif // ALIGN_CUDA_H
