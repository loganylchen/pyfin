/*
 * Common definitions and structures for event alignment
 * Shared between CPU (align.c) and GPU (align.cu) implementations
 */

#ifndef ALIGN_COMMON_H
#define ALIGN_COMMON_H

#include <stdint.h>
#include "event_detection_simple.h"

// Soft-clipping transition probabilities (from f5c)
// For RNA: events are reversed before alignment, so:
//   - Pre-flanking clips events at the 5' end of sequence (END of raw signal)
//   - Post-flanking clips events at the 3' end of sequence (START of raw signal)
#define TRANS_START_TO_CLIP 0.5f // Probability of entering clipping state
#define TRANS_CLIP_SELF 0.9f     // Probability of staying in clipping state

// Flags to modify the behaviour of the HMM (from f5c)
enum HMMAlignmentFlags
{
    HAF_ALLOW_PRE_CLIP = 1, // allow events to go unmatched before the aligning region
    HAF_ALLOW_POST_CLIP = 2 // allow events to go unmatched after the aligning region
};

// Alternative parameters for testing asymmetric clipping:
// #define TRANS_START_TO_CLIP_PRE 0.3f   // Less aggressive at 5' end
// #define TRANS_START_TO_CLIP_POST 0.7f  // More aggressive at 3' end

// HMM states for alignment (from f5c R9 profile)
typedef enum
{
    STATE_MATCH = 0,     // Event matches k-mer
    STATE_BAD_EVENT = 1, // Bad event (should be skipped)
    STATE_KMER_SKIP = 2, // K-mer with no event
    NUM_STATES = 3
} HMMState;

// Simplified model structure
typedef struct
{
    float level_mean;
    float level_stdv;
    float level_log_stdv;
} simple_model_t;

// Simplified scaling structure
typedef struct
{
    float scale;
    float shift;
    float var;
    float log_var;
} simple_scalings_t;

// Simplified aligned pair structure (for ABEA algorithm)
typedef struct
{
    int ref_pos;           // kmer index in sequence
    int read_pos;          // event index
    float log_probability; // Alignment probability
} simple_aligned_pair_t;

// Full event alignment structure (from f5c eventalign)
#define MAX_KMER_SIZE 16
typedef struct
{
    int32_t ref_position;               // Reference position (0-based)
    char ref_kmer[MAX_KMER_SIZE + 1];   // Reference k-mer
    int32_t event_idx;                  // Event index (in raw signal order)
    uint64_t signal_start;              // Start sample index in raw signal
    float signal_length;                // Length of event in samples
    char hmm_state;                     // 'M' = match, 'K' = kmer_skip, 'B' = bad_event
    uint8_t strand_idx;                 // 0=template, 1=complement (for DNA)
    char model_kmer[MAX_KMER_SIZE + 1]; // Model k-mer
    float event_mean;                   // Observed event mean
    float event_stdv;                   // Observed event stdv
    float event_duration;               // Event duration
    float model_mean;                   // Expected model mean
    float model_stdv;                   // Expected model stdv
    float scaled_model_mean;            // scale * model_mean + shift
    float scaled_model_stdv;            // model_stdv * var
} event_alignment_t;

// Function prototypes for alignment
#ifdef __cplusplus
extern "C"
{
#endif

    // CPU alignment functions

    // Simple ABEA alignment (fast, returns simple pairs)
    int32_t align_with_flanking_cpu(
        simple_aligned_pair_t **out_alignment,
        const char *sequence,
        int32_t seq_len,
        event_table events,
        simple_model_t *model,
        uint32_t kmer_size,
        simple_scalings_t scaling,
        uint32_t hmm_flags,
        uint32_t e_start);

    // Profile HMM eventalign (detailed, returns event_alignment_t)
    int32_t profile_hmm_align(
        event_alignment_t **out_alignment,
        const char *sequence,
        int32_t seq_len,
        event_table events,
        simple_model_t *model,
        uint32_t kmer_size,
        simple_scalings_t scaling,
        float events_per_base);

#ifdef __cplusplus
}
#endif

// Utility functions (inline for performance)
static inline uint32_t get_rank(char base)
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

static inline uint32_t get_kmer_rank(const char *str, uint32_t k)
{
    uint32_t r = 0;
    for (uint32_t i = 0; i < k; ++i)
    {
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

#endif // ALIGN_COMMON_H
