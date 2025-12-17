/*
 * CPU implementation of event-to-sequence alignment with soft-clipping
 * Full f5c HMM implementation with 3 states (MATCH, BAD_EVENT, KMER_SKIP)
 * Based on f5c/nanopolish eventalign algorithm
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "align_common.h"

// HMM transition types (from f5c)
typedef enum
{
    HMT_FROM_SAME_M = 0, // Match to match (same kmer)
    HMT_FROM_PREV_M,     // Match to match (next kmer)
    HMT_FROM_SAME_B,     // Bad event to match (same kmer)
    HMT_FROM_PREV_B,     // Bad event to match (next kmer)
    HMT_FROM_PREV_K,     // Kmer skip to match
    HMT_FROM_SOFT,       // Soft clip to match (flanking)
    HMT_NUM_MOVEMENT_TYPES
} HMMMovementType;

// HMM update scores structure
typedef struct
{
    float x[HMT_NUM_MOVEMENT_TYPES];
} HMMUpdateScores;

// Block transition probabilities (from f5c)
typedef struct
{
    // Transitions from match state
    float lp_mm_self; // Stay in same kmer
    float lp_mb;      // Match to bad event
    float lp_mk;      // Match to kmer skip
    float lp_mm_next; // Match to next kmer

    // Transitions from bad event state
    float lp_bb;      // Bad event to bad event
    float lp_bk;      // Bad event to kmer skip
    float lp_bm_next; // Bad event to next kmer match
    float lp_bm_self; // Bad event to same kmer match

    // Transitions from kmer skip state
    float lp_kk; // Kmer skip to kmer skip
    float lp_km; // Kmer skip to match
} BlockTransitions;

// Calculate log normal PDF for emission probability (from f5c)
static inline float log_normal_pdf(float x, float gp_mean, float gp_stdv, float gp_log_stdv)
{
    float log_inv_sqrt_2pi = -0.918938f;
    float a = (x - gp_mean) / gp_stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
}

// Calculate emission probability using f5c's scaling model
static inline float log_probability_match(simple_scalings_t scaling,
                                          simple_model_t *models,
                                          event_t *events,
                                          int event_idx,
                                          uint32_t kmer_rank)
{
    assert(kmer_rank < (1u << 10)); // Max for 5-mer

    float unscaled_level = events[event_idx].mean;
    float scaled_level = unscaled_level;

    // f5c scaling: model parameters are scaled to match events
    float gp_mean = scaling.scale * models[kmer_rank].level_mean + scaling.shift;
    float gp_stdv = models[kmer_rank].level_stdv * scaling.var;
    float gp_log_stdv = models[kmer_rank].level_log_stdv + scaling.log_var;

    return log_normal_pdf(scaled_level, gp_mean, gp_stdv, gp_log_stdv);
}

// Calculate transition probabilities for each kmer (from f5c)
static BlockTransitions *calculate_transitions(int32_t n_kmers, double events_per_base)
{
    BlockTransitions *transitions = (BlockTransitions *)malloc(n_kmers * sizeof(BlockTransitions));
    if (!transitions)
        return NULL;

    for (int32_t ki = 0; ki < n_kmers; ++ki)
    {
        // f5c transition model
        float p_stay = 1.0f - (1.0f / events_per_base);
        float p_skip = 0.0025f;
        float p_bad = 0.001f;
        float p_bad_self = p_bad;
        float p_skip_self = 0.3f;

        // Transitions from match state
        float p_mk = p_skip;
        float p_mb = p_bad;
        float p_mm_self = p_stay;
        float p_mm_next = 1.0f - p_mm_self - p_mk - p_mb;

        // Transitions from bad event state
        float p_bb = p_bad_self;
        float p_bk = (1.0f - p_bb) / 3.0f;
        float p_bm_next = (1.0f - p_bb) / 3.0f;
        float p_bm_self = (1.0f - p_bb) / 3.0f;

        // Transitions from kmer skip state
        float p_kk = p_skip_self;
        float p_km = 1.0f - p_kk;

        // Log-transform and store
        BlockTransitions *bt = &transitions[ki];
        bt->lp_mm_self = logf(p_mm_self);
        bt->lp_mb = logf(p_mb);
        bt->lp_mk = logf(p_mk);
        bt->lp_mm_next = logf(p_mm_next);

        bt->lp_bb = logf(p_bb);
        bt->lp_bk = logf(p_bk);
        bt->lp_bm_next = logf(p_bm_next);
        bt->lp_bm_self = logf(p_bm_self);

        bt->lp_kk = logf(p_kk);
        bt->lp_km = logf(p_km);
    }

    return transitions;
}

// Calculate pre-flanking probabilities (from f5c)
static float *make_pre_flanking(int32_t num_events)
{
    float *pre_flank = (float *)calloc(num_events + 1, sizeof(float));
    if (!pre_flank)
        return NULL;

    pre_flank[0] = logf(1.0f - TRANS_START_TO_CLIP);

    if (num_events > 0)
    {
        pre_flank[1] = logf(TRANS_START_TO_CLIP) + (-3.0f) + logf(1.0f - TRANS_CLIP_SELF);

        for (int i = 2; i <= num_events; ++i)
        {
            pre_flank[i] = logf(TRANS_CLIP_SELF) + (-3.0f) + pre_flank[i - 1];
        }
    }

    return pre_flank;
}

// Calculate post-flanking probabilities (from f5c)
// post_flank[i] = log probability of soft-clipping events i..n_events-1
static float *make_post_flanking(int32_t num_events)
{
    // Allocate num_events + 1 elements (rows 0..num_events)
    float *post_flank = (float *)calloc(num_events + 1, sizeof(float));
    if (!post_flank)
        return NULL;

    // post_flank[num_events] = log(1 - TRANS_START_TO_CLIP) = no clipping
    post_flank[num_events] = logf(1.0f - TRANS_START_TO_CLIP);

    if (num_events > 0)
    {
        // post_flank[num_events-1] = soft-clip last 1 event
        post_flank[num_events - 1] = logf(TRANS_START_TO_CLIP) + (-3.0f) + logf(1.0f - TRANS_CLIP_SELF);

        for (int i = num_events - 2; i >= 0; --i)
        {
            post_flank[i] = logf(TRANS_CLIP_SELF) + (-3.0f) + post_flank[i + 1];
        }
    }

    return post_flank;
}

// Update cell with maximum score and traceback (from f5c HMM logic)
static inline float update_cell(float *cell, HMMUpdateScores *scores, float emission,
                                int *best_prev_state, int *best_prev_block,
                                int curr_block, int curr_state)
{
    float max_score = -INFINITY;
    int best_type = -1;

    for (int i = 0; i < HMT_NUM_MOVEMENT_TYPES; ++i)
    {
        if (scores->x[i] > max_score)
        {
            max_score = scores->x[i];
            best_type = i;
        }
    }

    // Determine previous state and block based on movement type
    if (best_type == HMT_FROM_SAME_M)
    {
        *best_prev_state = STATE_MATCH;
        *best_prev_block = curr_block;
    }
    else if (best_type == HMT_FROM_PREV_M)
    {
        *best_prev_state = STATE_MATCH;
        *best_prev_block = curr_block - 1;
    }
    else if (best_type == HMT_FROM_SAME_B)
    {
        *best_prev_state = STATE_BAD_EVENT;
        *best_prev_block = curr_block;
    }
    else if (best_type == HMT_FROM_PREV_B)
    {
        *best_prev_state = STATE_BAD_EVENT;
        *best_prev_block = curr_block - 1;
    }
    else if (best_type == HMT_FROM_PREV_K)
    {
        *best_prev_state = STATE_KMER_SKIP;
        *best_prev_block = curr_block - 1;
    }
    else if (best_type == HMT_FROM_SOFT)
    {
        *best_prev_state = -1; // Start/flanking
        *best_prev_block = -1;
    }

    *cell = max_score + emission;
    return *cell;
}

// CPU implementation: Full f5c HMM with 3 states
int32_t align_with_flanking_cpu(
    simple_aligned_pair_t **out_alignment,
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

    // Estimate events per base for transition probabilities
    double events_per_base = (double)n_events / (double)n_kmers;
    if (events_per_base < 1.0)
        events_per_base = 1.0;

    // Calculate transition probabilities for each kmer
    BlockTransitions *transitions = calculate_transitions(n_kmers, events_per_base);
    if (!transitions)
        return 0;

    // Calculate pre and post flanking probabilities
    float *pre_flank = make_pre_flanking(n_events);
    float *post_flank = make_post_flanking(n_events);

    if (!pre_flank || !post_flank || !transitions)
    {
        free(transitions);
        free(pre_flank);
        free(post_flank);
        return 0;
    }

    // f5c 3-state HMM DP table: [event][kmer * NUM_STATES]
    // States: MATCH(0), BAD_EVENT(1), KMER_SKIP(2)
    int32_t num_blocks = n_kmers + 2; // +2 for start/end states
    int32_t cells_per_row = num_blocks * NUM_STATES;

    float **dp = (float **)malloc((n_events + 1) * sizeof(float *));
    int **traceback_state = (int **)malloc((n_events + 1) * sizeof(int *));
    int **traceback_kmer = (int **)malloc((n_events + 1) * sizeof(int *));

    if (!dp || !traceback_state || !traceback_kmer)
    {
        free(transitions);
        free(pre_flank);
        free(post_flank);
        free(dp);
        free(traceback_state);
        free(traceback_kmer);
        return 0;
    }

    for (int i = 0; i <= n_events; ++i)
    {
        dp[i] = (float *)malloc(cells_per_row * sizeof(float));
        traceback_state[i] = (int *)malloc(cells_per_row * sizeof(int));
        traceback_kmer[i] = (int *)malloc(cells_per_row * sizeof(int));

        if (!dp[i] || !traceback_state[i] || !traceback_kmer[i])
        {
            for (int j = 0; j <= i; ++j)
            {
                free(dp[j]);
                free(traceback_state[j]);
                free(traceback_kmer[j]);
            }
            free(dp);
            free(traceback_state);
            free(traceback_kmer);
            free(transitions);
            free(pre_flank);
            free(post_flank);
            return 0;
        }

        // Initialize all cells
        for (int c = 0; c < cells_per_row; ++c)
        {
            dp[i][c] = -INFINITY;
            traceback_state[i][c] = -1;
            traceback_kmer[i][c] = -1;
        }
    }

    // Initialize: start state transitions to first kmer
    float lp_sm = 0.0f; // Log prob from start to match
    float BAD_EVENT_PENALTY = 0.0f;

    // Row 0 is start state, initialize row 1 (first event)
    for (int32_t block = 1; block < num_blocks - 1; ++block)
    {
        int32_t kmer_idx = block - 1;
        int32_t curr_offset = NUM_STATES * block;
        uint32_t rank = get_kmer_rank(&sequence[kmer_idx], kmer_size);

        float emission = log_probability_match(scaling, model, events.event, 0, rank);

        // Can start from soft clip (pre-flanking)
        if (kmer_idx == 0)
        {
            dp[1][curr_offset + STATE_MATCH] = lp_sm + pre_flank[0] + emission;
            traceback_state[1][curr_offset + STATE_MATCH] = -1;
            traceback_kmer[1][curr_offset + STATE_MATCH] = -1;
        }
    }

    // Fill DP table with f5c's 3-state HMM
    for (int32_t row = 2; row <= n_events; ++row)
    {
        int32_t event_idx = row - 1;

        for (int32_t block = 1; block < num_blocks - 1; ++block)
        {
            int32_t kmer_idx = block - 1;
            BlockTransitions *bt = &transitions[kmer_idx];

            int32_t prev_block = block - 1;
            int32_t prev_offset = NUM_STATES * prev_block;
            int32_t curr_offset = NUM_STATES * block;

            uint32_t rank = get_kmer_rank(&sequence[kmer_idx], kmer_size);
            float lp_emission_m = log_probability_match(scaling, model, events.event, event_idx, rank);
            float lp_emission_b = BAD_EVENT_PENALTY;

            HMMUpdateScores scores;
            int prev_state_trace, prev_block_trace;

            // STATE_MATCH: Event matches this kmer
            scores.x[HMT_FROM_SAME_M] = bt->lp_mm_self + dp[row - 1][curr_offset + STATE_MATCH];
            scores.x[HMT_FROM_PREV_M] = bt->lp_mm_next + dp[row - 1][prev_offset + STATE_MATCH];
            scores.x[HMT_FROM_SAME_B] = bt->lp_bm_self + dp[row - 1][curr_offset + STATE_BAD_EVENT];
            scores.x[HMT_FROM_PREV_B] = bt->lp_bm_next + dp[row - 1][prev_offset + STATE_BAD_EVENT];
            scores.x[HMT_FROM_PREV_K] = bt->lp_km + dp[row - 1][prev_offset + STATE_KMER_SKIP];
            scores.x[HMT_FROM_SOFT] = (kmer_idx == 0) ? lp_sm + pre_flank[row - 1] : -INFINITY;

            update_cell(&dp[row][curr_offset + STATE_MATCH], &scores, lp_emission_m,
                        &prev_state_trace, &prev_block_trace, block, STATE_MATCH);
            traceback_state[row][curr_offset + STATE_MATCH] = prev_state_trace;
            traceback_kmer[row][curr_offset + STATE_MATCH] = prev_block_trace;

            // STATE_BAD_EVENT: This event should be ignored
            scores.x[HMT_FROM_SAME_M] = bt->lp_mb + dp[row - 1][curr_offset + STATE_MATCH];
            scores.x[HMT_FROM_PREV_M] = -INFINITY;
            scores.x[HMT_FROM_SAME_B] = bt->lp_bb + dp[row - 1][curr_offset + STATE_BAD_EVENT];
            scores.x[HMT_FROM_PREV_B] = -INFINITY;
            scores.x[HMT_FROM_PREV_K] = -INFINITY;
            scores.x[HMT_FROM_SOFT] = -INFINITY;

            update_cell(&dp[row][curr_offset + STATE_BAD_EVENT], &scores, lp_emission_b,
                        &prev_state_trace, &prev_block_trace, block, STATE_BAD_EVENT);
            traceback_state[row][curr_offset + STATE_BAD_EVENT] = prev_state_trace;
            traceback_kmer[row][curr_offset + STATE_BAD_EVENT] = prev_block_trace;

            // STATE_KMER_SKIP: No event for this kmer
            scores.x[HMT_FROM_SAME_M] = -INFINITY;
            scores.x[HMT_FROM_PREV_M] = bt->lp_mk + dp[row][prev_offset + STATE_MATCH];
            scores.x[HMT_FROM_SAME_B] = -INFINITY;
            scores.x[HMT_FROM_PREV_B] = bt->lp_bk + dp[row][prev_offset + STATE_BAD_EVENT];
            scores.x[HMT_FROM_PREV_K] = bt->lp_kk + dp[row][prev_offset + STATE_KMER_SKIP];
            scores.x[HMT_FROM_SOFT] = -INFINITY;

            update_cell(&dp[row][curr_offset + STATE_KMER_SKIP], &scores, 0.0f,
                        &prev_state_trace, &prev_block_trace, block, STATE_KMER_SKIP);
            traceback_state[row][curr_offset + STATE_KMER_SKIP] = prev_state_trace;
            traceback_kmer[row][curr_offset + STATE_KMER_SKIP] = prev_block_trace;
        }
    }

    // Find best ending: check post-flanking transitions
    float best_score = -INFINITY;
    int best_block = -1;
    int best_state = -1;

    for (int32_t block = 1; block < num_blocks - 1; ++block)
    {
        int32_t offset = NUM_STATES * block;
        for (int state = 0; state < NUM_STATES; ++state)
        {
            float score = dp[n_events][offset + state];
            if (state != STATE_KMER_SKIP) // Can't end on skip
            {
                score += post_flank[n_events - 1];
            }
            if (score > best_score)
            {
                best_score = score;
                best_block = block;
                best_state = state;
            }
        }
    }

    if (best_block < 0)
    {
        // No valid alignment
        for (int i = 0; i <= n_events; ++i)
        {
            free(dp[i]);
            free(traceback_state[i]);
            free(traceback_kmer[i]);
        }
        free(dp);
        free(traceback_state);
        free(traceback_kmer);
        free(transitions);
        free(pre_flank);
        free(post_flank);
        return 0;
    }

    // Traceback through 3-state HMM
    typedef struct
    {
        int event_idx;
        int kmer_idx;
        int state;
    } PathElement;

    PathElement *path = (PathElement *)malloc(n_events * sizeof(PathElement));
    int path_len = 0;

    int curr_row = n_events;
    int curr_block = best_block;
    int curr_state = best_state;

    while (curr_row > 0)
    {
        int32_t curr_offset = NUM_STATES * curr_block;
        int32_t kmer_idx = curr_block - 1;

        // Record if MATCH or BAD_EVENT (not KMER_SKIP)
        if (curr_state == STATE_MATCH || curr_state == STATE_BAD_EVENT)
        {
            path[path_len].event_idx = curr_row - 1;
            path[path_len].kmer_idx = kmer_idx;
            path[path_len].state = curr_state;
            path_len++;
        }

        // Traceback
        int prev_state = traceback_state[curr_row][curr_offset + curr_state];
        int prev_block = traceback_kmer[curr_row][curr_offset + curr_state];

        if (prev_state < 0 || prev_block < 0)
            break; // Reached start

        // KMER_SKIP doesn't consume event, others do
        if (curr_state == STATE_KMER_SKIP)
        {
            curr_block = prev_block;
            curr_state = prev_state;
            // Don't decrement row
        }
        else
        {
            curr_block = prev_block;
            curr_state = prev_state;
            curr_row--;
        }
    }

    // Reverse path
    for (int i = 0; i < path_len / 2; ++i)
    {
        PathElement temp = path[i];
        path[i] = path[path_len - 1 - i];
        path[path_len - 1 - i] = temp;
    }

    // Generate aligned pairs (only MATCH states)
    simple_aligned_pair_t *alignment = (simple_aligned_pair_t *)malloc(path_len * sizeof(simple_aligned_pair_t));
    int align_len = 0;

    if (!alignment)
    {
        free(path);
        for (int i = 0; i <= n_events; ++i)
        {
            free(dp[i]);
            free(traceback_state[i]);
            free(traceback_kmer[i]);
        }
        free(dp);
        free(traceback_state);
        free(traceback_kmer);
        free(transitions);
        free(pre_flank);
        free(post_flank);
        return 0;
    }

    for (int i = 0; i < path_len; ++i)
    {
        // Only output MATCH states (ignore BAD_EVENT)
        if (path[i].state == STATE_MATCH)
        {
            alignment[align_len].ref_pos = path[i].kmer_idx;
            alignment[align_len].read_pos = path[i].event_idx;
            align_len++;
        }
    }

    // Clean up
    free(path);
    for (int i = 0; i <= n_events; ++i)
    {
        free(dp[i]);
        free(traceback_state[i]);
        free(traceback_kmer[i]);
    }
    free(dp);
    free(traceback_state);
    free(traceback_kmer);
    free(transitions);
    free(pre_flank);
    free(post_flank);

    // Set output
    *out_alignment = alignment;
    return align_len;
}

//-----------------------------------------------------------------------------
// Profile HMM Eventalign (Full f5c implementation)
//-----------------------------------------------------------------------------

// Profile HMM fill function (Viterbi algorithm) - with soft-clipping like f5c
static void profile_hmm_fill_generic_r9(
    float **dp_matrix,
    uint8_t **path_matrix,
    const char *sequence,
    int32_t n_kmers,
    event_table events,
    int32_t n_events,
    simple_model_t *models,
    uint32_t kmer_size,
    simple_scalings_t scaling,
    BlockTransitions *transitions,
    float *pre_flank,
    float *post_flank)
{
    // Initialize all cells to -INFINITY
    for (int32_t ei = 0; ei <= n_events; ++ei)
    {
        for (int32_t ki = 0; ki < n_kmers; ++ki)
        {
            for (int s = 0; s < NUM_STATES; ++s)
            {
                int state_idx = ki * NUM_STATES + s;
                dp_matrix[ei][state_idx] = -INFINITY;
                path_matrix[ei][state_idx] = 0;
            }
        }
    }

    // Fill DP table - f5c style with pre_flank soft-clipping
    // Row 0 is before any events (initialization)
    // Rows 1..n_events correspond to events 0..n_events-1
    for (int32_t row = 1; row <= n_events; ++row)
    {
        int32_t event_idx = row - 1;

        for (int32_t ki = 0; ki < n_kmers; ++ki)
        {
            // Get kmer from sequence
            const char *kmer_ptr = sequence + ki;
            uint32_t kmer_rank = get_kmer_rank(kmer_ptr, kmer_size);

            // Get emission probability
            float lp_emission = log_probability_match(scaling, models, events.event, event_idx, kmer_rank);
            float lp_emission_b = 0.0f; // Bad event penalty

            // Calculate scores for each state
            for (int s = 0; s < NUM_STATES; ++s)
            {
                int curr_state_idx = ki * NUM_STATES + s;
                HMMUpdateScores scores;
                for (int i = 0; i < HMT_NUM_MOVEMENT_TYPES; ++i)
                    scores.x[i] = -INFINITY;

                if (s == STATE_MATCH)
                {
                    // Match state - event aligns to this kmer

                    // From same kmer match (stay) - previous row, same kmer
                    {
                        int prev_idx = ki * NUM_STATES + STATE_MATCH;
                        scores.x[HMT_FROM_SAME_M] = dp_matrix[row - 1][prev_idx] + transitions[ki].lp_mm_self;
                    }

                    // From previous kmer match (advance) - previous row, previous kmer
                    if (ki > 0)
                    {
                        int prev_idx = (ki - 1) * NUM_STATES + STATE_MATCH;
                        scores.x[HMT_FROM_PREV_M] = dp_matrix[row - 1][prev_idx] + transitions[ki - 1].lp_mm_next;
                    }

                    // From same kmer bad event - previous row, same kmer
                    {
                        int prev_idx = ki * NUM_STATES + STATE_BAD_EVENT;
                        scores.x[HMT_FROM_SAME_B] = dp_matrix[row - 1][prev_idx] + transitions[ki].lp_bm_self;
                    }

                    // From previous kmer bad event - previous row, previous kmer
                    if (ki > 0)
                    {
                        int prev_idx = (ki - 1) * NUM_STATES + STATE_BAD_EVENT;
                        scores.x[HMT_FROM_PREV_B] = dp_matrix[row - 1][prev_idx] + transitions[ki - 1].lp_bm_next;
                    }

                    // From previous kmer skip - previous row, previous kmer
                    if (ki > 0)
                    {
                        int prev_idx = (ki - 1) * NUM_STATES + STATE_KMER_SKIP;
                        scores.x[HMT_FROM_PREV_K] = dp_matrix[row - 1][prev_idx] + transitions[ki - 1].lp_km;
                    }

                    // From soft-clip (pre-flank) - only allowed for first kmer (ki == 0)
                    // This allows skipping the first N events with a penalty
                    // Note: emission is added later, so this just provides the start probability
                    if (ki == 0)
                    {
                        // pre_flank[row-1] gives log prob of skipping events 0..(row-2)
                        // and then transitioning from clip to match state
                        scores.x[HMT_FROM_SOFT] = pre_flank[row - 1];
                    }
                }
                else if (s == STATE_BAD_EVENT)
                {
                    // Bad event state - ignore this event

                    // From match to bad event - previous row, same kmer
                    int prev_idx = ki * NUM_STATES + STATE_MATCH;
                    scores.x[HMT_FROM_SAME_M] = dp_matrix[row - 1][prev_idx] + transitions[ki].lp_mb;

                    // From bad event to bad event - previous row, same kmer
                    prev_idx = ki * NUM_STATES + STATE_BAD_EVENT;
                    scores.x[HMT_FROM_SAME_B] = dp_matrix[row - 1][prev_idx] + transitions[ki].lp_bb;
                }
                else if (s == STATE_KMER_SKIP)
                {
                    // Kmer skip state - no event for this kmer
                    // Note: KMER_SKIP doesn't consume events, operates on same row
                    if (ki > 0)
                    {
                        int prev_idx = (ki - 1) * NUM_STATES + STATE_MATCH;
                        scores.x[HMT_FROM_PREV_M] = dp_matrix[row][prev_idx] + transitions[ki - 1].lp_mk;

                        prev_idx = (ki - 1) * NUM_STATES + STATE_BAD_EVENT;
                        scores.x[HMT_FROM_PREV_B] = dp_matrix[row][prev_idx] + transitions[ki - 1].lp_bk;

                        prev_idx = (ki - 1) * NUM_STATES + STATE_KMER_SKIP;
                        scores.x[HMT_FROM_PREV_K] = dp_matrix[row][prev_idx] + transitions[ki - 1].lp_kk;
                    }
                }

                // Find best path and score
                float best_score = -INFINITY;
                uint8_t best_movement = 0;
                for (int m = 0; m < HMT_NUM_MOVEMENT_TYPES; ++m)
                {
                    if (scores.x[m] > best_score)
                    {
                        best_score = scores.x[m];
                        best_movement = m;
                    }
                }

                // Add emission for MATCH and BAD_EVENT states (not KMER_SKIP)
                if (s == STATE_MATCH)
                {
                    dp_matrix[row][curr_state_idx] = best_score + lp_emission;
                }
                else if (s == STATE_BAD_EVENT)
                {
                    dp_matrix[row][curr_state_idx] = best_score + lp_emission_b;
                }
                else
                {
                    dp_matrix[row][curr_state_idx] = best_score; // KMER_SKIP has no emission
                }
                path_matrix[row][curr_state_idx] = best_movement;
            }
        }
    }
}

// Traceback through HMM to get alignment path
// Now uses post_flank to find best ending position (soft-clipping at end)
static int32_t profile_hmm_traceback_with_flanking(
    event_alignment_t **out_alignment,
    uint8_t **path_matrix,
    float **dp_matrix,
    const char *sequence,
    int32_t n_kmers,
    event_table events,
    int32_t n_events,
    simple_model_t *models,
    uint32_t kmer_size,
    simple_scalings_t scaling,
    float *post_flank)
{
    // Find best end position using post_flank
    // Search over all events at the last kmer, considering post-clipping penalty
    float best_score = -INFINITY;
    int best_event_row = n_events;
    int best_state = STATE_MATCH;

    for (int32_t row = 1; row <= n_events; ++row)
    {
        for (int s = 0; s < NUM_STATES; ++s)
        {
            int state_idx = (n_kmers - 1) * NUM_STATES + s;
            // post_flank[row] gives log prob of soft-clipping events row..n_events-1
            float score = dp_matrix[row][state_idx] + post_flank[row];
            if (score > best_score)
            {
                best_score = score;
                best_event_row = row;
                best_state = s;
            }
        }
    }

    // Allocate alignment array (worst case: all events)
    event_alignment_t *alignment = (event_alignment_t *)malloc(n_events * sizeof(event_alignment_t));
    if (!alignment)
        return -1;

    int align_idx = 0;
    int curr_kmer = n_kmers - 1;
    int curr_state = best_state;
    int curr_event = best_event_row - 1; // Convert row to event index

    // Traceback
    while (curr_event >= 0 && curr_kmer >= 0)
    {
        int state_idx = curr_kmer * NUM_STATES + curr_state;
        uint8_t movement = path_matrix[curr_event + 1][state_idx];

        // Record alignment for match states
        if (curr_state == STATE_MATCH && curr_event >= 0)
        {
            event_alignment_t *aln = &alignment[align_idx++];

            // Fill in alignment details
            aln->ref_position = curr_kmer;
            aln->event_idx = curr_event;
            aln->hmm_state = 'M';
            aln->strand_idx = 0;

            // Copy kmers with bounds checking
            int copy_len = (curr_kmer + kmer_size <= strlen(sequence)) ? kmer_size : strlen(sequence) - curr_kmer;
            if (copy_len > 0 && copy_len <= MAX_KMER_SIZE)
            {
                strncpy(aln->ref_kmer, sequence + curr_kmer, copy_len);
                aln->ref_kmer[copy_len] = '\0';
                strncpy(aln->model_kmer, sequence + curr_kmer, copy_len);
                aln->model_kmer[copy_len] = '\0';
            }
            else
            {
                aln->ref_kmer[0] = '\0';
                aln->model_kmer[0] = '\0';
            }

            // Event statistics
            if (curr_event >= 0 && curr_event < events.n)
            {
                aln->event_mean = events.event[curr_event].mean;
                aln->event_stdv = events.event[curr_event].stdv;
                aln->event_duration = events.event[curr_event].length;
                aln->signal_start = events.event[curr_event].start;
                aln->signal_length = events.event[curr_event].length;
            }
            else
            {
                aln->event_mean = 0.0f;
                aln->event_stdv = 0.0f;
                aln->event_duration = 0.0f;
                aln->signal_start = 0;
                aln->signal_length = 0.0f;
            }

            // Model statistics
            if (curr_kmer + kmer_size <= strlen(sequence))
            {
                uint32_t kmer_rank = get_kmer_rank(sequence + curr_kmer, kmer_size);
                aln->model_mean = models[kmer_rank].level_mean;
                aln->model_stdv = models[kmer_rank].level_stdv;
                aln->scaled_model_mean = scaling.scale * aln->model_mean + scaling.shift;
                aln->scaled_model_stdv = aln->model_stdv * scaling.var;
            }
            else
            {
                aln->model_mean = 0.0f;
                aln->model_stdv = 0.0f;
                aln->scaled_model_mean = 0.0f;
                aln->scaled_model_stdv = 0.0f;
            }
        }
        else if (curr_state == STATE_BAD_EVENT && curr_event >= 0 && curr_event < events.n)
        {
            // Record bad events too
            event_alignment_t *aln = &alignment[align_idx++];
            aln->ref_position = curr_kmer;
            aln->event_idx = curr_event;
            aln->hmm_state = 'B';
            aln->strand_idx = 0;

            int copy_len = (curr_kmer + kmer_size <= strlen(sequence)) ? kmer_size : strlen(sequence) - curr_kmer;
            if (copy_len > 0 && copy_len <= MAX_KMER_SIZE)
            {
                strncpy(aln->ref_kmer, sequence + curr_kmer, copy_len);
                aln->ref_kmer[copy_len] = '\0';
                strncpy(aln->model_kmer, sequence + curr_kmer, copy_len);
                aln->model_kmer[copy_len] = '\0';
            }
            else
            {
                aln->ref_kmer[0] = '\0';
                aln->model_kmer[0] = '\0';
            }

            aln->event_mean = events.event[curr_event].mean;
            aln->event_stdv = events.event[curr_event].stdv;
            aln->event_duration = events.event[curr_event].length;
            aln->signal_start = events.event[curr_event].start;
            aln->signal_length = events.event[curr_event].length;

            if (curr_kmer + kmer_size <= strlen(sequence))
            {
                uint32_t kmer_rank = get_kmer_rank(sequence + curr_kmer, kmer_size);
                aln->model_mean = models[kmer_rank].level_mean;
                aln->model_stdv = models[kmer_rank].level_stdv;
                aln->scaled_model_mean = scaling.scale * aln->model_mean + scaling.shift;
                aln->scaled_model_stdv = aln->model_stdv * scaling.var;
            }
            else
            {
                aln->model_mean = 0.0f;
                aln->model_stdv = 0.0f;
                aln->scaled_model_mean = 0.0f;
                aln->scaled_model_stdv = 0.0f;
            }
        }
        else if (curr_state == STATE_KMER_SKIP)
        {
            // Skip recording kmer skips for now (causes issues)
            // Just move to previous state
        }

        // Move to previous state based on movement type
        if (movement == HMT_FROM_SAME_M || movement == HMT_FROM_SAME_B)
        {
            // Stay at same kmer, previous event
            curr_event--;
            curr_state = (movement == HMT_FROM_SAME_M) ? STATE_MATCH : STATE_BAD_EVENT;
        }
        else if (movement == HMT_FROM_PREV_M || movement == HMT_FROM_PREV_B || movement == HMT_FROM_PREV_K)
        {
            // Move to previous kmer
            if (curr_state != STATE_KMER_SKIP)
                curr_event--; // Kmer skip doesn't consume event
            curr_kmer--;

            if (movement == HMT_FROM_PREV_M)
                curr_state = STATE_MATCH;
            else if (movement == HMT_FROM_PREV_B)
                curr_state = STATE_BAD_EVENT;
            else
                curr_state = STATE_KMER_SKIP;
        }
        else if (movement == HMT_FROM_SOFT)
        {
            // Came from soft-clip state - this is the start of the alignment
            // All events before this were soft-clipped (skipped)
            break;
        }
        else
        {
            break; // Unknown movement or start
        }
    }

    // Reverse alignment (we traced backwards)
    for (int i = 0; i < align_idx / 2; ++i)
    {
        event_alignment_t temp = alignment[i];
        alignment[i] = alignment[align_idx - 1 - i];
        alignment[align_idx - 1 - i] = temp;
    }

    *out_alignment = alignment;
    return align_idx;
}

// Main Profile HMM alignment function
int32_t profile_hmm_align(
    event_alignment_t **out_alignment,
    const char *sequence,
    int32_t seq_len,
    event_table events,
    simple_model_t *model,
    uint32_t kmer_size,
    simple_scalings_t scaling,
    float events_per_base)
{
    int32_t n_kmers = seq_len - kmer_size + 1;
    int32_t n_events = events.n;

    if (n_kmers <= 0 || n_events <= 0)
        return 0;

    // Calculate transition probabilities
    BlockTransitions *transitions = calculate_transitions(n_kmers, events_per_base);
    if (!transitions)
        return -1;

    // Allocate DP and path matrices
    int n_states = n_kmers * NUM_STATES;
    float **dp_matrix = (float **)malloc((n_events + 1) * sizeof(float *));
    uint8_t **path_matrix = (uint8_t **)malloc((n_events + 1) * sizeof(uint8_t *));

    if (!dp_matrix || !path_matrix)
    {
        free(transitions);
        free(dp_matrix);
        free(path_matrix);
        return -1;
    }

    for (int i = 0; i <= n_events; ++i)
    {
        dp_matrix[i] = (float *)malloc(n_states * sizeof(float));
        path_matrix[i] = (uint8_t *)malloc(n_states * sizeof(uint8_t));

        if (!dp_matrix[i] || !path_matrix[i])
        {
            for (int j = 0; j <= i; ++j)
            {
                free(dp_matrix[j]);
                free(path_matrix[j]);
            }
            free(dp_matrix);
            free(path_matrix);
            free(transitions);
            return -1;
        }
    }

    // Create pre-flanking and post-flanking arrays for soft-clipping
    float *pre_flank = make_pre_flanking(n_events);
    float *post_flank = make_post_flanking(n_events);

    if (!pre_flank || !post_flank)
    {
        for (int i = 0; i <= n_events; ++i)
        {
            free(dp_matrix[i]);
            free(path_matrix[i]);
        }
        free(dp_matrix);
        free(path_matrix);
        free(transitions);
        free(pre_flank);
        free(post_flank);
        return -1;
    }

    // Run Viterbi HMM with soft-clipping
    profile_hmm_fill_generic_r9(dp_matrix, path_matrix, sequence, n_kmers,
                                events, n_events, model, kmer_size,
                                scaling, transitions, pre_flank, post_flank);

    // Traceback to get alignment (using post_flank for soft-clipping at end)
    int32_t n_aligned = profile_hmm_traceback_with_flanking(out_alignment, path_matrix, dp_matrix,
                                                            sequence, n_kmers, events, n_events,
                                                            model, kmer_size, scaling, post_flank);

    // Cleanup
    for (int i = 0; i <= n_events; ++i)
    {
        free(dp_matrix[i]);
        free(path_matrix[i]);
    }
    free(dp_matrix);
    free(path_matrix);
    free(transitions);
    free(pre_flank);
    free(post_flank);

    return n_aligned;
}
