/* @file f5c.c
**
** f5c interface implementation
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

#include <sys/wait.h>
#include <unistd.h>

/* initialise a data batch */
db_t *init_db(int32_t batch_size)
{
    db_t *db = (db_t *)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->batch_size = batch_size;
    db->read_idx = 0;

    db->read_id = (char **)(malloc(sizeof(char *) * db->batch_size));
    MALLOC_CHK(db->read_id);
    db->read_len = (int32_t *)(malloc(sizeof(int32_t) * db->batch_size));
    MALLOC_CHK(db->read_len);

    db->sig = (signal_t **)malloc(sizeof(signal_t *) * db->batch_size);
    MALLOC_CHK(db->sig);

    db->et = (event_table *)malloc(sizeof(event_table) * db->batch_size);
    MALLOC_CHK(db->et);

    db->scalings =
        (scalings_t *)malloc(sizeof(scalings_t) * db->batch_size);
    MALLOC_CHK(db->scalings);

    db->event_align_pairs =
        (AlignedPair **)malloc(sizeof(AlignedPair *) * db->batch_size);
    MALLOC_CHK(db->event_align_pairs);
    db->n_event_align_pairs =
        (int32_t *)malloc(sizeof(int32_t) * db->batch_size);
    MALLOC_CHK(db->n_event_align_pairs);

    db->event_alignment = (event_alignment_t **)malloc(
        sizeof(event_alignment_t *) * db->batch_size);
    MALLOC_CHK(db->event_alignment);
    db->n_event_alignment =
        (int32_t *)malloc(sizeof(int32_t) * db->batch_size);
    MALLOC_CHK(db->n_event_alignment);

    db->events_per_base =
        (double *)malloc(sizeof(double) * db->batch_size);
    MALLOC_CHK(db->events_per_base);

    db->base_to_event_map =
        (index_pair_t **)malloc(sizeof(index_pair_t *) * db->batch_size);
    MALLOC_CHK(db->base_to_event_map);

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->batch_size);
    MALLOC_CHK(db->read_stat_flag);

    db->site_score_map = (std::map<int, ScoredSite> **)malloc(sizeof(std::map<int, ScoredSite> *) * db->batch_size);
    MALLOC_CHK(db->site_score_map);

    for (i = 0; i < db->batch_size; ++i)
    {
        db->site_score_map[i] = new std::map<int, ScoredSite>;
        NULL_CHK(db->site_score_map[i]);
    }

    db->total_reads = 0;

    // eventalign related

    db->eventalign_summary = (EventalignSummary *)malloc(sizeof(EventalignSummary) * db->batch_size);
    MALLOC_CHK(db->eventalign_summary);

    db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * db->batch_size);
    MALLOC_CHK(db->event_alignment_result);

    db->event_alignment_result_str = (char **)malloc(sizeof(char *) * db->batch_size);
    MALLOC_CHK(db->event_alignment_result_str);

    for (i = 0; i < db->batch_size; ++i)
    {
        db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        NULL_CHK(db->event_alignment_result[i]);
        (db->eventalign_summary[i]).num_events = 0; // done here in the same loop for efficiency
        db->event_alignment_result_str[i] = NULL;
    }

    return db;
}

/* completely free a data batch */
void free_db(db_t *db)
{
    if (db == NULL) {
        return;
    }

    // Free deep-copied strings
    if (db->read_id) {
        for (int32_t i = 0; i < db->batch_size; i++) {
            if (db->read_id[i]) {
                free(db->read_id[i]);
            }
        }
        free(db->read_id);
    }

    if (db->read_len) {
        free(db->read_len);
    }

    if (db->ref_sequence) {
        for (int32_t i = 0; i < db->batch_size; i++) {
            if (db->ref_sequence[i]) {
                free(db->ref_sequence[i]);
            }
        }
        free(db->ref_sequence);
    }

    if (db->ref_name) {
        for (int32_t i = 0; i < db->batch_size; i++) {
            if (db->ref_name[i]) {
                free(db->ref_name[i]);
            }
        }
        free(db->ref_name);
    }

    if (db->ref_len) {
        free(db->ref_len);
    }

    if (db->et) {
        // Free event data in each event_table
        for (int32_t i = 0; i < db->batch_size; i++) {
            if (db->et[i].event) {
                free(db->et[i].event);
            }
        }
        free(db->et);
    }

    // Free signal data
    if (db->sig) {
        for (int32_t i = 0; i < db->batch_size; i++) {
            if (db->sig[i]) {
                if (db->sig[i]->rawptr) {
                    free(db->sig[i]->rawptr);
                }
                free(db->sig[i]);
            }
        }
        free(db->sig);
    }

    if (db->scalings) {
        free(db->scalings);
    }

    if (db->event_align_pairs) {
        free(db->event_align_pairs);
    }

    if (db->n_event_align_pairs) {
        free(db->n_event_align_pairs);
    }

    if (db->event_alignment) {
        free(db->event_alignment);
    }

    if (db->n_event_alignment) {
        free(db->n_event_alignment);
    }

    if (db->events_per_base) {
        free(db->events_per_base);
    }

    if (db->base_to_event_map) {
        free(db->base_to_event_map);
    }

    if (db->read_stat_flag) {
        free(db->read_stat_flag);
    }

    if (db->site_score_map) {
        for (int32_t i = 0; i < db->batch_size; ++i) {
            delete db->site_score_map[i];
        }
        free(db->site_score_map);
    }

    // eventalign related
    if (db->eventalign_summary) {
        free(db->eventalign_summary);
    }

    if (db->event_alignment_result) {
        for (int32_t i = 0; i < db->batch_size; ++i) {
            delete db->event_alignment_result[i];
        }
        free(db->event_alignment_result);
    }

    if (db->event_alignment_result_str) {
        free(db->event_alignment_result_str);
    }

    free(db);
}

/* =============================================================================
 * init_db_from_python - Initialize db_t from Python-provided data
 *
 * This function creates a db_t structure populated with data provided from Python,
 * rather than reading from BAM/SLOW5 files.
 *
 * Args:
 *   batch_size: Number of reads in the batch
 *   read_ids: Array of read identifier strings (borrowed, not copied)
 *   read_ids_len: Array of read ID string lengths
 *   read_seqs: Array of read sequence strings (borrowed, not copied)
 *   read_lens: Array of read sequence lengths
 *   ref_seqs: Array of reference sequence strings (borrowed, not copied)
 *   ref_seqs_len: Array of reference sequence string lengths
 *   ref_names: Array of reference name strings (borrowed, not copied)
 *   ref_names_len: Array of reference name string lengths
 *   ref_lens: Array of reference sequence lengths
 *   ref_n: Number of reference sequences
 *   signals: Array of pointers to signal data (borrowed, not copied)
 *   signal_lens: Array of signal lengths
 *   signal_drifts: Array of signal drift values (optional, can be NULL)
 *   signal_scales: Array of signal scale values (optional, can be NULL)
 *   signal_shifts: Array of signal shift values (optional, can be NULL)
 *
 * Returns:
 *   Pointer to allocated db_t structure, or NULL on error
 * =============================================================================
 */
db_t *init_db_from_python(
    int32_t batch_size,
    char **read_ids,
    int32_t *read_ids_len,
    char **read_seqs,
    int32_t *read_lens,
    char **ref_seqs,
    int32_t *ref_seqs_len,
    char **ref_names,
    int32_t *ref_names_len,
    int32_t *ref_lens,
    int32_t ref_n,
    float **signals,
    uint64_t *signal_lens,
    float *signal_drifts,
    float *signal_scales,
    float *signal_shifts
)
{
    // Validate inputs
    if (batch_size <= 0) {
        fprintf(stderr, "Error: batch_size must be > 0\n");
        return NULL;
    }

    // Allocate db_t structure
    db_t *db = (db_t *)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->batch_size = batch_size;
    db->read_idx = 0;  // Starting index for this batch
    db->ref_n = ref_n;  // Number of reference sequences

    // Allocate and assign read_id array (deep copy strings from Python)
    db->read_id = (char **)(malloc(sizeof(char *) * batch_size));
    MALLOC_CHK(db->read_id);
    for (int32_t i = 0; i < batch_size; i++) {
        db->read_id[i] = strdup(read_ids[i]);  // Deep copy
        MALLOC_CHK(db->read_id[i]);
    }

    // Allocate and assign read_len array (copied values)
    db->read_len = (int32_t *)(malloc(sizeof(int32_t) * batch_size));
    MALLOC_CHK(db->read_len);
    memcpy(db->read_len, read_lens, sizeof(int32_t) * batch_size);

    // Allocate and assign ref_sequence array (deep copy strings from Python)
    db->ref_sequence = (char **)(malloc(sizeof(char *) * batch_size));
    MALLOC_CHK(db->ref_sequence);
    for (int32_t i = 0; i < batch_size; i++) {
        db->ref_sequence[i] = strdup(ref_seqs[i]);  // Deep copy
        MALLOC_CHK(db->ref_sequence[i]);
    }

    // Allocate and assign ref_name array (deep copy strings from Python)
    db->ref_name = (char **)(malloc(sizeof(char *) * batch_size));
    MALLOC_CHK(db->ref_name);
    for (int32_t i = 0; i < batch_size; i++) {
        db->ref_name[i] = strdup(ref_names[i]);  // Deep copy
        MALLOC_CHK(db->ref_name[i]);
    }

    // Allocate and assign ref_len array (copied values)
    db->ref_len = (int32_t *)(malloc(sizeof(int32_t) * batch_size));
    MALLOC_CHK(db->ref_len);
    memcpy(db->ref_len, ref_lens, sizeof(int32_t) * batch_size);

    // Allocate and initialize signal_t structures (deep copy signal data)
    db->sig = (signal_t **)malloc(sizeof(signal_t *) * batch_size);
    MALLOC_CHK(db->sig);
    for (int32_t i = 0; i < batch_size; i++) {
        db->sig[i] = (signal_t *)malloc(sizeof(signal_t));
        MALLOC_CHK(db->sig[i]);

        // Deep copy the signal data
        db->sig[i]->rawptr = (float *)malloc(signal_lens[i] * sizeof(float));
        MALLOC_CHK(db->sig[i]->rawptr);
        memcpy(db->sig[i]->rawptr, signals[i], signal_lens[i] * sizeof(float));
        db->sig[i]->nsample = signal_lens[i];      // Copied value
        db->sig[i]->digitisation = 8192.0f;
        db->sig[i]->offset = 0.0f;
        db->sig[i]->range = 1.0f;
        db->sig[i]->sample_rate = 4000.0f;

        // Use provided scaling values or defaults
        db->sig[i]->scale = (signal_scales != NULL) ? signal_scales[i] : 1.0f;
        db->sig[i]->shift = (signal_shifts != NULL) ? signal_shifts[i] : 0.0f;
        db->sig[i]->drift = (signal_drifts != NULL) ? signal_drifts[i] : 0.0f;
        db->sig[i]->var = 1.0f;  // Default variance
        db->sig[i]->scale_sd = 0.0f;
        db->sig[i]->var_sd = 0.0f;

#ifdef CACHED_LOG
        db->sig[i]->log_var = log(db->sig[i]->var);
        db->sig[i]->scaled_var = db->sig[i]->var;
        db->sig[i]->log_scaled_var = db->sig[i]->log_var;
#endif
    }

    // Allocate and initialize event_table array
    db->et = (event_table *)malloc(sizeof(event_table) * batch_size);
    MALLOC_CHK(db->et);
    for (int32_t i = 0; i < batch_size; i++) {
        db->et[i].n = 0;
        db->et[i].start = 0;
        db->et[i].end = 0;
        db->et[i].event = NULL;  // Will be filled by getevents
    }

    // Allocate and initialize scalings array
    db->scalings = (scalings_t *)malloc(sizeof(scalings_t) * batch_size);
    MALLOC_CHK(db->scalings);
    for (int32_t i = 0; i < batch_size; i++) {
        db->scalings[i].scale = (signal_scales != NULL) ? signal_scales[i] : 1.0f;
        db->scalings[i].shift = (signal_shifts != NULL) ? signal_shifts[i] : 0.0f;
        db->scalings[i].var = 1.0f;
        db->scalings[i].drift = (signal_drifts != NULL) ? signal_drifts[i] : 0.0f;
        db->scalings[i].scale_sd = 0.0f;
        db->scalings[i].var_sd = 0.0f;

#ifdef CACHED_LOG
        db->scalings[i].log_var = log(db->scalings[i].var);
#endif
    }

    // Allocate arrays for alignment data (will be filled later)
    db->event_align_pairs = (AlignedPair **)malloc(sizeof(AlignedPair *) * batch_size);
    MALLOC_CHK(db->event_align_pairs);
    for (int32_t i = 0; i < batch_size; i++) {
        db->event_align_pairs[i] = NULL;
    }

    db->n_event_align_pairs = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    MALLOC_CHK(db->n_event_align_pairs);
    memset(db->n_event_align_pairs, 0, sizeof(int32_t) * batch_size);

    db->event_alignment = (event_alignment_t **)malloc(sizeof(event_alignment_t *) * batch_size);
    MALLOC_CHK(db->event_alignment);
    for (int32_t i = 0; i < batch_size; i++) {
        db->event_alignment[i] = NULL;
    }

    db->n_event_alignment = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    MALLOC_CHK(db->n_event_alignment);
    memset(db->n_event_alignment, 0, sizeof(int32_t) * batch_size);

    db->events_per_base = (double *)malloc(sizeof(double) * batch_size);
    MALLOC_CHK(db->events_per_base);
    memset(db->events_per_base, 0, sizeof(double) * batch_size);

    db->base_to_event_map = (index_pair_t **)malloc(sizeof(index_pair_t *) * batch_size);
    MALLOC_CHK(db->base_to_event_map);
    for (int32_t i = 0; i < batch_size; i++) {
        db->base_to_event_map[i] = NULL;
    }

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    MALLOC_CHK(db->read_stat_flag);
    memset(db->read_stat_flag, 0, sizeof(int32_t) * batch_size);

    db->site_score_map = (std::map<int, ScoredSite> **)malloc(sizeof(std::map<int, ScoredSite> *) * batch_size);
    MALLOC_CHK(db->site_score_map);
    for (int32_t i = 0; i < batch_size; i++) {
        db->site_score_map[i] = new std::map<int, ScoredSite>;
        NULL_CHK(db->site_score_map[i]);
    }

    // Initialize statistics
    db->sum_bases = 0;
    for (int32_t i = 0; i < batch_size; i++) {
        db->sum_bases += read_lens[i];
    }
    db->total_reads = batch_size;

    // Eventalign related
    db->eventalign_summary = (EventalignSummary *)malloc(sizeof(EventalignSummary) * batch_size);
    MALLOC_CHK(db->eventalign_summary);
    memset(db->eventalign_summary, 0, sizeof(EventalignSummary) * batch_size);

    db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * batch_size);
    MALLOC_CHK(db->event_alignment_result);
    for (int32_t i = 0; i < batch_size; i++) {
        db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        NULL_CHK(db->event_alignment_result[i]);
    }

    db->event_alignment_result_str = (char **)malloc(sizeof(char *) * batch_size);
    MALLOC_CHK(db->event_alignment_result_str);
    memset(db->event_alignment_result_str, 0, sizeof(char *) * batch_size);

    return db;
}
