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

    free(db->read_id);
    free(db->read_len);
    free(db->read_idx);
    free(db->ref_n);
    free(db->ref_sequence);
    free(db->ref_name);
    free(db->ref_len);
    free(db->et);
    free(db->sig);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);
    for (i = 0; i < db->batch_size; ++i)
    {
        delete db->site_score_map[i];
    }
    free(db->site_score_map);
    // eventalign related
    if (db->eventalign_summary)
    {
        free(db->eventalign_summary);
    }
    if (db->event_alignment_result)
    {
        for (i = 0; i < db->batch_size; ++i)
        {
            delete db->event_alignment_result[i];
        }
        free(db->event_alignment_result);
    }
    if (db->event_alignment_result_str)
    {
        free(db->event_alignment_result_str);
    }

    free(db);
}
