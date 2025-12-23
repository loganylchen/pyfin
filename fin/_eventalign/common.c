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

ret_status_t load_db(db_t *db,
                     char **read_id, signal_t **sig, int32_t *read_len, int32_t read_num,
                     char **ref_sequence, char **ref_name, int32_t *ref_len, int32_t ref_num)
{

    int32_t result = 0;
    db->read_idx = 0;
    db->sum_bases = 0;
    db->total_reads = 0;

    ret_status_t status = {0, 0};
    int32_t i = 0;

    while (db->read_idx < db->batch_size)
    {
        i = db->read_idx;
        db->read_idx[i] = core->read_index;

        if (result < 0)
        {
            break;
        }
        else
        {
            if ((record->core.flag & BAM_FUNMAP) == 0 && record->core.qual >= core->opt.min_mapq)
            {
                // printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);

                if (!(core->opt.flag & F5C_SECONDARY_YES))
                {
                    if ((record->core.flag & BAM_FSECONDARY))
                    {
                        db->skip_sec_reads++;
                        continue;
                    }
                }

                db->total_reads++; // candidate read

                std::string qname = bam_get_qname(record);
                std::string fast5_path_str;
                t = realtime();
                // todo : make efficient (redudantly accessed below, can be combined with it?)
                int64_t read_length = core->readbb->get_read_sequence(qname).size();
                if (!(core->opt.flag & F5C_RD_SLOW5))
                { // if not slow5
                    fast5_path_str = core->readbb->get_signal_path(qname);
                }
                core->db_fasta_time += realtime() - t;

                // skipping ultra-long-reads
                if (core->ultra_long_tmp != NULL && read_length > core->opt.ultra_thresh)
                {
                    db->ultra_long_skipped++;
                    int ret_wr = sam_write1(core->ultra_long_tmp, core->m_hdr, record);
                    NEG_CHK(ret_wr);
                    continue;
                }

                if (!(core->opt.flag & F5C_RD_SLOW5) && fast5_path_str == "")
                {
                    handle_bad_fast5(core, db, fast5_path_str, qname);
                    continue;
                }

                int8_t read_status = 0;
                if (core->opt.flag & F5C_RD_RAW_DUMP)
                {
                    t = realtime();
                    read_status = read_from_fast5_dump(core, db, i);
                    double rt = realtime() - t;
                    core->db_fast5_read_time += rt;
                    core->db_fast5_time += rt;
                }
                else if (core->opt.flag & F5C_RD_SLOW5)
                {
                    read_status = 1; // we do this later with multiple threaded
                }
                else
                {
                    read_status = read_from_fast5_files(core, db, qname, fast5_path_str, i);
                }
                if (read_status == 1)
                {
                    db->n_bam_rec++;
                    status.num_bases += read_length;
                }
            }
            else
            {
                if (record->core.qual < core->opt.min_mapq)
                {
                    db->skip_mapq_reads++;
                }
                else
                {
                    db->unmapped_reads++;
                }
            }
        }
    }
    // fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

    // get ref sequences (todo can make efficient by taking the the start and end of the sorted bam)
    for (i = 0; i < db->n_bam_rec; i++)
    {
        bam1_t *record = db->bam_rec[i];
        char *ref_name = core->m_hdr->target_name[record->core.tid];
        // printf("refname : %s\n",ref_name);
        int32_t ref_start_pos = record->core.pos;
        int32_t ref_end_pos = bam_endpos(record);
        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        int32_t fetched_len = 0;
        t = realtime();
        char *refseq = faidx_fetch_seq(core->fai, ref_name, ref_start_pos, ref_end_pos, &fetched_len); // todo : error handle?
        core->db_fasta_time += realtime() - t;
        db->fasta_cache[i] = refseq;
        // printf("seq : %s\n",db->fasta_cache[i]);

        // get the read in ASCII
        std::string qname = bam_get_qname(db->bam_rec[i]);

        t = realtime();
        std::string read_seq = core->readbb->get_read_sequence(qname);
        core->db_fasta_time += realtime() - t;

        db->read[i] = (char *)malloc(read_seq.size() + 1); // todo : is +1 needed? do errorcheck
        strcpy(db->read[i], read_seq.c_str());
        db->read_len[i] = strlen(db->read[i]);
        if (core->opt.flag & F5C_RNA)
        {
            replace_char(db->read[i], 'U', 'T');
        }
        db->sum_bases += db->read_len[i];

        db->read_stat_flag[i] = 0; // reset the flag
    }
    // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);
    if (core->opt.verbosity > 1)
    {
        STDERR("Average read len %.0f", db->sum_bases / (float)db->n_bam_rec);
    }
    status.num_reads = db->n_bam_rec;
    assert(status.num_bases == db->sum_bases);

    // read the slow5 batch
    if (core->opt.flag & F5C_RD_SLOW5)
    {
        t = realtime();
        pthread_db(core, db, read_slow5_single);
        core->db_fast5_time += realtime() - t;
    }

    double load_end = realtime();
    core->load_db_time += (load_end - load_start);

    return status;
}